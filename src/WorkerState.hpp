#ifndef LAGRANGE_WORKER_STATE_H
#define LAGRANGE_WORKER_STATE_H

#include <omp.h>
#include <pthread.h>

#include <atomic>
#include <logger.hpp>
#include <memory>
#include <mutex>
#include <vector>

#include "Goal.hpp"
#include "Operation.hpp"
#include "Workspace.hpp"

namespace lagrange {

/**
 * Barrier class that wraps either POSIX barriers (when available)
 * or a custom implementation for platforms like macOS that lack them.
 */
class Barrier {
 public:
  explicit Barrier(size_t count);
  ~Barrier();

  // Non-copyable, non-movable
  Barrier(const Barrier&) = delete;
  Barrier& operator=(const Barrier&) = delete;
  Barrier(Barrier&&) = delete;
  Barrier& operator=(Barrier&&) = delete;

  void wait();

 private:
#if !defined _POSIX_BARRIERS || _POSIX_BARRIERS < 0 \
    || defined LAGRANGE_USER_BARRIER
  // Custom implementation for macOS and other platforms without POSIX barriers
  size_t _needed;
  size_t _called;
  pthread_mutex_t _mutex;
  pthread_cond_t _cond;
#else
  // POSIX barrier for systems that support it
  pthread_barrier_t _barrier;
#endif
};

enum class WorkerMode : uint8_t {
  ComputeForward,
  ComputeReverse,
  ComputeLH,
  ComputeStateGoal,
  ComputeSplitGoal,
  StreamSplitGoal,
  Halt,
  Unknown
};

struct WorkerContext {
 public:
  WorkerContext(std::vector<std::shared_ptr<SplitOperation>>& fwb,
                std::vector<std::shared_ptr<ReverseSplitOperation>>& rwb,
                std::vector<LLHGoal>& lhg,
                std::vector<StateLHGoal>& stg,
                std::vector<SplitLHGoal>& slg) :
      _forward_work_buffer{fwb},
      _reverse_work_buffer{rwb},
      _lh_goal{lhg},
      _state_lh_work_buffer{stg},
      _split_lh_work_buffer{slg} {}

  WorkerContext(const WorkerContext&) = delete;
  WorkerContext& operator=(const WorkerContext&) = delete;

  auto activeThreads() -> size_t { return _total_threads - _finished_threads; }

  auto acquireLock() -> std::unique_lock<std::mutex> {
    return std::unique_lock<std::mutex>{_context_lock};
  }

  auto acquireLockDefer() -> std::unique_lock<std::mutex> {
    return std::unique_lock<std::mutex>{_context_lock, std::defer_lock};
  }

  void setTotalThreads(size_t threads) {
    auto lock = acquireLock();
    _total_threads = threads;
  }

  void lockAndAddFinishedThread() {
    auto lock = acquireLock();
    _finished_threads++;
  }

  void setRunMode(WorkerMode wm) {
    auto lock = acquireLock();
    _mode = wm;
  }

  void initBarrier() {
    auto lock = acquireLock();
    _barrier = std::make_unique<Barrier>(_total_threads);
  }

  void startIter(WorkerMode wm) {
    auto lock = acquireLock();
    _mode = wm;
    _start_index = 0;
    _finished_threads = 0;
  }

  WorkerMode getMode() {
    auto lock = acquireLock();
    return _mode;
  }

  void barrier() { _barrier->wait(); }

  std::vector<std::shared_ptr<SplitOperation>>& _forward_work_buffer;
  std::vector<std::shared_ptr<ReverseSplitOperation>>& _reverse_work_buffer;
  std::vector<LLHGoal>& _lh_goal;
  std::vector<StateLHGoal>& _state_lh_work_buffer;
  std::vector<SplitLHGoal>& _split_lh_work_buffer;
  size_t _start_index;
  size_t _finished_threads;
  size_t _total_threads;
  std::atomic<WorkerMode> _mode;
  std::mutex _context_lock;
  std::unique_ptr<Barrier> _barrier;
};

class WorkerState {
 public:
  WorkerState(size_t i) : _tid{i} {}

  WorkerState(const WorkerState&) = delete;
  auto operator=(const WorkerState&) -> WorkerState& = delete;

  WorkerState(WorkerState&& ts) noexcept { *this = std::move(ts); }

  auto operator=(WorkerState&& ts) noexcept -> WorkerState& {
    _tid = ts._tid;
    _assigned_threads = ts._assigned_threads;
    return *this;
  }

  void work(WorkerContext& tc, const std::shared_ptr<Workspace>& ws) {
#ifdef MKL_ENABLED
    mkl_set_num_threads(_assigned_threads);
#else
    openblas_set_num_threads(_assigned_threads);
#endif
    omp_set_num_threads(_assigned_threads);
    while (true) {
      tc.barrier();

      WorkerMode mode = WorkerMode::Unknown;
      {
        auto lock = tc.acquireLock();
        mode = tc._mode.load();
      }

      switch (mode) {
        case WorkerMode::ComputeForward:
          LOG_DEBUG("Worker {} computing forward operations", _tid);
          work(tc, tc._forward_work_buffer, ws);
          break;

        case WorkerMode::ComputeReverse:
          LOG_DEBUG("Worker {} computing reverse operations", _tid);
          work(tc, tc._reverse_work_buffer, ws);
          break;

        case WorkerMode::ComputeSplitGoal:
          LOG_DEBUG("Worker {} computing split operations", _tid);
          workGoal(tc, tc._split_lh_work_buffer, ws);
          break;

        case WorkerMode::ComputeStateGoal:
          LOG_DEBUG("Worker {} computing state operations", _tid);
          workGoal(tc, tc._state_lh_work_buffer, ws);
          break;

        case WorkerMode::ComputeLH:
          if (masterThread()) { workGoal(tc, tc._lh_goal, ws); }
          break;

        case WorkerMode::Halt:
          LOG_DEBUG("Worker {} halteng", _tid);
          return;

        default:
          LOG(ERROR, "Found a mode that doesn't exist");
          throw std::runtime_error{"Found a mode that doesn't exist"};
      }
      tc.lockAndAddFinishedThread();
      tc.barrier();
      if (masterThread()) { return; }
    }
  }

  void work(WorkerMode tm,
            WorkerContext& tc,
            const std::shared_ptr<Workspace>& ws) {
    tc.startIter(tm);
    work(tc, ws);
  }

  [[nodiscard]] auto masterThread() const -> bool { return _tid == 0; }

  [[nodiscard]] auto threadID() const -> size_t { return _tid; }

  void setAssignedThreads(size_t at) { _assigned_threads = at; }

  void assign_threads() const {
#ifdef MKL_ENABLED
    mkl_set_num_threads(_assigned_threads);
#else
    openblas_set_num_threads(_assigned_threads);
#endif
  }

 private:
  template <typename T>
  void workGoal(WorkerContext& wc,
                std::vector<T>& work_buffer,
                const std::shared_ptr<Workspace>& workspace) {
    if (!masterThread()) { return; }

    auto lock = wc.acquireLock();

    for (auto& w : work_buffer) {
      while (!w.ready(workspace)) {}
      w.eval(workspace);
    }
  }

  template <typename T>
  void work(WorkerContext& wc,
            std::vector<std::shared_ptr<T>>& work_buffer,
            const std::shared_ptr<Workspace>& workspace) {
    for (auto w = findWork(wc, work_buffer, workspace); w != nullptr;
         w = findWork(wc, work_buffer, workspace)) {
      /* We only lock on greater than 2 here because we have 2 copies of the
       * shared pointer here. One in the vector, and one returned by value from
       * the find_work function */
      std::unique_lock<std::mutex> lock(w->getLock(), std::defer_lock);
      if (lock.try_lock()) { w->eval(workspace); }
    }
  }

  template <typename T>
  auto findWorkGoal(WorkerContext& tc,
                    std::vector<T>& work_buffer,
                    const std::shared_ptr<Workspace>& workspace) ->
      typename std::vector<T>::iterator {
    auto lock = tc.acquireLock();

    if (work_buffer.size() - tc._start_index == 0
        || tc.activeThreads() > work_buffer.size()) {
      return work_buffer.end();
    }

    size_t local_index = tc._start_index++;
    while (!work_buffer[local_index].ready(workspace)) {}
    return work_buffer.begin() + local_index;
  }

  template <typename T>
  auto findWork(WorkerContext& tc,
                std::vector<std::shared_ptr<T>>& work_buffer,
                const std::shared_ptr<Workspace>& workspace)
      -> std::shared_ptr<T> {
    assert(!work_buffer.empty());
    auto lock = tc.acquireLock();

    if (work_buffer.size() - tc._start_index == 0) { return {}; }

    size_t local_index = tc._start_index;
    uint8_t iter = 0;

    /* Spin until we find a ready operation */
    while (iter < 100) {
      if (local_index >= work_buffer.size()) {
        local_index = tc._start_index;
        iter += 1;
      }

      if (work_buffer[local_index]->ready(workspace)) {
        auto op = work_buffer[local_index];
        std::swap(work_buffer[local_index], work_buffer[tc._start_index]);
        tc._start_index++;
        return op;
      }

      local_index++;
    }

    /* This should never be hit, but put it here just in case */
    return {};
  }

  size_t _tid{};

  size_t _assigned_threads{1};
};
}  // namespace lagrange

#endif
