#ifndef LAGRANGE_WORKER_STATE_H
#define LAGRANGE_WORKER_STATE_H

#include <atomic>
#include <cctype>
#include <chrono>
#include <functional>
#include <memory>
#include <mutex>
#include <vector>

#include "Common.h"
#include "Operation.h"
#include "Quarantine.h"
#include "Utils.h"
#include "Workspace.h"

enum class WorkerMode {
  ComputeForward,
  ComputeReverse,
  ComputeLH,
  ComputeStateGoal,
  ComputeSplitGoal,
  Halt
};

struct WorkerContext {
 public:
  WorkerContext(std::vector<std::shared_ptr<SplitOperation>>& fwb,
                std::vector<std::shared_ptr<ReverseSplitOperation>>& rwb,
                std::vector<LLHGoal>& lhg, std::vector<StateLHGoal>& stg,
                std::vector<SplitLHGoal>& slg)
      : _forward_work_buffer{fwb},
        _reverse_work_buffer{rwb},
        _lh_goal{lhg},
        _state_lh_work_buffer{stg},
        _split_lh_work_buffer{slg} {}

  std::vector<std::shared_ptr<SplitOperation>>& _forward_work_buffer;
  std::vector<std::shared_ptr<ReverseSplitOperation>>& _reverse_work_buffer;
  std::vector<LLHGoal>& _lh_goal;
  std::vector<StateLHGoal>& _state_lh_work_buffer;
  std::vector<SplitLHGoal>& _split_lh_work_buffer;
};

class WorkerState {
 public:
  WorkerState() : _tid{_total_threads++} {}

  WorkerState(const WorkerState&) = delete;
  auto operator=(const WorkerState&) -> WorkerState& = delete;

  WorkerState(WorkerState&& ts) noexcept { *this = std::move(ts); }
  auto operator=(WorkerState&& ts) noexcept -> WorkerState& {
    _tid = ts._tid;
    _assigned_threads = ts._assigned_threads;
    ++_total_threads;
    return *this;
  }

  ~WorkerState() {
    if (_total_threads > 0) { _total_threads--; }
  }

  void set_mode(WorkerMode tm) const {
    if (master_thread()) { _mode = tm; }
  }

  void work(WorkerContext& tc, const std::shared_ptr<Workspace>& ws) {
#ifdef MKL_ENABLED
    mkl_set_num_threads(_assigned_threads);
#else
    openblas_set_num_threads(_assigned_threads);
#endif
    while (true) {
      barrier();
      if (master_thread()) {
        _start_index = 0;
        _finished_threads = 0;
      }
      switch (_mode) {
        case WorkerMode::ComputeForward:
          work(tc._forward_work_buffer, ws);
          break;

        case WorkerMode::ComputeReverse:
          work(tc._reverse_work_buffer, ws);
          break;

        case WorkerMode::ComputeSplitGoal:
          work_goal(tc._split_lh_work_buffer, ws);
          break;

        case WorkerMode::ComputeStateGoal:
          work_goal(tc._state_lh_work_buffer, ws);
          break;

        case WorkerMode::ComputeLH:
          if (master_thread()) { work_goal(tc._lh_goal, ws); }
          break;

        case WorkerMode::Halt:
          return;

        default:
          throw std::runtime_error{"Found a mode that doesn't exist"};
      }
      _finished_threads += 1;
      if (master_thread()) { return; }
    }
  }

  inline void work(WorkerMode tm, WorkerContext& tc,
                   const std::shared_ptr<Workspace>& ws) {
    set_mode(tm);
    work(tc, ws);
  }

  auto master_thread() const -> bool { return _tid == 0; }
  auto thread_id() const -> size_t { return _tid; }

  void halt_threads() {
    set_mode(WorkerMode::Halt);
    barrier();
  }

  void set_assigned_threads(size_t at) { _assigned_threads = at; }

  void assign_threads() const {
#ifdef MKL_ENABLED
    mkl_set_num_threads(_assigned_threads);
#else
    openblas_set_num_threads(_assigned_threads);
#endif
  }

 private:
  template <typename T>
  void work_goal(std::vector<T>& work_buffer,
                 const std::shared_ptr<Workspace>& workspace) {
    if (!master_thread()) { return; }

    for (auto& w : work_buffer) {
      while (!w.ready(workspace)) {}
      w.eval(workspace);
    }
  }

  template <typename T>
  void work(std::vector<std::shared_ptr<T>>& work_buffer,
            const std::shared_ptr<Workspace>& workspace) {
    for (auto w = find_work(work_buffer, workspace); w != nullptr;
         w = find_work(work_buffer, workspace)) {
      /* We only lock on greater than 2 here because we have 2 copies of the
       * shared pointer here. One in the vector, and one returned by value from
       * the find_work function
       */
      if (w.use_count() > 2) {
        std::lock_guard<std::mutex> lock(w->getLock());
        w->eval(workspace);
      } else {
        w->eval(workspace);
      }
    }
  }

  void barrier() const {
    static volatile int wait_flag = 0;
    // static std::atomic<size_t> barrier_threads{0};
    static size_t barrier_threads = 0;

    // if (_total_threads == 1) { return; }

    // barrier_threads++;
    __sync_fetch_and_add(&barrier_threads, 1);

    if (master_thread()) {
      while (barrier_threads < _total_threads) {}
      barrier_threads = 0;
      wait_flag = wait_flag == 0 ? 1 : 0;
    } else {
      static thread_local volatile int local_wait_flag = 0;
      while (local_wait_flag == wait_flag) {}
      local_wait_flag = local_wait_flag == 0 ? 1 : 0;
    }
  }

  template <typename T>
  auto find_work_goal(std::vector<T>& work_buffer,
                      const std::shared_ptr<Workspace>& workspace) ->
      typename std::vector<T>::iterator {
    // std::lock_guard<std::mutex> work_lock(_work_buffer_mutex);
    if (work_buffer.size() - _start_index == 0 ||
        active_threads() > work_buffer.size()) {
      return work_buffer.end();
    }

    size_t local_index = _start_index++;
    while (!work_buffer[local_index].ready(workspace)) {}
    return work_buffer.begin() + local_index;
  }

  template <typename T>
  auto find_work(std::vector<std::shared_ptr<T>>& work_buffer,
                 const std::shared_ptr<Workspace>& workspace)
      -> std::shared_ptr<T> {
    assert(!work_buffer.empty());
    // auto t1 = std::chrono::high_resolution_clock::now();
    std::lock_guard<std::mutex> work_lock(_work_buffer_mutex);
    /*
    std::cout << "[thread: " << _tid
              << "] looking for work, starting index: " << _start_index
              << std::endl;
              */

    if (work_buffer.size() - _start_index == 0 ||
        active_threads() > work_buffer.size()) {
      return {};
    }

    size_t local_index = _start_index;

    /* Spin lock until we find a ready operation */
    while (true) {
      if (local_index >= work_buffer.size()) { local_index = _start_index; }

      if (work_buffer[local_index]->ready(workspace)) {
        // if (local_index != _start_index) {
        std::swap(work_buffer[local_index], work_buffer[_start_index]);
        //}
        auto op = work_buffer[_start_index];
        _start_index++;
        /*
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> d = t2 - t1;
        std::cout << "Took " << d.count() << " seconds to find work"
                  << std::endl;
                  */
        return op;
      }

      local_index++;
    }

    /* This should never be hit, but put it here just in case */
    return {};
  }

  static inline auto active_threads() -> size_t {
    return _total_threads - _finished_threads;
  }

  size_t _tid{};

  size_t _assigned_threads{0};

  static std::shared_ptr<SplitOperation> _forward_work_buffer;
  static std::shared_ptr<ReverseSplitOperation> _backwards_work_buffer;
  static std::mutex _work_buffer_mutex;
  static std::mutex _io_lock;
  static std::atomic_size_t _total_threads;
  static std::atomic_size_t _finished_threads;
  static std::atomic_size_t _start_index;
  static std::atomic<WorkerMode> _mode;
};

#endif
