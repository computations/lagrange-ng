#ifndef LAGRANGE_THREAD_STATE_H__
#define LAGRANGE_THREAD_STATE_H__

#include <atomic>
#include <cctype>
#include <functional>
#include <memory>
#include <mutex>

#include "Common.h"
#include "Operation.h"
#include "Utils.h"
#include "Workspace.h"

class ThreadState {
 public:
  ThreadState() : _tid{_total_threads++} {}
  ThreadState(const ThreadState&) = delete;
  ThreadState& operator=(const ThreadState&) = delete;

  ThreadState(ThreadState&& ts) { *this = std::move(ts); }
  ThreadState& operator=(ThreadState&& ts) {
    _tid = ts._tid;
    ++_total_threads;
    return *this;
  }

  ~ThreadState() {
    if (_total_threads > 0) {
      _total_threads--;
    }
  }

  template <typename T>
  void work(std::vector<T>& work_buffer,
            const std::shared_ptr<Workspace>& workspace) {
    start_thread();

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

    end_thread();
  }

 private:
  template <typename T>
  T find_work(std::vector<T>& work_buffer,
              const std::shared_ptr<Workspace>& workspace) {
    std::lock_guard<std::mutex> work_lock(_work_buffer_mutex);
    size_t local_index = _start_index;

    if (work_buffer.size() - _start_index == 0 ||
        _total_threads > work_buffer.size()) {
      return {};
    }

    /* Spin lock until we find a ready operation */
    while (true) {
      if (local_index >= work_buffer.size()) {
        local_index = _start_index;
      }

      if (work_buffer[local_index]->ready(workspace)) {
        std::swap(work_buffer[local_index], work_buffer[_start_index]);
        auto op = work_buffer[_start_index];
        _start_index++;
        return op;
      }

      local_index++;
    }

    /* This should never be hit, but put it here just in case */
    return {};
  }

  void start_thread() {}
  void end_thread() {
    _finished_threads++;
    if (_finished_threads == _total_threads) {
      _start_index = 0;
      _finished_threads = 0;
    }
  }

  size_t _tid;

  static std::mutex _work_buffer_mutex;
  static std::mutex _io_lock;
  static std::atomic_size_t _total_threads;
  static std::atomic_size_t _finished_threads;
  static std::atomic_size_t _start_index;
};

#endif
