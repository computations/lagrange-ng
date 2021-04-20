#ifndef LAGRANGE_THREAD_STATE_H__
#define LAGRANGE_THREAD_STATE_H__

#include <atomic>
#include <functional>
#include <memory>
#include <mutex>

#include "Common.h"
#include "Operation.h"
#include "Utils.h"
#include "Workspace.h"

class ThreadState {
 public:
  template <typename T>
  void work(std::vector<T>& work_buffer,
            const std::shared_ptr<Workspace>& workspace) {
    start_thread();

    for (auto w = find_work(work_buffer, workspace); w != nullptr;
         w = find_work(work_buffer, workspace)) {
      w->eval(workspace);
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
      if (local_index > work_buffer.size()) {
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

  void start_thread() { _total_threads++; }
  void end_thread() {
    _total_threads--;
    if (_total_threads == 0) {
      _start_index = 0;
    }
  }

  size_t _tid;

  static std::mutex _work_buffer_mutex;
  static std::atomic_size_t _total_threads;
  static std::atomic_size_t _start_index;
};

#endif
