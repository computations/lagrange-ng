#include "ThreadState.h"
std::mutex ThreadState::_io_lock;
std::mutex ThreadState::_work_buffer_mutex;
std::atomic_size_t ThreadState::_total_threads;
std::atomic_size_t ThreadState::_finished_threads;
std::atomic_size_t ThreadState::_start_index;
std::atomic<ThreadMode> ThreadState::_thread_mode;
