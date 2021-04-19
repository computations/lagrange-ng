#include "ThreadState.h"
std::mutex ThreadState::_work_buffer_mutex;
std::atomic_size_t ThreadState::_total_threads;
std::atomic_size_t ThreadState::_start_index;
