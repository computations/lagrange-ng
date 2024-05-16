#include "WorkerState.hpp"

std::mutex WorkerState::_io_lock;
std::mutex WorkerState::_work_buffer_mutex;
std::atomic_size_t WorkerState::_total_threads;
std::atomic_size_t WorkerState::_finished_threads;
std::atomic_size_t WorkerState::_start_index;
std::atomic<WorkerMode> WorkerState::_mode;
