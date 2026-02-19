#include "WorkerState.hpp"

namespace lagrange {

#if !defined _POSIX_BARRIERS || _POSIX_BARRIERS < 0 \
    || defined LAGRANGE_USER_BARRIER

Barrier::Barrier(size_t count) : _needed(count), _called(0) {
  pthread_mutex_init(&_mutex, nullptr);
  pthread_cond_init(&_cond, nullptr);
}

Barrier::~Barrier() {
  LOG_ASSERT(!pthread_mutex_destroy(&_mutex));
  LOG_ASSERT(!pthread_cond_destroy(&_cond));
}

void Barrier::wait() {
  pthread_mutex_lock(&_mutex);
  _called++;
  if (_called == _needed) {
    _called = 0;
    pthread_cond_broadcast(&_cond);
  } else {
    pthread_cond_wait(&_cond, &_mutex);
  }
  pthread_mutex_unlock(&_mutex);
}

#else

Barrier::Barrier(size_t count) {
  LOG_ASSERT(pthread_barrier_init(
      &_barrier, nullptr, static_cast<unsigned int>(count)));
}

Barrier::~Barrier() { pthread_barrier_destroy(&_barrier); }

void Barrier::wait() { pthread_barrier_wait(&_barrier); }

#endif

};  // namespace lagrange
