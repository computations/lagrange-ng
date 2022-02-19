#ifndef LAGRANGE_QUARANTINE_H__
#pragma GCC system_header

#ifdef MKL_ENABLED
#include <mkl.h>
#else
#include "cblas.h"
#include "lapack.h"
#endif

#endif
