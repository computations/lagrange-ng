#ifndef LAGRANGE_QUARANTINE_H
#pragma GCC system_header

#ifdef MKL_ENABLED
#include <mkl.h>
#else
#include "openblas/cblas.h"
#include "openblas/lapacke.h"
#endif

#endif
