/* Common.h
 *
 * Created On: 27 Oct 2020
 * Author: Ben Bettisworth
 */
#ifndef LAGRANGE_COMMON_H
#define LAGRANGE_COMMON_H

#include <atomic>
#include <cassert>
#include <cstdint>
#include <span>

#ifndef MKL_ENABLED
  #include <complex>

  #define COMPLEX_DATA_TYPE double
  #define lapack_complex_double std::complex<COMPLEX_DATA_TYPE>
  #define lapack_complex_double_real(z) \
    (reinterpret_cast<COMPLEX_DATA_TYPE *>(&z)[0])
  #define lapack_complex_double_imag(z) \
    (reinterpret_cast<COMPLEX_DATA_TYPE *>(&z)[1])
#else
  #define lapack_complex_double_real(z) (z.real)
  #define lapack_complex_double_imag(z) (z.imag)
#endif

#include "Quarantine.hpp"

namespace lagrange {

using LagrangeFloat = double;
using LagrangeComplex = lapack_complex_double;

using LagrangeMatrixBase = LagrangeFloat;
using StateReturnList = std::span<LagrangeMatrixBase>;

using LagrangeMatrix = LagrangeMatrixBase *;
using LagrangeConstMatrix = const LagrangeMatrixBase *const;
using LagrangeComplexMatrix = LagrangeMatrixBase *;
using LagrangeColVector = LagrangeMatrixBase *;
using LagrangeConstColVector = const LagrangeMatrixBase *const;
using LagrangeComplexColVector = LagrangeMatrixBase *;

using Range = uint64_t;
using RangeMask = Range;

using ClockTick = uint64_t;
using Clock = std::atomic<uint64_t>;

using OpID = uint64_t;

constexpr double lagrange_scaling_factor = 0x1p256;

constexpr double lagrange_scale_threshold = 1.0 / lagrange_scaling_factor;

constexpr double lagrange_scaling_factor_log =
    177.445678223345993274051579646766185760498046875;
/* = std::log(lagrange_scaling_factor) */

enum class LagrangeOperationMode { OPTIMIZE, EVALUATE };

enum class LagrangeEXPMComputationMode { PADE, KRYLOV, ADAPTIVE };

enum class LagrangeResultEvaluationMode { IMMEDIATE, STREAMING };
}  // namespace lagrange
#endif  // LAGRANGE_COMMON_H
