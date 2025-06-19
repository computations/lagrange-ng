/* Common.h
 *
 * Created On: 27 Oct 2020
 * Author: Ben Bettisworth
 */
#ifndef LAGRANGE_COMMON_H
#define LAGRANGE_COMMON_H

#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>

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

using Range = uint64_t;
using RangeMask = Range;

using ClockTick = uint64_t;
using Clock = std::atomic<uint64_t>;

using OpID = uint64_t;

struct NodeReservation {
  NodeReservation() :
      _top_clv{std::numeric_limits<size_t>::max()},
      _bot1_clv{std::numeric_limits<size_t>::max()},
      _bot2_clv{std::numeric_limits<size_t>::max()},
      _top_rclv{std::numeric_limits<size_t>::max()},
      _bot1_rclv{std::numeric_limits<size_t>::max()},
      _bot2_rclv{std::numeric_limits<size_t>::max()} {}

  size_t _top_clv;
  size_t _bot1_clv;
  size_t _bot2_clv;

  size_t _top_rclv;
  size_t _bot1_rclv;
  size_t _bot2_rclv;
};

struct PeriodDerivative {
  double d_dispersion;
  double d_extinction;

  auto norm() const -> double {
    return d_dispersion * d_dispersion + d_extinction * d_extinction;
  }
};

struct PeriodParams {
  double dispersion_rate;
  double extinction_rate;
  std::shared_ptr<double[]> adjustment_matrix = nullptr;
  size_t regions;

  void applyDerivative(const PeriodDerivative &d) {
    dispersion_rate += d.d_dispersion;
    extinction_rate += d.d_extinction;
    if (dispersion_rate < 0) { dispersion_rate = 0.0; }
    if (extinction_rate < 0) { extinction_rate = 0.0; }
  }

  auto toString() const -> std::string {
    std::ostringstream os;
    os << "(disp: " << dispersion_rate << ", ext: " << extinction_rate << ")";
    return os.str();
  }

  inline auto getDispersionRate(size_t from, size_t to) const -> double {
    return dispersion_rate
           * (adjustment_matrix != nullptr
                  ? adjustment_matrix[from * regions + to]
                  : 1.0);
  }

  inline auto getExtinctionRate() const -> double { return extinction_rate; }

  double operator[](size_t i) const {
    if (i == 0) { return dispersion_rate; }
    if (i == 1) { return extinction_rate; }
    return std::numeric_limits<double>::quiet_NaN();
  }

  double &operator[](size_t i) {
    assert(i < 2);
    if (i == 0) { return dispersion_rate; }
    return extinction_rate;
  }
};

using LagrangeFloat = double;
using LagrangeComplex = lapack_complex_double;

using LagrangeMatrixBase = LagrangeFloat;

using LagrangeMatrix = LagrangeMatrixBase *;
using LagrangeConstMatrix = const LagrangeMatrixBase *const;
using LagrangeComplexMatrix = LagrangeMatrixBase *;
using LagrangeColVector = LagrangeMatrixBase *;
using LagrangeConstColVector = const LagrangeMatrixBase *const;
using LagrangeComplexColVector = LagrangeMatrixBase *;

constexpr double lagrange_scaling_factor = 0x1p256;

constexpr double lagrange_scale_threshold = 1.0 / lagrange_scaling_factor;

constexpr double lagrange_scaling_factor_log =
    177.445678223345993274051579646766185760498046875;
/* = std::log(lagrange_scaling_factor) */

enum class LagrangeOperationMode { OPTIMIZE, EVALUATE };

enum class LagrangeEXPMComputationMode { PADE, KRYLOV, ADAPTIVE };
}  // namespace lagrange
#endif  // LAGRANGE_COMMON_H
