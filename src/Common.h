/* Common.h
 *
 * Created On: 27 Oct 2020
 * Author: Ben Bettisworth
 */
#ifndef LAGRANGE_COMMON_H
#define LAGRANGE_COMMON_H

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

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

#include "Quarantine.h"

using lagrange_dist_t = uint64_t;

using lagrange_clock_tick_t = uint64_t;
using lagrange_clock_t = std::atomic<uint64_t>;

using lagrange_op_id_t = uint64_t;

struct node_reservation_t {
  node_reservation_t()
      : _top_clv{std::numeric_limits<size_t>::max()},
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

struct period_derivative_t {
  double d_dispersion;
  double d_extinction;

  auto norm() const -> double {
    return d_dispersion * d_dispersion + d_extinction * d_extinction;
  }
};

struct period_t {
  double dispersion_rate;
  double extinction_rate;
  std::shared_ptr<std::vector<std::vector<double>>> adjustment_matrix = nullptr;

  void applyDerivative(const period_derivative_t &d) {
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
    return dispersion_rate * (adjustment_matrix != nullptr
                                  ? (*adjustment_matrix)[from][to]
                                  : 1.0);
  }

  inline auto getExtinctionRate() const -> double { return extinction_rate; }
};

using lagrange_float_t = double;
using lagrange_complex_t = lapack_complex_double;

using lagrange_matrix_base_t = lagrange_float_t;

using lagrange_matrix_t = lagrange_matrix_base_t *;
using lagrange_const_matrix_t = const lagrange_matrix_base_t *const;
using lagrange_complex_matrix_t = lagrange_matrix_base_t *;
using lagrange_col_vector_t = lagrange_matrix_base_t *;
using lagrange_const_col_vector_t = const lagrange_matrix_base_t *const;
using lagrange_complex_col_vector_t = lagrange_matrix_base_t *;

constexpr double lagrange_scaling_factor = 0x1p256;

constexpr double lagrange_scale_threshold = 1.0 / lagrange_scaling_factor;

constexpr double lagrange_scaling_factor_log =
    177.445678223345993274051579646766185760498046875;
/* = std::log(lagrange_scaling_factor) */

enum class lagrange_operation_mode { OPTIMIZE, EVALUATE };

enum class lagrange_expm_computation_mode { PADE, KRYLOV };

#endif  // LAGRANGE_COMMON_H
