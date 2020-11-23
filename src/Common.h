/* Common.h
 *
 * Created On: 27 Oct 2020
 * Author: Ben Bettisworth
 */
#ifndef LAGRANGE_COMMON_H__
#define LAGRANGE_COMMON_H__

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

typedef uint64_t lagrange_dist_t;

typedef uint64_t lagrange_clock_t;

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

  double norm() const {
    return d_dispersion * d_dispersion + d_extinction * d_extinction;
  }
};

struct period_t {
  double dispersion_rate;
  double extinction_rate;
  std::shared_ptr<std::vector<std::vector<double>>> adjustment_matrix = nullptr;

  void applyDerivative(const period_derivative_t& d) {
    dispersion_rate += d.d_dispersion;
    extinction_rate += d.d_extinction;
    if (dispersion_rate < 0) {
      dispersion_rate = 0.0;
    }
    if (extinction_rate < 0) {
      extinction_rate = 0.0;
    }
  }

  std::string toString() const {
    std::ostringstream os;
    os << "(disp: " << dispersion_rate << ", ext: " << extinction_rate << ")";
    return os.str();
  }

  inline double getDispersionRate(size_t from, size_t to) const {
    return dispersion_rate * (adjustment_matrix != nullptr
                                  ? (*adjustment_matrix)[from][to]
                                  : 1.0);
  }

  inline double getExtinctionRate() const { return extinction_rate; }
};

#endif  // LAGRANGE_COMMON_H__
