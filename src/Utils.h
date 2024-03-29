/*
 * Utils.h
 *
 *  Created on: Mar 10, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef UTILS_H
#define UTILS_H

#include <stdexcept>
#include <string>
#include <vector>

#include "Common.h"

struct lagrange_region_split_t {
  lagrange_dist_t left;
  lagrange_dist_t right;
};

inline auto lagrange_bextr(lagrange_dist_t a, size_t i) -> uint64_t {
  return (a >> i) & 1ULL;
}

inline auto lagrange_popcount(lagrange_dist_t a) -> size_t {
  return static_cast<size_t>(__builtin_popcountll(a));
}

constexpr inline auto lagrange_clz(lagrange_dist_t a) -> size_t {
  return static_cast<size_t>(__builtin_clzll(a));
}

inline auto convert_vector_to_lagrange_dist(const std::vector<int> &vec_dist)
    -> lagrange_dist_t {
  lagrange_dist_t ret = 0;
  size_t area_count = vec_dist.size();

  for (size_t i = 0; i < area_count; ++i) {
    ret |= static_cast<lagrange_dist_t>(vec_dist[i]) << i;
  }

  return ret;
}

/* Computes floor(log2) + 1 actually */
constexpr inline auto lagrange_fast_log2(size_t x) -> size_t {
  constexpr size_t BITS_IN_BYTE = 8;
  return sizeof(x) * BITS_IN_BYTE - lagrange_clz(x | 1);
}

inline auto lagrange_compute_best_dist(
    const lagrange_const_col_vector_t &dist_lhs, size_t states)
    -> lagrange_dist_t {
  lagrange_dist_t best_dist = 0;

  double best_lh = dist_lhs[0];

  for (lagrange_dist_t i = 1; i < states; i++) {
    double cur_lh = dist_lhs[i];
    if (best_lh < cur_lh) {
      best_dist = i;
      best_lh = cur_lh;
    }
  }

  return best_dist;
}

auto lagrange_compute_restricted_state_count(size_t regions, size_t max_areas)
    -> size_t;

auto compute_index_from_dist(lagrange_dist_t i, size_t max_areas) -> size_t;

auto lagrange_convert_dist_string(lagrange_dist_t dist,
                                  const std::vector<std::string> &names)
    -> std::string;

auto lagrange_convert_dist_string_to_dist(const std::string &dist,
                                          const std::vector<std::string> &names)
    -> lagrange_dist_t;

auto lagrange_convert_dist_binary_string_to_dist(const std::string &dist)
    -> lagrange_dist_t;

auto lagrange_parse_size_t(const std::string &str) -> size_t;

constexpr inline auto next_dist(lagrange_dist_t d, uint32_t n)
    -> lagrange_dist_t {
  d += 1;
  while (static_cast<size_t>(__builtin_popcountll(d)) > n) { d += 1; }
  return d;
}

/* Returns true if the dist "passes" the check. That is, if there are no common
 * bits set */
constexpr inline auto check_excl_dist(lagrange_dist_t dist,
                                      lagrange_dist_t excl_dist) -> bool {
  return !(dist & excl_dist);
}

constexpr inline auto check_incl_dist(lagrange_dist_t dist,
                                      lagrange_dist_t incl_dist) -> bool {
  return (dist & incl_dist) == incl_dist;
}

constexpr inline auto next_dist(lagrange_dist_t d, uint32_t n, size_t index,
                                lagrange_dist_t excl_area_mask = 0,
                                lagrange_dist_t incl_area_mask = 0)
    -> std::pair<lagrange_dist_t, size_t> {
  auto next = next_dist(d, n);
  index += 1;
  if (check_excl_dist(next, excl_area_mask) &&
      check_incl_dist(next, incl_area_mask)) {
    return {next, index};
  }
  return next_dist(next, n, index, excl_area_mask, incl_area_mask);
}

std::vector<std::string> lagrange_convert_dist_to_list(
    lagrange_dist_t dist, const std::vector<std::string> &names);

std::string get_file_extension(const std::string &filename);

template <typename T>
class lagrange_option_t {
 public:
  lagrange_option_t() = default;
  explicit lagrange_option_t(const T &val) : _value{val}, _has_value{true} {}

  lagrange_option_t(const lagrange_option_t<T> &o) = default;

  auto operator=(const T &v) -> lagrange_option_t<T> & {
    _value = v;
    _has_value = true;
    return *this;
  }

  auto get() -> T & {
    if (_has_value) { return _value; }
    throw std::runtime_error{"lagrange_option_t has no value when used"};
  }

  auto get() const -> const T & {
    if (_has_value) { return _value; }
    throw std::runtime_error{"lagrange_option_t has no value when used"};
  }

  auto get(const T &def) const -> T {
    if (_has_value) { return _value; }
    return def;
  }

  auto has_value() const -> bool { return _has_value; }

  operator bool() const { return has_value(); }

 private:
  T _value;
  bool _has_value{false};
};

namespace std {
template <>
struct hash<std::pair<size_t, double>> {
  auto operator()(std::pair<size_t, double> const &p) const -> std::size_t {
    return std::hash<size_t>{}(p.first) ^ std::hash<double>{}(p.second);
  }
};

class lagrange_util_dist_index_conversion_exception : public std::exception {};

}  // namespace std
#endif /* UTILS_H_ */
