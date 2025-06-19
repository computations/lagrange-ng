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

#include <logger.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#include "Common.hpp"

namespace lagrange {

struct RegionSplit {
  Range left;
  Range right;
};

inline auto lagrange_bextr(Range a, size_t i) -> uint64_t {
  return (a >> i) & 1ULL;
}

inline auto lagrange_popcount(Range a) -> size_t {
  return static_cast<size_t>(__builtin_popcountll(a));
}

constexpr auto lagrange_clz(Range a) -> size_t {
  assert(a != 0);
  return static_cast<size_t>(__builtin_clzll(a));
}

inline auto convert_vector_to_lagrange_dist(const std::vector<int> &vec_dist)
    -> Range {
  Range ret = 0;
  size_t area_count = vec_dist.size();

  for (size_t i = 0; i < area_count; ++i) {
    ret |= static_cast<Range>(vec_dist[i]) << i;
  }

  return ret;
}

/* Computes floor(log2) + 1 actually */
constexpr auto lagrange_fast_log2(size_t x) -> size_t {
  constexpr size_t BITS_IN_BYTE = 8;
  return sizeof(x) * BITS_IN_BYTE - lagrange_clz(x | 1);
}

constexpr auto lagrange_compute_region_mask(Range a) -> Range {
  return lagrange_fast_log2(a) - 1;
}

inline auto lagrange_compute_best_dist(const LagrangeConstColVector &dist_lhs,
                                       size_t states) -> Range {
  Range best_dist = 0;

  double best_lh = dist_lhs[0];

  for (Range i = 1; i < states; i++) {
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

auto compute_index_from_dist(Range i, size_t max_areas) -> size_t;

auto lagrange_convert_dist_string(Range dist,
                                  const std::vector<std::string> &names)
    -> std::string;

auto convert_dist_string_to_dist(const std::string &dist,
                                 const std::vector<std::string> &names)
    -> Range;

auto convert_dist_binary_string_to_dist(const std::string &dist) -> Range;

auto lagrange_parse_size_t(const std::string &str) -> size_t;

constexpr auto next_dist(Range d, uint32_t n) -> Range {
  d += 1;
  while (static_cast<size_t>(__builtin_popcountll(d)) > n) { d += 1; }
  return d;
}

/* Returns true if the dist "passes" the check. That is, if there are no common
 * bits set */
constexpr auto check_excl_dist(Range dist, Range excl_dist) -> bool {
  return (dist & excl_dist) == 0u;
}

constexpr auto check_incl_dist(Range dist, Range incl_dist) -> bool {
  return (dist & incl_dist) == incl_dist;
}

constexpr auto next_dist(Range d,
                         uint32_t n,
                         size_t index,
                         Range excl_area_mask = 0,
                         Range incl_area_mask = 0) -> std::pair<Range, size_t> {
  auto next = next_dist(d, n);
  index += 1;
  if (check_excl_dist(next, excl_area_mask)
      && check_incl_dist(next, incl_area_mask)) {
    return {next, index};
  }
  return next_dist(next, n, index, excl_area_mask, incl_area_mask);
}

/* Produces a dist -> index map */
inline auto invert_dist_map(size_t regions, size_t max_areas)
    -> std::vector<size_t> {
  size_t max_state = 1ULL << regions;
  std::vector<size_t> ret(max_state, std::numeric_limits<size_t>::max());

  size_t index = 0;
  Range dist = 0;

  for (index = 0, dist = 0; dist < max_state;
       ++index, dist = next_dist(dist, max_areas)) {
    ret[dist] = index;
  }

  return ret;
}

auto lagrange_convert_dist_to_list(Range dist,
                                   const std::vector<std::string> &names)
    -> std::vector<std::string>;

auto get_file_extension(const std::string &filename) -> std::string;

template <typename... Ts>
void lagrange_assert(bool condition,
                     std::basic_format_string<char> msg,
                     Ts... args) {
  if (!condition) {
    LOG_ERROR(msg, args...);
    abort();
  }
}

template <typename T>
class Option {
 public:
  Option() = default;

  explicit Option(const T &val) : _value{val}, _has_value{true} {}

  Option(const Option<T> &o) = default;
  Option& operator=(const Option<T> &o) = default;

  auto operator=(const T &v) -> Option<T> & {
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

  [[nodiscard]] auto hasValue() const -> bool { return _has_value; }

  operator bool() const { return hasValue(); }

 private:
  T _value;
  bool _has_value{false};
};

class UtilDistIndexConversionException : public std::exception {};
}  // namespace lagrange

namespace std {
template <>
struct hash<std::pair<size_t, double>> {
  auto operator()(std::pair<size_t, double> const &p) const -> std::size_t {
    return std::hash<size_t>{}(p.first) ^ std::hash<double>{}(p.second);
  }
};

}  // namespace std
#endif /* UTILS_H_ */
