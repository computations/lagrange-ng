/*
 * Utils.h
 *
 *  Created on: Mar 10, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <stdexcept>
#include <string>
#include <vector>

#include "Common.h"

struct lagrange_region_split_t {
  lagrange_dist_t left;
  lagrange_dist_t right;
};

void Tokenize(const std::string &str, std::vector<std::string> &tokens,
              const std::string &delimiters = " ");
void TrimSpaces(std::string &str);

inline uint64_t lagrange_bextr(lagrange_dist_t a, size_t i) {
  return (a >> i) & 1ull;
}

inline size_t lagrange_popcount(lagrange_dist_t a) {
  return static_cast<size_t>(__builtin_popcountll(a));
}

constexpr inline size_t lagrange_clz(lagrange_dist_t a) {
  return static_cast<size_t>(__builtin_clzll(a));
}

inline lagrange_dist_t convert_vector_to_lagrange_dist(
    const std::vector<int> &vec_dist) {
  lagrange_dist_t ret = 0;
  size_t area_count = vec_dist.size();

  for (size_t i = 0; i < area_count; ++i) {
    ret |= static_cast<lagrange_dist_t>(vec_dist[i]) << i;
  }

  return ret;
}

constexpr inline size_t lagrange_fast_log2(size_t x) {
  return sizeof(x) * 8 - lagrange_clz(x | 1);
}

inline lagrange_dist_t lagrange_compute_best_dist(
    const lagrange_col_vector_t &dist_lhs, size_t states) {
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

size_t lagrange_compute_restricted_state_count(size_t regions,
                                               size_t max_areas);

size_t compute_index_from_dist(lagrange_dist_t i, size_t max_areas);

std::string lagrange_convert_dist_string(lagrange_dist_t dist,
                                         const std::vector<std::string> &names);

size_t lagrange_parse_size_t(const std::string &str);

template <typename T>
class lagrange_option_t {
 public:
  lagrange_option_t() : _has_value{false} {}
  lagrange_option_t(const T &val) : _value{val}, _has_value{true} {}

  T &get() {
    if (_has_value) { return _value; }
    throw std::runtime_error{"lagrange_option_t has no value when used"};
  }

  const T &get() const {
    if (_has_value) { return _value; }
    throw std::runtime_error{"lagrange_option_t has no value when used"};
  }

  bool has_value() const { return _has_value; }

  operator bool() const { return has_value(); }

 private:
  T _value;
  bool _has_value;
};

namespace std {
template <>
struct hash<std::pair<size_t, double>> {
  std::size_t operator()(std::pair<size_t, double> const &p) const {
    return std::hash<size_t>{}(p.first) ^ std::hash<double>{}(p.second);
  }
};

class lagrange_util_dist_index_conversion_exception : public std::exception {};

}  // namespace std
#endif /* UTILS_H_ */
