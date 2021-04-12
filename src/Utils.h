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

#include <map>
#include <string>
#include <vector>

#include "Common.h"

using namespace std;

struct lagrange_region_split_t {
  lagrange_dist_t left;
  lagrange_dist_t right;
};

void Tokenize(const string &str, vector<string> &tokens,
              const string &delimiters = " ");
void TrimSpaces(string &str);

inline uint64_t lagrange_bextr(lagrange_dist_t a, size_t i) {
  return (a >> i) & 1ull;
}

inline size_t lagrange_popcount(lagrange_dist_t a) {
  return __builtin_popcountll(a);
}

inline size_t lagrange_clz(lagrange_dist_t a) { return __builtin_clzll(a); }

inline lagrange_dist_t convert_vector_to_lagrange_dist(
    const vector<int> &vec_dist) {
  lagrange_dist_t ret = 0;
  size_t area_count = vec_dist.size();

  for (size_t i = 0; i < area_count; ++i) {
    ret |= static_cast<lagrange_dist_t>(vec_dist[i]) << i;
  }

  return ret;
}

constexpr inline size_t lagrange_fast_log2(size_t x) {
  return __builtin_clzll(x);
}

inline lagrange_dist_t lagrange_compute_best_dist(
    const lagrange_col_vector_t &dist_lhs) {
  lagrange_dist_t best_dist = 0;
  double best_lh = dist_lhs[0];
  for (lagrange_dist_t i = 1; i < static_cast<unsigned long>(dist_lhs.size());
       i++) {
    if (best_lh < dist_lhs[i]) {
      best_dist = i;
      best_lh = dist_lhs[i];
    }
  }
  return best_dist;
}

std::string lagrange_convert_dist_string(lagrange_dist_t dist,
                                         const std::vector<std::string> &names);

namespace std {
template <>
struct hash<std::pair<size_t, double>> {
  std::size_t operator()(std::pair<size_t, double> const &p) const {
    return std::hash<size_t>{}(p.first) ^ std::hash<double>{}(p.second);
  }
};

}  // namespace std
#endif /* UTILS_H_ */
