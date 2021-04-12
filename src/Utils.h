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

inline double bli_getiv_real(obj_t *obj, size_t idx) {
  double val;
  bli_getijm(idx, 0, obj, &val, nullptr);
  return val;
}

inline void bli_setiv_real(obj_t *obj, double val, size_t idx) {
  bli_setijm(val, 0, idx, 0, obj);
}

inline double bli_getijm_real(obj_t *obj, size_t i, size_t j) {
  double val;
  bli_getijm(i, j, obj, &val, nullptr);
  return val;
}

inline void bli_setijm_real(obj_t *obj, double val, size_t i, size_t j) {
  bli_setijm(val, 0, i, j, obj);
}

inline void bli_setv_real(obj_t *obj, double val) {
  obj_t alpha;
  bli_obj_create_1x1(BLIS_DOUBLE, &alpha);
  bli_setiv_real(&alpha, val, 0);

  bli_setv(&alpha, obj);

  bli_obj_free(&alpha);
}

inline obj_t bli_obj_wrap_scalar(double &alpha) {
  obj_t alpha_obj;
  bli_obj_create_1x1_with_attached_buffer(BLIS_DOUBLE, &alpha, &alpha_obj);
  return alpha_obj;
}

inline void bli_setm_scalar(double alpha, obj_t *obj) {
  auto a_obj = bli_obj_wrap_scalar(alpha);
  bli_setm(&a_obj, obj);
}

inline double bli_inf_normm(obj_t *A) {
  double norm = 0.0;
  auto norm_obj = bli_obj_wrap_scalar(norm);
  bli_normim(A, &norm_obj);
  return norm;
}

inline double bli_dotv_scalar(obj_t *a, obj_t *b) {
  double rho = 0.0;
  auto obj_rho = bli_obj_wrap_scalar(rho);

  bli_dotv(a, b, &obj_rho);

  return rho;
}

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
  double best_lh = bli_getiv_real(dist_lhs, 0);
  size_t states = bli_obj_length(dist_lhs);
  for (lagrange_dist_t i = 1; i < states; i++) {
    double cur_lh = bli_getiv_real(dist_lhs, i);
    if (best_lh < cur_lh) {
      best_dist = i;
      best_lh = cur_lh;
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
