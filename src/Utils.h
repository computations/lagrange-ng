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

inline double bli_getiv_real(obj_t *obj, size_t idx) {
  double val;
  double ival;
  bli_getijm(idx, 0, obj, &val, &ival);
  return val;
}

inline void bli_setiv_real(obj_t *obj, double val, size_t idx) {
  bli_setijm(val, 0, idx, 0, obj);
}

inline double bli_getijm_real(obj_t *obj, size_t i, size_t j) {
  double val;
  double ival = 0.0;
  bli_getijm(i, j, obj, &val, &ival);
  return val;
}

inline void bli_setijm_real(obj_t *obj, double val, size_t i, size_t j) {
  bli_setijm(val, 0, i, j, obj);
}

inline void bli_setv_real(obj_t *obj, double val) {
  obj_t alpha;
  bli_obj_create_1x1_with_attached_buffer(BLIS_DOUBLE, &val, &alpha);
  bli_setiv_real(&alpha, val, 0);

  bli_setv(&alpha, obj);
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

inline void lagrange_set_bli_vector(obj_t *obj, std::vector<double> vec) {
  bli_obj_create(BLIS_DOUBLE, vec.size(), 1, 0, 0, obj);
  for (size_t i = 0; i < vec.size(); i++) { bli_setiv_real(obj, vec[i], i); }
}

inline void lagrange_set_bli_matrix(obj_t *obj, size_t rows, size_t cols,
                                    std::vector<double> vec) {
  bli_obj_create(BLIS_DOUBLE, rows, cols, 0, 0, obj);
  if (rows * cols != vec.size()) {
    throw std::runtime_error{
        "Dims did not match the provided vector for set_matrix"};
  }
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      bli_setijm_real(obj, vec[i * rows + j], i, j);
    }
  }
}

inline uint64_t lagrange_bextr(lagrange_dist_t a, size_t i) {
  return (a >> i) & 1ull;
}

inline size_t lagrange_popcount(lagrange_dist_t a) {
  return __builtin_popcountll(a);
}

inline size_t lagrange_clz(lagrange_dist_t a) { return __builtin_clzll(a); }

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

}  // namespace std
#endif /* UTILS_H_ */
