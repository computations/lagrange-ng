/* Operation.cpp
 *
 * Created by Ben Bettisworth
 * 2020-10-27
 */

#include <blaze/Math.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/IdentityMatrix.h>

#include <array>
#include <memory>
#include <stdexcept>
#include <utility>

#include "Common.h"
#include "Operation.h"
#include "RateModel.h"
#include "Utils.h"
#include "Workspace.h"

constexpr inline size_t fast_log2(size_t x) { return __builtin_clzll(x); }

inline std::vector<lagrange_region_split_t> generate_splits(uint64_t state,
                                                            size_t regions) {
  uint64_t valid_region_mask = (1ull << regions) - 1;
  if (state == 0) {
    return {};
  }

  if (lagrange_popcount(state) == 1) {
    return {{state, state}};
  }

  std::vector<lagrange_region_split_t> ret;
  ret.reserve(regions);
  for (size_t i = 0; i < regions; ++i) {
    uint64_t x = 1ull << i;
    if ((state & x) == 0) {
      continue;
    }

    ret.push_back({x, state});
    ret.push_back({state, x});

    uint64_t y = (x ^ state) & valid_region_mask;
    ret.push_back({x, y});
    if (lagrange_popcount(y) > 1) {
      ret.push_back({y, x});
    }
  }
  return ret;
}

inline void weighted_combine(const lagrange_col_vector_t &c1,
                             const lagrange_col_vector_t &c2,
                             lagrange_col_vector_t &dest) {
  if (c1.size() != c2.size() && dest.size() == c1.size()) {
    throw std::runtime_error{"The vectors to combine are not equal sizes"};
  }
  size_t states = c1.size();
  size_t regions = fast_log2(states);

  for (size_t i = 0; i < states; i++) {
    auto splits = generate_splits(i, regions);
    for (auto p : splits) {
      dest[i] += c1[p.left] * c2[p.right];
    }
    if (splits.size() != 0) {
      dest[i] /= static_cast<double>(splits.size());
    }
  }
}

void MakeRateMatrixOperation::eval(Workspace &ws) const {
  auto &rm = ws.rate_matrix(_rate_matrix_index);
  for (lagrange_dist_t dist = 0; dist < ws.states(); dist++) {
    for (lagrange_dist_t i = 0; i < ws.regions(); i++) {
      if (lagrange_bextr(dist, i) != 0) {
        continue;
      }
      lagrange_dist_t gain_dist = dist | (1ul << i);
      rm(dist, gain_dist) = _dispersion_rate;
      rm(gain_dist, dist) = _extinction_rate;
    }
  }
}

void ExpmOperation::eval(std::shared_ptr<Workspace> ws) {
  lagrange_matrix_t A(ws->rate_matrix(_rate_matrix_index));

  size_t rows = A.rows();
  int scale_exp = std::max(0, 1 + static_cast<int>(blaze::linfNorm(A) * _t));
  A /= std::pow(2.0, scale_exp) * _t;
  // q is a magic parameter that controls the number of iterations of the loop
  // higher is more accurate, with each increase of q decreasing error by 4
  // orders of magnitude. Anything above 12 is probably snake oil.
  constexpr int q = 3;
  double c = 0.5;
  double sign = -1.0;
  blaze::IdentityMatrix<double, blaze::columnMajor> I(rows);

  lagrange_matrix_t X = _transposed ? blaze::trans(A) : A;
  lagrange_matrix_t N = I + c * (_transposed ? blaze::trans(A) : A);
  lagrange_matrix_t D = I - c * (_transposed ? blaze::trans(A) : A);

  // Using fortran indexing, and we started an iteration ahead to skip some
  // setup
  for (int i = 2; i <= q; ++i) {
    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    X = A * X;
    N += c * X;
    sign *= -1.0;
    D += sign * c * X;
  }
  A = blaze::inv(D) * N;
  for (int i = 0; i < scale_exp; ++i) {
    A *= A;
  }
  ws->prob_matrix(_prob_matrix_index) = std::move(A);
}

void DispersionOperation::eval(std::shared_ptr<Workspace> ws) {
  if (_expm_op != nullptr) {
    _expm_op->eval(ws);
  }
  ws->clv(_top_clv) = ws->prob_matrix(_prob_matrix_index) * ws->clv(_bot_clv);
}

void SplitOperation::eval(std::shared_ptr<Workspace> ws) const {
  for (auto &op : _lbranch_ops) {
    op->eval(ws);
  }
  for (auto &op : _rbranch_ops) {
    op->eval(ws);
  }

  auto &parent_clv = ws->clv(_parent_clv_index);
  auto &lchild_clv = ws->clv(_lbranch_clv_index);
  auto &rchild_clv = ws->clv(_rbranch_clv_index);

  weighted_combine(lchild_clv, rchild_clv, parent_clv);
}

void ReverseSplitOperation::eval(std::shared_ptr<Workspace> ws) const {
  for (auto &op : _branch_ops) {
    op->eval(ws);
  }

  if (_eval_clvs) {
    auto &bot_clv = ws->clv(_bot_clv_index);
    auto &ltop_clv = ws->clv(_ltop_clv_index);
    auto &rtop_clv = ws->clv(_rtop_clv_index);

    weighted_combine(ltop_clv, rtop_clv, bot_clv);
  }
}

double LHGoal::eval(std::shared_ptr<Workspace> ws) const {
  return blaze::dot(ws->clv(_root_clv_index),
                    ws->get_base_frequencies(_prior_index));
}
