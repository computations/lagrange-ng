#include "Common.h"
#include "Operation.h"
#include "RateModel.h"
#include "Utils.h"
#include <array>
#include <blaze/math/IdentityMatrix.h>
#include <memory>
#include <utility>

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
  lagrange_matrix_t X = A;
  lagrange_matrix_t N = I + c * A;
  lagrange_matrix_t D = I - c * A;

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

std::vector<lagrange_region_split_t>
SplitOperation::generate_splits(uint64_t state, size_t regions) {
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

void SplitOperation::eval(std::shared_ptr<Workspace> ws) {
  for (size_t i = 0; i < _lbranch_ops_count; i++) {
    _lbranch_ops.get()[i].eval(ws);
  }
  for (size_t i = 0; i < _rbranch_ops_count; i++) {
    _rbranch_ops.get()[i].eval(ws);
  }

  auto &parent_clv = ws->clv(_parent_clv_index);
  auto &lchild_clv = ws->clv(_lbranch_clv_index);
  auto &rchild_clv = ws->clv(_rbranch_clv_index);
  for (size_t i = 0; i < ws->states(); i++) {
    for (auto p : generate_splits(i, ws->regions())) {
      parent_clv[i] += rchild_clv[p.left] * lchild_clv[p.right];
    }
  }
}
