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
#include <cstddef>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "AncSplit.h"
#include "Common.h"
#include "Operation.h"
#include "RateModel.h"
#include "Utils.h"
#include "Workspace.h"

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
                             const std::vector<lagrange_dist_t> excl_dists,
                             lagrange_col_vector_t &dest) {
  if (c1.size() != c2.size() && dest.size() == c1.size()) {
    throw std::runtime_error{"The vectors to combine are not equal sizes"};
  }
  size_t states = c1.size();
  size_t regions = lagrange_fast_log2(states);

  size_t idx_excl = 0;

  dest = 0.0;

  for (size_t i = 0; i < states; i++) {
    while (idx_excl < excl_dists.size() && i > excl_dists[idx_excl]) {
      idx_excl++;
    }
    if (idx_excl < excl_dists.size() && i == excl_dists[idx_excl]) {
      idx_excl++;
      continue;
    }
    auto splits = generate_splits(i, regions);
    for (auto p : splits) {
      dest[i] += c1[p.left] * c2[p.right];
    }
    if (splits.size() != 0) {
      dest[i] /= static_cast<double>(splits.size());
    }
  }
}

inline std::string make_tabs(size_t tabLevel) {
  std::ostringstream tabs;
  tabs << "|";
  for (size_t i = 0; i < tabLevel; ++i) {
    tabs << " |";
  }
  tabs << " ";
  return tabs.str();
}

inline std::string opening_line(const std::string &tabs) {
  std::ostringstream line;

  line << tabs.substr(0, tabs.size() - 2) << "┌";
  size_t cur_len = line.str().size();
  for (size_t i = 0; i < (80 - cur_len); ++i) {
    line << "─";
  }

  return line.str();
}

inline std::string closing_line(const std::string &tabs) {
  std::ostringstream line;

  line << tabs.substr(0, tabs.size() - 2) << "└";
  size_t cur_len = line.str().size();
  for (size_t i = 0; i < (80 - cur_len); ++i) {
    line << "─";
  }

  return line.str();
}

void MakeRateMatrixOperation::eval(std::shared_ptr<Workspace> ws) {
  if (_last_update < _last_execution) {
    return;
  }
  auto &rm = ws->rate_matrix(_rate_matrix_index);
  rm = 0.0;
  auto period = ws->get_period_params(_period_index);
  double dispersion_rate = period.dispersion_rate;
  double extinction_rate = period.extinction_rate;
  for (lagrange_dist_t dist = 0; dist < ws->states(); dist++) {
    for (lagrange_dist_t i = 0; i < ws->regions(); i++) {
      if (lagrange_bextr(dist, i) != 0) {
        continue;
      }
      lagrange_dist_t gain_dist = dist | (1ul << i);
      rm(gain_dist, dist) = extinction_rate;
      if (dist != 0) {
        rm(dist, gain_dist) = dispersion_rate;
      }
    }
  }

  for (size_t i = 0; i < rm.rows(); i++) {
    double sum = 0;
    for (size_t j = 0; j < rm.columns(); j++) {
      sum += rm(i, j);
    }
    rm(i, i) = -sum;
  }

  _last_execution = ws->advance_clock();
}

void MakeRateMatrixOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                          std::ostream &os,
                                          size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "MakeRateMatrixOperation:\n"
     << tabs << "Rate Matrix (index: " << _rate_matrix_index << "):\n";

  auto &rm = ws->rate_matrix(_rate_matrix_index);

  for (size_t i = 0; i < rm.rows(); ++i) {
    auto row = blaze::row(rm, i);
    os << tabs << std::setprecision(10) << row;
  }

  os << tabs << "Period index: " << _period_index << "\n"
     << tabs << ws->get_period_params(_period_index).toString() << "\n"
     << tabs << "_last_execution: " << _last_execution << "\n"
     << tabs << "_last_update: " << _last_update << "\n";
  os << closing_line(tabs);
}

std::string MakeRateMatrixOperation::printStatus(
    const std::shared_ptr<Workspace> &ws, size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void ExpmOperation::eval(std::shared_ptr<Workspace> ws) {
  if (_rate_matrix_op != nullptr &&
      _last_execution > _rate_matrix_op->last_update()) {
    return;
  }
  _rate_matrix_op->eval(ws);

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
  _last_execution = ws->advance_clock();
  ws->prob_matrix(_prob_matrix_index) = A;
}

void ExpmOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                std::ostream &os, size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "ExpmOperation:\n"
     << tabs << "Rate Matrix (index: " << _rate_matrix_index << "):\n";

  auto &rm = ws->rate_matrix(_rate_matrix_index);
  for (size_t i = 0; i < rm.rows(); ++i) {
    auto row = blaze::row(rm, i);
    os << tabs << std::setprecision(10) << row;
  }

  os << tabs << "Prob Matrix (index: " << _prob_matrix_index << "):\n";

  auto &pm = ws->prob_matrix(_prob_matrix_index);

  for (size_t i = 0; i < pm.rows(); ++i) {
    auto row = blaze::row(pm, i);
    os << tabs << std::setprecision(10) << row;
  }

  os << tabs << "t: " << std::setprecision(10) << _t << "\n";
  os << tabs << "_last_execution: " << _last_execution;
  if (_rate_matrix_op != nullptr) {
    os << "\n" << _rate_matrix_op->printStatus(ws, tabLevel + 1) << "\n";
  }
  os << closing_line(tabs);
}

std::string ExpmOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                       size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void DispersionOperation::eval(std::shared_ptr<Workspace> ws) {
  if (_expm_op != nullptr) {
    _expm_op->eval(ws);
  }

  ws->clv(_top_clv) = ws->prob_matrix(_prob_matrix_index) * ws->clv(_bot_clv);
}

void DispersionOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                      std::ostream &os, size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);

  os << opening_line(tabs) << "\n";
  os << tabs << "DispersionOperation:\n";
  os << tabs << "Top clv (index: " << _top_clv << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_top_clv));
  os << tabs << "Bot clv (index: " << _bot_clv << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_bot_clv));

  os << tabs << "Prob Matrix (index: " << _prob_matrix_index << "):\n";

  auto &pm = ws->prob_matrix(_prob_matrix_index);

  for (size_t i = 0; i < pm.rows(); ++i) {
    auto row = blaze::row(pm, i);
    os << tabs << std::setprecision(10) << row;
  }

  if (_expm_op != nullptr) {
    os << _expm_op->printStatus(ws, tabLevel + 1);
  }
  os << "\n";
  os << closing_line(tabs);
}

std::string DispersionOperation::printStatus(
    const std::shared_ptr<Workspace> &ws, size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
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

  weighted_combine(lchild_clv, rchild_clv, _excl_dists, parent_clv);
}

void SplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                 std::ostream &os, size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "SplitOperation:\n";

  os << tabs << "Lbranch clv (index: " << _lbranch_clv_index
     << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_lbranch_clv_index));
  os << tabs << "Rbranch clv (index: " << _rbranch_clv_index
     << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_rbranch_clv_index));
  os << tabs << "Parent clv (index: " << _parent_clv_index
     << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_parent_clv_index));

  if (_excl_dists.size() != 0) {
    os << tabs << "Excluded dists:\n";
    for (size_t i = 0; i < _excl_dists.size(); ++i) {
      os << tabs << _excl_dists[i] << "\n";
    }
  }

  os << tabs << "Left Branch ops:\n";
  if (_lbranch_ops.size() != 0) {
    for (auto &op : _lbranch_ops) {
      os << op->printStatus(ws, tabLevel + 1);
    }
  }
  os << "\n" << tabs << "Right Branch ops:\n";
  if (_rbranch_ops.size() != 0) {
    for (auto &op : _rbranch_ops) {
      os << op->printStatus(ws, tabLevel + 1);
    }
  }
  os << "\n";
  os << closing_line(tabs);
}

std::string SplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                        size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void ReverseSplitOperation::eval(std::shared_ptr<Workspace> ws) const {
  for (auto &op : _branch_ops) {
    op->eval(ws);
  }

  if (_eval_clvs) {
    auto &bot_clv = ws->clv(_bot_clv_index);
    auto &ltop_clv = ws->clv(_ltop_clv_index);
    auto &rtop_clv = ws->clv(_rtop_clv_index);

    weighted_combine(ltop_clv, rtop_clv, _excl_dists, bot_clv);
  }
}

void ReverseSplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                        std::ostream &os,
                                        size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "SplitOperation:\n";

  os << tabs << "Bot clv (index: " << _bot_clv_index
     << "): " << std::setprecision(10) << blaze::trans(ws->clv(_bot_clv_index));
  os << tabs << "Ltop clv (index: " << _ltop_clv_index
     << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_ltop_clv_index));
  os << tabs << "Rtop clv (index: " << _rtop_clv_index
     << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_rtop_clv_index));

  if (_excl_dists.size() != 0) {
    os << tabs << "Excluded dists:\n";
    for (size_t i = 0; i < _excl_dists.size(); ++i) {
      os << tabs << _excl_dists[i] << "\n";
    }
  }

  os << tabs << "Eval CLVS: " << _eval_clvs << "\n";

  os << tabs << "Branch ops:\n";
  if (_branch_ops.size() != 0) {
    for (auto &op : _branch_ops) {
      os << op->printStatus(ws, tabLevel + 1);
    }
  }
  os << "\n";
  os << closing_line(tabs);
}

std::string ReverseSplitOperation::printStatus(
    const std::shared_ptr<Workspace> &ws, size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

double LHGoal::eval(std::shared_ptr<Workspace> ws) const {
  return blaze::dot(ws->clv(_root_clv_index),
                    ws->get_base_frequencies(_prior_index));
}

lagrange_col_vector_t StateLHGoal::eval(std::shared_ptr<Workspace> ws) const {
  lagrange_col_vector_t tmp;
  weighted_combine(ws->clv(_lchild_clv_index), ws->clv(_rchild_clv_index),
                   /*excl_dists=*/{}, tmp);
  return ws->clv(_parent_clv_index) * tmp;
}

std::unordered_map<lagrange_dist_t, std::vector<AncSplit>> SplitLHGoal::eval(
    std::shared_ptr<Workspace> ws) const {
  std::unordered_map<lagrange_dist_t, std::vector<AncSplit>> ret;

  auto &parent_clv = ws->clv(_parent_clv_index);
  auto &lchild_clv = ws->clv(_lchild_clv_index);
  auto &rchild_clv = ws->clv(_rchild_clv_index);
  for (lagrange_dist_t dist = 0; dist < ws->states(); dist++) {
    std::vector<AncSplit> anc_split_vec;
    auto splits = generate_splits(dist, ws->regions());
    double weight = 1.0 / splits.size();
    for (auto sp : splits) {
      AncSplit anc_split(dist, sp.left, sp.right, weight);
      double lh = parent_clv[dist] * lchild_clv[sp.left] *
                  rchild_clv[sp.right] * weight;
      anc_split.setLikelihood(lh);
      anc_split_vec.push_back(anc_split);
    }
    ret[dist] = anc_split_vec;
  }
  return ret;
}
