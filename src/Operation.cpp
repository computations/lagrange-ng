/* Operation.cpp
 *
 * Created by Ben Bettisworth
 * 2020-10-27
 */

#include <array>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "AncSplit.h"
#include "Common.h"
#include "Operation.h"
#include "Utils.h"
#include "Workspace.h"

inline void generate_splits(uint64_t state, size_t regions,
                            std::vector<lagrange_region_split_t> &results) {
  results.clear();
  uint64_t valid_region_mask = (1ull << regions) - 1;
  if (state == 0) { return; }

  if (lagrange_popcount(state) == 1) {
    results.push_back({state, state});
    return;
  }

  results.reserve(regions);
  for (size_t i = 0; i < regions; ++i) {
    uint64_t x = 1ull << i;
    if ((state & x) == 0) { continue; }

    results.push_back({x, state});
    results.push_back({state, x});

    uint64_t y = (x ^ state) & valid_region_mask;
    results.push_back({x, y});
    if (lagrange_popcount(y) > 1) { results.push_back({y, x}); }
  }
}

inline void weighted_combine(const lagrange_col_vector_t &c1,
                             const lagrange_col_vector_t &c2,
                             const std::vector<lagrange_dist_t> excl_dists,
                             lagrange_col_vector_t &dest, size_t c1_scale,
                             size_t c2_scale, size_t &scale_count) {
  if (c1.size() != c2.size() && dest.size() == c1.size()) {
    throw std::runtime_error{"The vectors to combine are not equal sizes"};
  }
  size_t states = c1.size();
  size_t regions = lagrange_fast_log2(states);

  size_t idx_excl = 0;

  dest = 0.0;

  scale_count = c1_scale + c2_scale;
  bool scale = true;

  std::vector<lagrange_region_split_t> splits;

  for (size_t i = 0; i < states; i++) {
    while (idx_excl < excl_dists.size() && i > excl_dists[idx_excl]) {
      idx_excl++;
    }
    if (idx_excl < excl_dists.size() && i == excl_dists[idx_excl]) {
      idx_excl++;
      continue;
    }
    generate_splits(i, regions, splits);
    for (auto &p : splits) { dest[i] += c1[p.left] * c2[p.right]; }
    if (splits.size() != 0) { dest[i] /= static_cast<double>(splits.size()); }
    if (dest[i] < lagrange_scale_threshold) {
      scale &= true;
    } else {
      scale &= false;
    }
  }
  if (scale) {
    dest *= lagrange_scaling_factor;
    scale_count += 1;
  }
}

inline void reverse_weighted_combine(
    const lagrange_col_vector_t &c1, const lagrange_col_vector_t &c2,
    const std::vector<lagrange_dist_t> excl_dists,
    lagrange_col_vector_t &dest) {
  if (c1.size() != c2.size() && dest.size() == c1.size()) {
    throw std::runtime_error{"The vectors to combine are not equal sizes"};
  }
  size_t states = c1.size();
  size_t regions = lagrange_fast_log2(states);

  size_t idx_excl = 0;

  dest = 0.0;

  std::vector<lagrange_region_split_t> splits;

  for (size_t i = 0; i < states; i++) {
    while (idx_excl < excl_dists.size() && i > excl_dists[idx_excl]) {
      idx_excl++;
    }
    if (idx_excl < excl_dists.size() && i == excl_dists[idx_excl]) {
      idx_excl++;
      continue;
    }

    generate_splits(i, regions, splits);
    if (splits.size() == 0) { continue; }

    double weight = 1.0 / static_cast<double>(splits.size());
    for (auto &p : splits) { dest[p.left] += c1[i] * c2[p.right] * weight; }
  }
}

inline std::string make_tabs(size_t tabLevel) {
  std::ostringstream tabs;
  tabs << "|";
  for (size_t i = 0; i < tabLevel; ++i) { tabs << " |"; }
  tabs << " ";
  return tabs.str();
}

inline std::string opening_line(const std::string &tabs) {
  std::ostringstream line;

  line << tabs.substr(0, tabs.size() - 2) << "┌";
  size_t cur_len = line.str().size();
  for (size_t i = 0; i < (80 - cur_len); ++i) { line << "─"; }

  return line.str();
}

inline std::string closing_line(const std::string &tabs) {
  std::ostringstream line;

  line << tabs.substr(0, tabs.size() - 2) << "└";
  size_t cur_len = line.str().size();
  for (size_t i = 0; i < (80 - cur_len); ++i) { line << "─"; }

  return line.str();
}

void MakeRateMatrixOperation::eval(std::shared_ptr<Workspace> ws) {
  // std::lock_guard<std::mutex> t_lock(*_lock);
  auto &rm = ws->rate_matrix(_rate_matrix_index);
  rm = 0.0;
  auto &period = ws->get_period_params(_period_index);
  for (lagrange_dist_t dist = 0; dist < ws->states(); dist++) {
    for (lagrange_dist_t i = 0; i < ws->regions(); i++) {
      if (lagrange_bextr(dist, i) != 0) { continue; }

      lagrange_dist_t gain_dist = dist | (1ul << i);
      rm(gain_dist, dist) = period.getExtinctionRate();

      if (dist == 0) { continue; }

      for (size_t j = 0; j < ws->regions(); ++j) {
        rm(dist, gain_dist) +=
            period.getDispersionRate(j, i) * lagrange_bextr(dist, j);
      }
    }
  }

  for (size_t i = 0; i < rm.rows(); i++) {
    double sum = 0;

    for (size_t j = 0; j < rm.columns(); j++) { sum += rm(i, j); }

    rm(i, i) = -sum;
  }

  _last_execution = ws->advance_clock();
  ws->update_rate_matrix_clock(_rate_matrix_index);
}

void MakeRateMatrixOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                          std::ostream &os,
                                          size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "MakeRateMatrixOperation:\n"
     << tabs << "Rate Matrix (index: " << _rate_matrix_index
     << " update: " << ws->last_update_rate_matrix(_rate_matrix_index)
     << "):\n";

  auto &rm = ws->rate_matrix(_rate_matrix_index);

  for (size_t i = 0; i < rm.rows(); ++i) {
    auto row = blaze::row(rm, i);
    os << tabs << std::setprecision(10) << row;
  }

  os << tabs << "Period index: " << _period_index << "\n"
     << tabs << ws->get_period_params(_period_index).toString() << "\n"
     << tabs << "_last_execution: " << _last_execution << "\n"
     << tabs << "_last_update: " << _last_update << "\n"
     << tabs << "correct: " << std::boolalpha << correct(ws) << "\n";
  os << closing_line(tabs);
}

std::string MakeRateMatrixOperation::printStatus(
    const std::shared_ptr<Workspace> &ws, size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void ExpmOperation::eval(std::shared_ptr<Workspace> ws) {
  if ((_rate_matrix_op != nullptr) &&
      (_last_execution > _rate_matrix_op->last_update())) {
    return;
  }

  lagrange_matrix_t A(ws->rate_matrix(_rate_matrix_index));

  size_t rows = A.rows();
  // We place an arbitrary limit on the size of scale exp because if it is too
  // large we run into numerical issues.
  int scale_exp =
      std::min(30, std::max(0, 1 + static_cast<int>(blaze::linfNorm(A) * _t)));

  A /= std::pow(2.0, scale_exp) / _t;
  // q is a magic parameter that controls the number of iterations of the loop
  // higher is more accurate, with each increase of q decreasing error by 4
  // orders of magnitude. Anything above 12 is probably snake oil.
  constexpr int q = 3;
  double c = 0.5;
  double sign = -1.0;
  blaze::IdentityMatrix<double, blaze::columnMajor> I(rows);

#if 0
  lagrange_matrix_t X = A;
  lagrange_matrix_t N = I + c * X;
  lagrange_matrix_t D = I - c * X;
#endif

  lagrange_matrix_t X_2 = A;
  lagrange_matrix_t N = I + c * X_2;
  lagrange_matrix_t D = I - c * X_2;
  lagrange_matrix_t X_1;
  // Using fortran indexing, and we started an iteration ahead to skip some
  // setup. Furhthermore, we are going to unroll the loop to allow us to skip
  // some assignments.
#if 0
  for (int i = 2; i <= q; i++) {
    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    std::cout << "C: " << c << std::endl;
    X = A * X;
    N += c * X;
    sign *= -1.0;
    D += sign * c * X;
  }
#endif

  for (int i = 2; i <= q;) {
    c = c * (q - i + 1) / (i * (2 * q - i + 1));

    X_1 = A * X_2;
    N += c * X_1;
    sign *= -1.0;
    D += sign * c * X_1;

    i += 1;

    if (i > q) { break; }

    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    X_2 = A * X_1;
    N += c * X_2;
    sign *= -1.0;
    D += sign * c * X_2;

    i += 1;
  }

  X_1 = blaze::solve(D, N);
  auto &rX_1 = X_1;
  auto &rX_2 = X_2;
  for (int i = 0; i < scale_exp; ++i) {
    rX_2 = rX_1 * rX_1;
    std::swap(rX_1, rX_2);
  }

#if 0
  X = blaze::solve(D, N);
  for (int i = 0; i < scale_exp; ++i) {
    X = X;
  }
#endif

  if (_transposed) {
    blaze::transpose(X_1);
    blaze::row(X_1, 0) = 0.0;
    X_1(0, 0) = 1.0;
  }

  _last_execution = ws->advance_clock();
  ws->update_prob_matrix(_prob_matrix_index, X_1);
}

void ExpmOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                std::ostream &os, size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "ExpmOperation:\n"
     << tabs << "Rate Matrix (index: " << _rate_matrix_index << "):\n";

  os << tabs << "Prob Matrix (index: " << _prob_matrix_index
     << " update: " << ws->last_update_prob_matrix(_prob_matrix_index)
     << "):\n";

  auto &pm = ws->prob_matrix(_prob_matrix_index);

  for (size_t i = 0; i < pm.rows(); ++i) {
    auto row = blaze::row(pm, i);
    os << tabs << std::setprecision(10) << row;
  }

  os << tabs << "t: " << std::setprecision(16) << _t << "\n";
  os << tabs << "_last_execution: " << _last_execution << "\n";
  os << tabs << "correct:  " << std::boolalpha << correct(ws);
  if (_rate_matrix_op != nullptr) {
    os << "\n" << _rate_matrix_op->printStatus(ws, tabLevel + 1) << "\n";
  } else {
    auto &rm = ws->rate_matrix(_rate_matrix_index);
    for (size_t i = 0; i < rm.rows(); ++i) {
      auto row = blaze::row(rm, i);
      os << tabs << std::setprecision(10) << row;
    }
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
    if (_expm_op.use_count() > 1) {
      std::lock_guard<std::mutex>(_expm_op->getLock());
      _expm_op->eval(ws);
    } else {
      _expm_op->eval(ws);
    }
  }

  if (ws->last_update_clv(_bot_clv) < ws->last_update_clv(_top_clv)) {
    /* we have already computed stuff, return */
    return;
  }

  _last_execution = ws->advance_clock();

  ws->update_clv(
      std::move(ws->prob_matrix(_prob_matrix_index) * ws->clv(_bot_clv)),
      _top_clv);
  ws->clv_scalar(_top_clv) = ws->clv_scalar(_bot_clv);
}

void DispersionOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                      std::ostream &os, size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);

  os << opening_line(tabs) << "\n";
  os << tabs << "DispersionOperation:\n";
  os << tabs << "Top clv (index: " << _top_clv
     << " update: " << ws->last_update_clv(_top_clv)
     << "): " << std::setprecision(10) << blaze::trans(ws->clv(_top_clv));
  os << tabs << "Bot clv (index: " << _bot_clv
     << " update: " << ws->last_update_clv(_bot_clv)
     << "): " << std::setprecision(10) << blaze::trans(ws->clv(_bot_clv));
  os << tabs << "Last Executed: " << _last_execution << "\n";
  os << tabs << "correct: " << std::boolalpha << correct(ws) << "\n";
  os << tabs << "Prob Matrix (index: " << _prob_matrix_index << ")\n";

  if (_expm_op != nullptr) { os << _expm_op->printStatus(ws, tabLevel + 1); }
  os << "\n";
  os << closing_line(tabs);
}

std::string DispersionOperation::printStatus(
    const std::shared_ptr<Workspace> &ws, size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void SplitOperation::eval(std::shared_ptr<Workspace> ws) {
  for (auto &op : _lbranch_ops) {
    if (op.use_count() > 1) {
      std::lock_guard<std::mutex> lock(op->getLock());
      op->eval(ws);
    } else {
      op->eval(ws);
    }
  }
  for (auto &op : _rbranch_ops) {
    if (op.use_count() > 1) {
      std::lock_guard<std::mutex> lock(op->getLock());
      op->eval(ws);
    } else {
      op->eval(ws);
    }
  }

  auto &parent_clv = ws->clv_ref(_parent_clv_index);
  auto &lchild_clv = ws->clv(_lbranch_clv_index);
  auto &rchild_clv = ws->clv(_rbranch_clv_index);

  weighted_combine(lchild_clv, rchild_clv, _excl_dists, parent_clv,
                   ws->clv_scalar(_lbranch_clv_index),
                   ws->clv_scalar(_rbranch_clv_index),
                   ws->clv_scalar(_parent_clv_index));

  _last_execution = ws->advance_clock();
  ws->update_clv_clock(_parent_clv_index);
}

void SplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                 std::ostream &os, size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "SplitOperation:\n";

  os << tabs << "Lbranch clv (index: " << _lbranch_clv_index
     << ", scalar: " << ws->clv_scalar(_lbranch_clv_index)
     << ", update: " << ws->last_update_clv(_lbranch_clv_index)
     << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_lbranch_clv_index));
  os << tabs << "Rbranch clv (index: " << _rbranch_clv_index
     << ", scalar: " << ws->clv_scalar(_rbranch_clv_index)
     << ", update: " << ws->last_update_clv(_rbranch_clv_index)
     << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_rbranch_clv_index));
  os << tabs << "Parent clv (index: " << _parent_clv_index
     << ", scalar: " << ws->clv_scalar(_parent_clv_index)
     << ", update: " << ws->last_update_clv(_parent_clv_index)
     << "):  " << std::setprecision(10)
     << blaze::trans(ws->clv(_parent_clv_index));
  os << tabs << "Last Executed: " << _last_execution << "\n";
  os << tabs << "correct: " << std::boolalpha << correct(ws) << "\n";

  if (_excl_dists.size() != 0) {
    os << tabs << "Excluded dists:\n";
    for (size_t i = 0; i < _excl_dists.size(); ++i) {
      os << tabs << _excl_dists[i] << "\n";
    }
  }

  os << tabs << "Left Branch ops:\n";
  if (_lbranch_ops.size() != 0) {
    for (auto &op : _lbranch_ops) { os << op->printStatus(ws, tabLevel + 1); }
  }
  os << "\n" << tabs << "Right Branch ops:\n";
  if (_rbranch_ops.size() != 0) {
    for (auto &op : _rbranch_ops) { os << op->printStatus(ws, tabLevel + 1); }
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

void ReverseSplitOperation::eval(std::shared_ptr<Workspace> ws) {
  for (auto &op : _branch_ops) {
    if (op.use_count() > 1) {
      std::lock_guard<std::mutex> lock(op->getLock());
      op->eval(ws);
    } else {
      op->eval(ws);
    }
  }

  if (_eval_clvs) {
    lagrange_col_vector_t tmp(ws->states());
    auto &ltop_clv = ws->clv(_ltop_clv_index);
    auto &rtop_clv = ws->clv(_rtop_clv_index);

    reverse_weighted_combine(ltop_clv, rtop_clv, _excl_dists, tmp);

    ws->update_clv(std::move(tmp), _bot_clv_index);
  }

  _last_execution = ws->advance_clock();
}

void ReverseSplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                        std::ostream &os,
                                        size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "ReverseSplitOperation:\n";

  os << tabs << "Bot clv (index: " << _bot_clv_index
     << ", update: " << ws->last_update_clv(_bot_clv_index)
     << "): " << std::setprecision(10) << blaze::trans(ws->clv(_bot_clv_index));
  os << tabs << "Ltop clv (index: " << _ltop_clv_index
     << ", update: " << ws->last_update_clv(_ltop_clv_index)
     << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_ltop_clv_index));
  os << tabs << "Rtop clv (index: " << _rtop_clv_index
     << ", update: " << ws->last_update_clv(_rtop_clv_index)
     << "): " << std::setprecision(10)
     << blaze::trans(ws->clv(_rtop_clv_index));

  if (_excl_dists.size() != 0) {
    os << tabs << "Excluded dists:\n";
    for (size_t i = 0; i < _excl_dists.size(); ++i) {
      os << tabs << _excl_dists[i] << "\n";
    }
  }

  os << tabs << "Eval CLVS: " << _eval_clvs << "\n";

  if (_branch_ops.size() != 0) {
    os << tabs << "Branch ops:\n";
    for (auto &op : _branch_ops) { os << op->printStatus(ws, tabLevel + 1); }
    os << "\n";
  }
  os << closing_line(tabs);
}

std::string ReverseSplitOperation::printStatus(
    const std::shared_ptr<Workspace> &ws, size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void LLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  _result = std::log(blaze::dot(ws->clv(_root_clv_index),
                                ws->get_base_frequencies(_prior_index))) -
            lagrange_scaling_factor_log * ws->clv_scalar(_root_clv_index);
}

void StateLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  lagrange_col_vector_t tmp(ws->states());
  tmp = 0.0;
  size_t tmp_scalar = 0;
  weighted_combine(ws->clv(_lchild_clv_index), ws->clv(_rchild_clv_index),
                   /*excl_dists=*/{}, tmp, ws->clv_scalar(_lchild_clv_index),
                   ws->clv_scalar(_rchild_clv_index), tmp_scalar);

  tmp_scalar += ws->clv_scalar(_parent_clv_index);

  _result = blaze::log(ws->clv(_parent_clv_index) * tmp) -
            tmp_scalar * lagrange_scaling_factor_log;
}

void SplitLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  std::unordered_map<lagrange_dist_t, std::vector<AncSplit>> ret;

  auto &parent_clv = ws->clv(_parent_clv_index);
  auto &lchild_clv = ws->clv(_lchild_clv_index);
  auto &rchild_clv = ws->clv(_rchild_clv_index);
  std::vector<lagrange_region_split_t> splits;
  for (lagrange_dist_t dist = 0; dist < ws->states(); dist++) {
    std::vector<AncSplit> anc_split_vec;
    generate_splits(dist, ws->regions(), splits);
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
  _result = ret;
}
