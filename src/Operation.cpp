/* Operation.cpp
 *
 * Created by Ben Bettisworth
 * 2020-10-27
 */

#undef NDEBUG
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "AncSplit.h"
#include "Arnoldi.h"
#include "Common.h"
#include "Operation.h"
#include "Utils.h"
#include "Workspace.h"

static auto operator<<(std::ostream &os,
                       std::tuple<double *, size_t> vector_tuple)
    -> std::ostream & {
  double *v = nullptr;
  size_t n = 0;
  std::tie(v, n) = vector_tuple;
  os << "(";
  for (size_t i = 0; i < n; i++) {
    os << v[i];
    if (i != n - 1) { os << ", "; }
  }
  os << ")";
  return os;
}

inline void generate_splits(uint64_t state, size_t regions, size_t max_areas,
                            std::vector<lagrange_region_split_t> &results) {
  assert(regions > 0);
  results.clear();
  uint64_t valid_region_mask = (1ULL << regions) - 1;
  if (state == 0) { return; }

  if (lagrange_popcount(state) == 1) {
    results.push_back({state, state});
    return;
  }

  // results.reserve(regions);
  for (size_t i = 0; i < regions; i++) {
    uint64_t x = 1ULL << i;
    if ((state & x) == 0) { continue; }

    results.push_back({x, state});
    results.push_back({state, x});

    uint64_t y = (x ^ state) & valid_region_mask;
    results.push_back({x, y});
    if (lagrange_popcount(y) > 1 && lagrange_popcount(y) <= max_areas) {
      results.push_back({y, x});
    }
  }
}

inline void join_splits(
    lagrange_dist_t i, size_t dist_index, size_t regions, size_t max_areas,
    std::vector<lagrange_region_split_t> &splits,
    const lagrange_const_col_vector_t &c1,
    const lagrange_const_col_vector_t &c2, lagrange_col_vector_t dest,
    bool &scale, const std::function<size_t(lagrange_dist_t)> &dist_map) {
  generate_splits(i, regions, max_areas, splits);
  double sum = 0.0;
  for (auto &p : splits) {
    sum += c1[dist_map(p.left)] * c2[dist_map(p.right)];
  }

  if (!splits.empty()) { sum /= static_cast<double>(splits.size()); }

  scale &= sum < lagrange_scale_threshold;

  dest[dist_index] = sum;
}

/* Produces a dist -> index map */
static std::vector<size_t> invert_dist_map(size_t regions, size_t max_areas) {
  size_t states = 1ULL << regions;
  std::vector<size_t> ret(states, std::numeric_limits<size_t>::max());

  size_t index = 0;
  lagrange_dist_t dist = 0;
  for (index = 0, dist = 0; dist < states;
       ++index, dist = next_dist(dist, max_areas)) {
    ret[dist] = index;
  }

  return ret;
}

inline void weighted_combine(const lagrange_const_col_vector_t &c1,
                             const lagrange_const_col_vector_t &c2,
                             size_t states, size_t regions, size_t max_areas,
                             lagrange_col_vector_t dest, size_t c1_scale,
                             size_t c2_scale, size_t &scale_count) {
  assert(states != 0);
  // size_t regions = lagrange_fast_log2(states);

  bool scale = true;

  std::vector<lagrange_region_split_t> splits;

  if (max_areas == states) {
    const auto identity_func = [](lagrange_dist_t d) -> size_t { return d; };
    for (size_t i = 0; i < states; i++) {
      join_splits(i, i, regions, max_areas, splits, c1, c2, dest, scale,
                  identity_func);
    }
  } else {
    lagrange_dist_t dist = 0;
    size_t index = 0;
    const auto dist_map = invert_dist_map(regions, max_areas);
    const auto dist_map_func = [&dist_map](lagrange_dist_t d) -> size_t {
      return dist_map.at(d);
    };
    for (dist = 0, index = 0; index < states;
         dist = next_dist(dist, max_areas), index++) {
      join_splits(dist, index, regions, max_areas, splits, c1, c2, dest, scale,
                  dist_map_func);
    }
  }

  scale_count = c1_scale + c2_scale;

  if (scale) {
    for (size_t i = 0; i < states; i++) { dest[i] *= lagrange_scaling_factor; }

    scale_count += 1;
  }
}

inline void reverse_join_splits(
    lagrange_dist_t i, size_t regions, size_t max_areas,
    std::vector<lagrange_region_split_t> &splits,
    const lagrange_const_col_vector_t &c1,
    const lagrange_const_col_vector_t &c2, lagrange_col_vector_t dest,
    const std::function<size_t(lagrange_dist_t)> &dist_map) {
  generate_splits(i, regions, max_areas, splits);
  if (splits.empty()) { return; }

  const size_t total_states = 1ULL << regions;

  double weight = 1.0 / static_cast<double>(splits.size());
  for (auto &p : splits) {
    size_t l_index = dist_map(p.left);
    size_t r_index = dist_map(p.right);
    size_t i_index = dist_map(i);
    assert(l_index < total_states);
    assert(r_index < total_states);
    assert(i_index < total_states);
    dest[l_index] += c1[i_index] * c2[r_index] * weight;
  }
  // for (auto &p : splits) { dest[p.left] += c1[i] * c2[p.right] * weight; }
}

inline void reverse_weighted_combine(const lagrange_const_col_vector_t &c1,
                                     const lagrange_const_col_vector_t &c2,
                                     size_t states, size_t regions,
                                     size_t max_areas,
                                     lagrange_col_vector_t dest) {
  assert(states != 0);
  // size_t regions = lagrange_fast_log2(states);

  std::vector<lagrange_region_split_t> splits;

  if (max_areas == regions) {
    const auto identity_func = [](lagrange_dist_t d) -> size_t { return d; };
    for (size_t i = 0; i < states; i++) {
      reverse_join_splits(i, regions, max_areas, splits, c1, c2, dest,
                          identity_func);
    }
  } else {
    const auto dist_map = invert_dist_map(regions, max_areas);
    const auto dist_map_func = [&dist_map](lagrange_dist_t d) -> size_t {
      return dist_map.at(d);
    };
    for (lagrange_dist_t i = 0; i < states; i = next_dist(i, max_areas)) {
      reverse_join_splits(i, regions, max_areas, splits, c1, c2, dest,
                          dist_map_func);
    }
  }
}

inline auto make_tabs(size_t tabLevel) -> std::string {
  std::ostringstream tabs;
  tabs << "|";
  for (size_t i = 0; i < tabLevel; ++i) { tabs << " |"; }
  tabs << " ";
  return tabs.str();
}

inline auto boarder_line(const std::string &tabs,
                         const std::string &corner_char) -> std::string {
  std::ostringstream line;

  line << tabs.substr(0, tabs.size() - 2) << corner_char;
  size_t cur_len = line.str().size();
  for (size_t i = 0; i < (80 - cur_len); ++i) { line << "─"; }

  return line.str();
}

inline auto opening_line(const std::string &tabs) -> std::string {
  return boarder_line(tabs, "┌");
}

inline auto closing_line(const std::string &tabs) -> std::string {
  return boarder_line(tabs, "└");
}

void MakeRateMatrixOperation::eval(const std::shared_ptr<Workspace> &ws) {
  auto &rm = ws->rate_matrix(_rate_matrix_index);

  for (size_t i = 0; i < ws->matrix_size(); i++) { rm[i] = 0.0; }

  const auto &period = ws->get_period_params(_period_index);

  size_t source_index = 0;
  lagrange_dist_t source_dist = 0;

  for (source_index = 0, source_dist = 0;
       source_dist < ws->restricted_state_count();
       source_dist = next_dist(source_dist, ws->max_areas()), ++source_index) {
    size_t dest_index = 0;
    lagrange_dist_t dest_dist = 0;

    for (dest_index = 0, dest_dist = 0;
         dest_dist < ws->restricted_state_count();
         dest_dist =
             next_dist(dest_dist, static_cast<uint32_t>(ws->max_areas())),
        ++dest_index) {
      if (lagrange_popcount(source_dist ^ dest_dist) != 1) { continue; }

      /* Source is "gaining" a region, so we add */
      if (source_dist < dest_dist && source_dist != 0) {
        double sum = 0.0;
        size_t i = lagrange_fast_log2(source_dist ^ dest_dist);
        for (size_t j = 0; j < ws->regions(); ++j) {
          sum +=
              period.getDispersionRate(i, j) * lagrange_bextr(source_dist, j);
        }

        rm[ws->compute_matrix_index(source_index, dest_index)] = sum;

      }
      /* Otherwise, source is loosing a region, so we just set the value */
      else if (source_dist != 0) {
        rm[ws->compute_matrix_index(source_index, dest_index)] =
            period.getExtinctionRate();
      }
    }
  }

  for (size_t i = 0; i < ws->restricted_state_count(); i++) {
    double sum = 0.0;
    for (size_t j = 0; j < ws->restricted_state_count(); j++) {
      sum += rm[ws->compute_matrix_index(i, j)];
    }
    assert(sum >= 0);
    rm[ws->compute_matrix_index(i, i)] = -sum;
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

  const auto &rm = ws->rate_matrix(_rate_matrix_index);

  for (size_t i = 0; i < ws->restricted_state_count(); ++i) {
    os << tabs << std::setprecision(10)
       << std::make_tuple(rm + ws->compute_matrix_index(i, 0),
                          ws->restricted_state_count())
       << "\n";
  }

  os << tabs << "Period index: " << _period_index << "\n"
     << tabs << ws->get_period_params(_period_index).toString() << "\n"
     << tabs << "_last_execution: " << _last_execution << "\n"
     << tabs << "_last_update: " << _last_update << "\n";
  os << closing_line(tabs);
}

auto MakeRateMatrixOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                          size_t tabLevel) const
    -> std::string {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void ExpmOperation::eval(const std::shared_ptr<Workspace> &ws) {
  if ((_rate_matrix_op != nullptr) &&
      (_last_execution > _rate_matrix_op->last_update())) {
    return;
  }

  if (_last_execution == 0) {
    _A.reset(new lagrange_matrix_base_t[ws->matrix_size()]);
    _lapack_work_buffer.reset(
        new lagrange_matrix_base_t[ws->restricted_state_count()]);
  }

  int ret = 0;
  int rows = static_cast<int>(ws->matrix_rows());
  int leading_dim = static_cast<int>(ws->leading_dimension());

  assert(rows > 0);
  assert(leading_dim > 0);

  for (size_t i = 0; i < ws->matrix_size(); i++) {
    _A.get()[i] = ws->rate_matrix(_rate_matrix_index)[i];
  }

  // We place an arbitrary limit on the size of scale exp because if it is too
  // large we run into numerical issues.
#ifdef MKL_ENABLED
  double inf_norm =
      LAPACKE_dlange(CblasRowMajor, 'I', rows, rows, _A.get(), leading_dim);
#else
  double inf_norm = LAPACK_dlange("I", &rows, &rows, _A.get(), &leading_dim,
                                  _lapack_work_buffer.get());
#endif
  assert(inf_norm > 0.0);
  int At_norm = static_cast<int>(inf_norm * _t);
  int scale_exp = std::min(30, std::max(0, 1 + At_norm));

  double Ascal = _t / std::pow(2.0, scale_exp);
  assert(Ascal > 0.0);
  cblas_dscal(ws->matrix_size(), Ascal, _A.get(), 1);

  // q is a magic parameter that controls the number of iterations of the loop
  // higher is more accurate, with each increase of q decreasing error by 4
  // orders of magnitude. Anything above 12 is probably snake oil.
  constexpr int q = 3;
  double c = 0.5;
  double sign = -1.0;

  if (_last_execution == 0) {
    _X_1.reset(new lagrange_matrix_base_t[ws->matrix_size()]);
    _X_2.reset(new lagrange_matrix_base_t[ws->matrix_size()]);
    _N.reset(new lagrange_matrix_base_t[ws->matrix_size()]);
    _D.reset(new lagrange_matrix_base_t[ws->matrix_size()]);
  }

  cblas_dcopy(ws->matrix_size(), _A.get(), 1, _X_1.get(), 1);

  for (size_t i = 0; i < ws->matrix_size(); i++) {
    _X_2.get()[i] = 0.0;
    _N.get()[i] = 0.0;
    _D.get()[i] = 0.0;
  }

  for (size_t i = 0; i < static_cast<size_t>(rows); i++) {
    _N.get()[ws->compute_matrix_index(i, i)] = 1.0;
    _D.get()[ws->compute_matrix_index(i, i)] = 1.0;
  }

  cblas_daxpy(ws->matrix_size(), c, _X_1.get(), 1, _N.get(), 1);
  cblas_daxpy(ws->matrix_size(), sign * c, _X_1.get(), 1, _D.get(), 1);

  // Using fortran indexing, and we started an iteration ahead to skip some
  // setup. Furhthermore, we are going to unroll the loop to allow us to skip
  // some assignments.
#if 0
  for (int i = 2; i <= q; i++) {
    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    std::cout << "C: " << c << std::endl;
    X = A * X;
    _N += c * X;
    sign *= -1.0;
    _D += sign * c * X;
  }
#endif

  for (int i = 2; i <= q;) {
    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows,
                1.0, _A.get(), leading_dim, _X_1.get(), leading_dim, 0.0,
                _X_2.get(), leading_dim);
    cblas_daxpy(ws->matrix_size(), c, _X_2.get(), 1, _N.get(), 1);
    cblas_daxpy(ws->matrix_size(), sign * c, _X_2.get(), 1, _D.get(), 1);

    i += 1;

    if (i > q) { break; }

    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows,
                1.0, _A.get(), leading_dim, _X_2.get(), leading_dim, 0.0,
                _X_1.get(), leading_dim);
    cblas_daxpy(ws->matrix_size(), c, _X_1.get(), 1, _N.get(), 1);
    cblas_daxpy(ws->matrix_size(), sign * c, _X_1.get(), 1, _D.get(), 1);

    i += 1;
  }

  {
    int *ipiv = (int *)malloc(sizeof(int) * static_cast<size_t>(rows));
    assert(ipiv != nullptr);

#ifdef MKL_ENABLED
    ret = LAPACKE_dgesv(CblasRowMajor, rows, rows, _D.get(), leading_dim, ipiv,
                        _N.get(), leading_dim);
    assert(ret == 0);
#else
    int info = 0;
    LAPACK_dgesv(&rows, &rows, _D.get(), &leading_dim, ipiv, _N.get(),
                 &leading_dim, &info);
    assert(info == 0);
#endif

    free(ipiv);
  }

  auto &r1 = _N;
  auto &r2 = _D;
  for (int i = 0; i < scale_exp; ++i) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows,
                1.0, r1.get(), leading_dim, r1.get(), leading_dim, 0.0,
                r2.get(), leading_dim);
    std::swap(r1, r2);
  }

  if (_transposed) {
#ifdef MKL_ENABLED
    mkl_dimatcopy('r', 't', rows, rows, 1.0, r1.get(), leading_dim,
                  leading_dim);
    for (size_t i = 0; i < static_cast<size_t>(rows); i++) {
      r1.get()[ws->compute_matrix_index(i, 0)] = 0.0;
    }
    r1.get()[0] = 1.0;
#else
    cblas_domatcopy(CblasRowMajor, CblasTrans, rows, rows, 1.0, r1.get(),
                    leading_dim, r2.get(), leading_dim);
    for (int i = 0; i < rows; i++) {
      r2.get()[ws->compute_matrix_index(i, 0)] = 0.0;
    }
    r2.get()[0] = 1.0;
    std::swap(r1, r2);
#endif
  }

  ws->update_prob_matrix(_prob_matrix_index, r1.get());

  _last_execution = ws->advance_clock();
  _execution_count += 1;
}

void ExpmOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                std::ostream &os, size_t tabLevel) const {
  if (_arnoldi_mode) { return; }
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "ExpmOperation:\n"
     << tabs << "Rate Matrix (index: " << _rate_matrix_index << "):\n";

  os << tabs << "Prob Matrix (index: " << _prob_matrix_index
     << " update: " << ws->last_update_prob_matrix(_prob_matrix_index)
     << "):\n";

  const auto &pm = ws->prob_matrix(_prob_matrix_index);

  for (size_t i = 0; i < ws->restricted_state_count(); ++i) {
    os << tabs << std::setprecision(10)
       << std::make_tuple(pm + ws->compute_matrix_index(i, 0),
                          ws->restricted_state_count())
       << "\n";
  }

  os << tabs << "t: " << std::setprecision(16) << _t << "\n";
  os << tabs << "_last_execution: " << _last_execution;
  if (_rate_matrix_op != nullptr) {
    os << "\n" << _rate_matrix_op->printStatus(ws, tabLevel + 1) << "\n";
  } else {
    auto &rm = ws->rate_matrix(_rate_matrix_index);
    for (size_t i = 0; i < ws->restricted_state_count(); ++i) {
      os << tabs << std::setprecision(10)
         << std::make_tuple(rm + ws->compute_matrix_index(i, 0),
                            ws->restricted_state_count())
         << "\n";
    }
  }
  os << closing_line(tabs);
}

auto ExpmOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                size_t tabLevel) const -> std::string {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void DispersionOperation::eval(const std::shared_ptr<Workspace> &ws) {
  if (ws->last_update_clv(_bot_clv) < ws->last_update_clv(_top_clv)) {
    /* we have already computed this operation, return */
    return;
  }

  if (_expm_op != nullptr) {
    if (_expm_op->isArnoldiMode()) {
      expm::multiply_arnoldi_chebyshev(ws, _expm_op->rate_matrix(), _bot_clv,
                                       _top_clv, _expm_op->transposed(),
                                       _expm_op->get_t());
    } else {
      if (_expm_op.use_count() > 1) {
        std::lock_guard<std::mutex> expm_guard(_expm_op->getLock());
        _expm_op->eval(ws);
      } else {
        _expm_op->eval(ws);
      }
    }
  }

  if (_expm_op != nullptr && !_expm_op->isArnoldiMode()) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, ws->restricted_state_count(),
                ws->restricted_state_count(), 1.0,
                ws->prob_matrix(_prob_matrix_index), ws->leading_dimension(),
                ws->clv(_bot_clv), 1, 0.0, ws->clv(_top_clv), 1);
  }
  _last_execution = ws->advance_clock();

  ws->clv_scalar(_top_clv) = ws->clv_scalar(_bot_clv);
  ws->update_clv_clock(_top_clv);
}

void DispersionOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                      std::ostream &os, size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);

  os << opening_line(tabs) << "\n";
  os << tabs << "DispersionOperation:\n";
  os << tabs << "Top clv (index: " << _top_clv
     << " update: " << ws->last_update_clv(_top_clv)
     << "): " << std::setprecision(10)
     << std::make_tuple(ws->clv(_top_clv), ws->restricted_state_count())
     << "\n";
  os << tabs << "Bot clv (index: " << _bot_clv
     << " update: " << ws->last_update_clv(_bot_clv)
     << "): " << std::setprecision(10)
     << std::make_tuple(ws->clv(_bot_clv), ws->restricted_state_count())
     << "\n";
  os << tabs << "Last Executed: " << _last_execution << "\n";
  os << tabs << "Prob Matrix (index: " << _prob_matrix_index << ")\n";

  if (_expm_op != nullptr) { os << _expm_op->printStatus(ws, tabLevel + 1); }
  os << "\n";
  os << closing_line(tabs);
}

auto DispersionOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                      size_t tabLevel) const -> std::string {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

static inline void eval_branch_ops(
    const std::vector<std::shared_ptr<DispersionOperation>> &branch_ops,
    const std::shared_ptr<Workspace> &ws) {
  for (const auto &op : branch_ops) {
    if (op.use_count() > 1) {
      std::lock_guard<std::mutex> lock(op->getLock());
      op->eval(ws);
    } else {
      op->eval(ws);
    }
  }
}

void SplitOperation::eval(const std::shared_ptr<Workspace> &ws) {
  eval_branch_ops(_lbranch_ops, ws);
  eval_branch_ops(_rbranch_ops, ws);

  const auto &parent_clv = ws->clv(_parent_clv_index);
  const auto &lchild_clv = ws->clv(_lbranch_clv_index);
  const auto &rchild_clv = ws->clv(_rbranch_clv_index);

  weighted_combine(
      lchild_clv, rchild_clv, ws->restricted_state_count(), ws->regions(),
      ws->max_areas(), parent_clv, ws->clv_scalar(_lbranch_clv_index),
      ws->clv_scalar(_rbranch_clv_index), ws->clv_scalar(_parent_clv_index));

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
     << "): " << std::setprecision(10) << ws->clv_size_tuple(_lbranch_clv_index)
     << "\n";
  os << tabs << "Rbranch clv (index: " << _rbranch_clv_index
     << ", scalar: " << ws->clv_scalar(_rbranch_clv_index)
     << ", update: " << ws->last_update_clv(_rbranch_clv_index)
     << "): " << std::setprecision(10) << ws->clv_size_tuple(_rbranch_clv_index)
     << "\n";
  os << tabs << "Parent clv (index: " << _parent_clv_index
     << ", scalar: " << ws->clv_scalar(_parent_clv_index)
     << ", update: " << ws->last_update_clv(_parent_clv_index)
     << "):  " << std::setprecision(10) << ws->clv_size_tuple(_parent_clv_index)
     << "\n";
  os << tabs << "Last Executed: " << _last_execution << "\n";

  os << tabs << "Left Branch ops:\n";
  if (!_lbranch_ops.empty()) {
    for (const auto &op : _lbranch_ops) {
      os << op->printStatus(ws, tabLevel + 1);
    }
  }
  os << "\n" << tabs << "Right Branch ops:\n";
  if (!_rbranch_ops.empty()) {
    for (const auto &op : _rbranch_ops) {
      os << op->printStatus(ws, tabLevel + 1);
    }
  }
  os << "\n";
  os << closing_line(tabs);
}

auto SplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                 size_t tabLevel) const -> std::string {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void ReverseSplitOperation::eval(const std::shared_ptr<Workspace> &ws) {
  eval_branch_ops(_branch_ops, ws);

  if (_eval_clvs) {
    const auto &ltop_clv = ws->clv(_ltop_clv_index);
    const auto &rtop_clv = ws->clv(_rtop_clv_index);

    reverse_weighted_combine(ltop_clv, rtop_clv, ws->restricted_state_count(),
                             ws->regions(), ws->max_areas(),
                             ws->clv(_bot_clv_index));

    for (size_t i = 0; i < ws->clv_size(); ++i) {
      assert(!std::isnan(ws->clv(_bot_clv_index)[i]));
    }

    _last_execution = ws->advance_clock();
    ws->update_clv_clock(_bot_clv_index);
  }
}

void ReverseSplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                        std::ostream &os,
                                        size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "ReverseSplitOperation:\n";

  os << tabs << "Bot clv (index: " << _bot_clv_index
     << ", update: " << ws->last_update_clv(_bot_clv_index)
     << "): " << std::setprecision(10) << ws->clv_size_tuple(_bot_clv_index)
     << "\n";
  os << tabs << "Ltop clv (index: " << _ltop_clv_index
     << ", update: " << ws->last_update_clv(_ltop_clv_index)
     << "): " << std::setprecision(10) << ws->clv_size_tuple(_ltop_clv_index)
     << "\n";
  os << tabs << "Rtop clv (index: " << _rtop_clv_index
     << ", update: " << ws->last_update_clv(_rtop_clv_index)
     << "): " << std::setprecision(10) << ws->clv_size_tuple(_rtop_clv_index)
     << "\n";

  if (!_excl_dists.empty()) {
    os << tabs << "Excluded dists:\n";
    for (unsigned long _excl_dist : _excl_dists) {
      os << tabs << _excl_dist << "\n";
    }
  }

  os << tabs << "Eval CLVS: " << _eval_clvs << "\n";

  if (!_branch_ops.empty()) {
    os << tabs << "Branch ops:\n";
    for (const auto &op : _branch_ops) {
      os << op->printStatus(ws, tabLevel + 1);
    }
    os << "\n";
  }
  os << closing_line(tabs);
}

auto ReverseSplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                        size_t tabLevel) const -> std::string {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void LLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  double rho =
      cblas_ddot(ws->restricted_state_count(), ws->clv(_root_clv_index), 1,
                 ws->get_base_frequencies(_prior_index), 1);
  assert(rho > 0.0);
  _result = std::log(rho) -
            lagrange_scaling_factor_log * ws->clv_scalar(_root_clv_index);

  _last_execution = ws->advance_clock();
}

auto LLHGoal::ready(const std::shared_ptr<Workspace> &ws) const -> bool {
  return ws->last_update_clv(_root_clv_index) > _last_execution;
}

void StateLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  if (_last_execution == 0) {
    _result.reset(new lagrange_matrix_base_t[ws->restricted_state_count()]);
    _states = ws->restricted_state_count();
    for (lagrange_dist_t i = 0; i < _states; ++i) { _result[i] = NAN; }
  }

  size_t tmp_scalar = 0;

  weighted_combine(ws->clv(_lchild_clv_index), ws->clv(_rchild_clv_index),
                   _states, ws->regions(), ws->max_areas(), _result.get(),
                   ws->clv_scalar(_lchild_clv_index),
                   ws->clv_scalar(_rchild_clv_index), tmp_scalar);

  tmp_scalar += ws->clv_scalar(_parent_clv_index);
  auto result = _result.get();

  for (size_t i = 0; i < ws->restricted_state_count(); ++i) {
    double tmp_val = result[i];
    double parent_val = ws->clv(_parent_clv_index)[i];
    result[i] = std::log(tmp_val * parent_val) -
                tmp_scalar * lagrange_scaling_factor_log;
    assert(!std::isnan(result[i]));
  }
  _last_execution = ws->advance_clock();
}

auto StateLHGoal::ready(const std::shared_ptr<Workspace> &ws) const -> bool {
  return ws->last_update_clv(_lchild_clv_index) > _last_execution &&
         ws->last_update_clv(_rchild_clv_index) > _last_execution &&
         ws->last_update_clv(_parent_clv_index) > _last_execution;
}

void SplitLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  std::unordered_map<lagrange_dist_t, std::vector<AncSplit>> ret;

  const auto &parent_clv = ws->clv(_parent_clv_index);
  const auto &lchild_clv = ws->clv(_lchild_clv_index);
  const auto &rchild_clv = ws->clv(_rchild_clv_index);

  std::vector<lagrange_region_split_t> splits;

  for (lagrange_dist_t dist = 0; dist < ws->restricted_state_count(); dist++) {
    std::vector<AncSplit> anc_split_vec;
    generate_splits(dist, ws->regions(), ws->max_areas(), splits);
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

auto SplitLHGoal::ready(const std::shared_ptr<Workspace> &ws) const -> bool {
  return ws->last_update_clv(_lchild_clv_index) > _last_execution &&
         ws->last_update_clv(_rchild_clv_index) > _last_execution &&
         ws->last_update_clv(_parent_clv_index) > _last_execution;
}
