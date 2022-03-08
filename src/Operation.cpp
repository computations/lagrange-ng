/* Operation.cpp
 *
 * Created by Ben Bettisworth
 * 2020-10-27
 */

#include <array>
#include <cmath>
#include <cstddef>
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
#include "Common.h"
#include "Operation.h"
#include "Utils.h"
#include "Workspace.h"

inline lagrange_dist_t next_dist(lagrange_dist_t d, uint32_t n) {
  d += 1;
  while (static_cast<size_t>(__builtin_popcountll(d)) > n) { d++; }
  return d;
}

std::ostream &operator<<(std::ostream &os,
                         std::tuple<double *, size_t> vector_tuple) {
  double *v;
  size_t n;
  std::tie(v, n) = vector_tuple;
  os << "(";
  for (size_t i = 0; i < n; i++) {
    os << v[i];
    if (i != n - 1) { os << ", "; }
  }
  os << ")";
  return os;
}

inline void generate_splits(uint64_t state, size_t regions,
                            std::vector<lagrange_region_split_t> &results) {
  results.clear();
  uint64_t valid_region_mask = (1ull << regions) - 1;
  if (state == 0) { return; }

  if (lagrange_popcount(state) == 1) {
    results.push_back({state, state});
    return;
  }

  // results.reserve(regions);
  for (size_t i = 0; i < regions; i++) {
    uint64_t x = 1ull << i;
    if ((state & x) == 0) { continue; }

    results.push_back({x, state});
    results.push_back({state, x});

    uint64_t y = (x ^ state) & valid_region_mask;
    results.push_back({x, y});
    if (lagrange_popcount(y) > 1) { results.push_back({y, x}); }
  }
}

inline void weighted_combine(const lagrange_const_col_vector_t &c1,
                             const lagrange_const_col_vector_t &c2,
                             size_t states, size_t max_areas,
                             lagrange_col_vector_t dest, size_t c1_scale,
                             size_t c2_scale, size_t &scale_count) {
  size_t regions = lagrange_fast_log2(states);

  bool scale = true;

  std::vector<lagrange_region_split_t> splits;

  if (max_areas == states) {
    for (size_t i = 0; i < states; i++) {
      generate_splits(i, regions, splits);
      double sum = 0.0;
      for (auto &p : splits) { sum += c1[p.left] * c2[p.right]; }

      if (splits.size() != 0) { sum /= static_cast<double>(splits.size()); }

      scale &= sum < lagrange_scale_threshold;

      dest[i] = sum;
    }
  } else {
    for (lagrange_dist_t i = 0; i < states; i = next_dist(i, max_areas)) {
      generate_splits(i, regions, splits);
      double sum = 0.0;
      for (auto &p : splits) { sum += c1[p.left] * c2[p.right]; }

      if (splits.size() != 0) { sum /= static_cast<double>(splits.size()); }

      scale &= sum < lagrange_scale_threshold;

      dest[i] = sum;
    }
  }

  scale_count = c1_scale + c2_scale;

  if (scale) {
    for (size_t i = 0; i < states; i++) { dest[i] *= lagrange_scaling_factor; }

    scale_count += 1;
  }
}

inline void reverse_weighted_combine(const lagrange_const_col_vector_t &c1,
                                     const lagrange_const_col_vector_t &c2,
                                     size_t states, size_t max_areas,
                                     lagrange_col_vector_t dest) {
  size_t regions = lagrange_fast_log2(states);

  std::vector<lagrange_region_split_t> splits;

  if (max_areas == regions) {
    for (size_t i = 0; i < states; i++) {
      generate_splits(i, regions, splits);
      if (splits.size() == 0) { continue; }

      double weight = 1.0 / static_cast<double>(splits.size());
      for (auto &p : splits) { dest[p.left] += c1[i] * c2[p.right] * weight; }
    }
  } else {
    for (lagrange_dist_t i = 0; i < states; i = next_dist(i, max_areas)) {
      generate_splits(i, regions, splits);
      if (splits.size() == 0) { continue; }

      double weight = 1.0 / static_cast<double>(splits.size());
      for (auto &p : splits) { dest[p.left] += c1[i] * c2[p.right] * weight; }
    }
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

void MakeRateMatrixOperation::eval(const std::shared_ptr<Workspace> &ws) {
  auto &rm = ws->rate_matrix(_rate_matrix_index);

  for (size_t i = 0; i < ws->matrix_size(); i++) { rm[i] = 0.0; }

  auto &period = ws->get_period_params(_period_index);

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

std::string MakeRateMatrixOperation::printStatus(
    const std::shared_ptr<Workspace> &ws, size_t tabLevel) const {
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

  int rows = static_cast<int>(ws->matrix_rows());
  int leading_dim = static_cast<int>(ws->leading_dimension());

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
  int At_norm = static_cast<int>(inf_norm * _t);
  int scale_exp = std::min(30, std::max(0, 1 + At_norm));

  double Ascal = _t / std::pow(2.0, scale_exp);
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

    /*
    bli_gemm_ex(&BLIS_ONE, &_A, &_X_1, &BLIS_ZERO, &_X_2, blis_context,
                blis_runtime);
    bli_axpym_ex(&c_obj, &_X_2, &_N, blis_context, blis_runtime);
    bli_axpym_ex(&signed_c_obj, &_X_2, &_D, blis_context, blis_runtime);
    */

    i += 1;

    if (i > q) { break; }

    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows,
                1.0, _A.get(), leading_dim, _X_2.get(), leading_dim, 0.0,
                _X_1.get(), leading_dim);
    cblas_daxpy(ws->matrix_size(), c, _X_1.get(), 1, _N.get(), 1);
    cblas_daxpy(ws->matrix_size(), sign * c, _X_1.get(), 1, _D.get(), 1);

    /*
    bli_gemm_ex(&BLIS_ONE, &_A, &_X_2, &BLIS_ZERO, &_X_1, blis_context,
                blis_runtime);
    bli_axpym_ex(&c_obj, &_X_1, &_N, blis_context, blis_runtime);
    bli_axpym_ex(&signed_c_obj, &_X_1, &_D, blis_context, blis_runtime);
    */

    i += 1;
  }

  {
    int *ipiv = (int *)malloc(sizeof(int) * static_cast<size_t>(rows));

#ifdef MKL_ENABLED
    LAPACKE_dgesv(CblasRowMajor, rows, rows, _D.get(), leading_dim, ipiv,
                  _N.get(), leading_dim);
#else
    int info = 0;
    LAPACK_dgesv(&rows, &rows, _D.get(), &leading_dim, ipiv, _N.get(),
                 &leading_dim, &info);
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

std::string ExpmOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                       size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void DispersionOperation::eval(const std::shared_ptr<Workspace> &ws) {
  if (_expm_op != nullptr) {
    if (_expm_op.use_count() > 1) {
      std::lock_guard<std::mutex> expm_guard(_expm_op->getLock());
      _expm_op->eval(ws);
    } else {
      _expm_op->eval(ws);
    }
  }

  if (ws->last_update_clv(_bot_clv) < ws->last_update_clv(_top_clv)) {
    /* we have already computed this operation, return */
    return;
  }

  cblas_dgemv(CblasRowMajor, CblasNoTrans, ws->restricted_state_count(),
              ws->restricted_state_count(), 1.0,
              ws->prob_matrix(_prob_matrix_index), ws->leading_dimension(),
              ws->clv(_bot_clv), 1, 0.0, ws->clv(_top_clv), 1);

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

std::string DispersionOperation::printStatus(
    const std::shared_ptr<Workspace> &ws, size_t tabLevel) const {
  std::ostringstream os;
  printStatus(ws, os, tabLevel);
  return os.str();
}

void SplitOperation::eval(const std::shared_ptr<Workspace> &ws) {
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

  const auto &parent_clv = ws->clv(_parent_clv_index);
  const auto &lchild_clv = ws->clv(_lbranch_clv_index);
  const auto &rchild_clv = ws->clv(_rbranch_clv_index);

  weighted_combine(
      lchild_clv, rchild_clv, ws->restricted_state_count(), ws->max_areas(),
      parent_clv, ws->clv_scalar(_lbranch_clv_index),
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

void ReverseSplitOperation::eval(const std::shared_ptr<Workspace> &ws) {
  for (auto &op : _branch_ops) {
    if (op.use_count() > 1) {
      std::lock_guard<std::mutex> lock(op->getLock());
      op->eval(ws);
    } else {
      op->eval(ws);
    }
  }

  if (_eval_clvs) {
    auto &ltop_clv = ws->clv(_ltop_clv_index);
    auto &rtop_clv = ws->clv(_rtop_clv_index);

    reverse_weighted_combine(ltop_clv, rtop_clv, ws->restricted_state_count(),
                             ws->max_areas(), ws->clv(_bot_clv_index));

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
  double rho =
      cblas_ddot(ws->restricted_state_count(), ws->clv(_root_clv_index), 1,
                 ws->get_base_frequencies(_prior_index), 1);
  _result = std::log(rho) -
            lagrange_scaling_factor_log * ws->clv_scalar(_root_clv_index);

  _last_execution = ws->advance_clock();
}

bool LLHGoal::ready(const std::shared_ptr<Workspace> &ws) const {
  return ws->last_update_clv(_root_clv_index) > _last_execution;
}

void StateLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  if (_last_execution == 0) {
    _result.reset(new lagrange_matrix_base_t[ws->restricted_state_count()]);
    _states = ws->restricted_state_count();
  }

  size_t tmp_scalar = 0;

  weighted_combine(ws->clv(_lchild_clv_index), ws->clv(_rchild_clv_index),
                   ws->restricted_state_count(), ws->max_areas(), _result.get(),
                   ws->clv_scalar(_lchild_clv_index),
                   ws->clv_scalar(_rchild_clv_index), tmp_scalar);

  tmp_scalar += ws->clv_scalar(_parent_clv_index);

  for (size_t i = 0; i < ws->restricted_state_count(); ++i) {
    double tmp_val = _result.get()[i];
    double parent_val = ws->clv(_parent_clv_index)[i];
    _result.get()[i] = std::log(tmp_val * parent_val) -
                       tmp_scalar * lagrange_scaling_factor_log;
  }
  _last_execution = ws->advance_clock();
}

bool StateLHGoal::ready(const std::shared_ptr<Workspace> &ws) const {
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

bool SplitLHGoal::ready(const std::shared_ptr<Workspace> &ws) const {
  return ws->last_update_clv(_lchild_clv_index) > _last_execution &&
         ws->last_update_clv(_rchild_clv_index) > _last_execution &&
         ws->last_update_clv(_parent_clv_index) > _last_execution;
}
