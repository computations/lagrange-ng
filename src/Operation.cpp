/* Operation.cpp
 *
 * Created by Ben Bettisworth
 * 2020-10-27
 */

#include <array>
#include <cstddef>
#include <iomanip>
#include <iostream>
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

  // results.reserve(regions);
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
                             lagrange_col_vector_t dest, size_t c1_scale,
                             size_t c2_scale, size_t &scale_count) {
  size_t states = bli_obj_length(c1);
  size_t regions = lagrange_fast_log2(states);

  size_t idx_excl = 0;

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
    double sum = 0.0;
    for (auto &p : splits) {
      double c1_val = bli_getiv_real(c1, p.left);
      double c2_val = bli_getiv_real(c2, p.right);

      sum += c1_val * c2_val;
    }
    if (splits.size() != 0) { sum /= static_cast<double>(splits.size()); }
    if (sum < lagrange_scale_threshold) {
      scale &= true;
    } else {
      scale &= false;
    }
    bli_setiv_real(dest, sum, i);
  }

  scale_count = c1_scale + c2_scale;
  if (scale) {
    double scal_copy = lagrange_scaling_factor;
    auto scal_factor_obj = bli_obj_wrap_scalar(scal_copy);
    bli_scalv(&scal_factor_obj, dest);
    scale_count += 1;
  }
}

inline void reverse_weighted_combine(
    const lagrange_col_vector_t &c1, const lagrange_col_vector_t &c2,
    const std::vector<lagrange_dist_t> excl_dists, lagrange_col_vector_t dest) {
  size_t states = bli_obj_length(c1);
  size_t regions = lagrange_fast_log2(states);

  size_t idx_excl = 0;

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
    for (auto &p : splits) {
      double c1_val = bli_getiv_real(c1, i);
      double c2_val = bli_getiv_real(c2, p.right);
      double dest_val = bli_getiv_real(dest, p.left);

      dest_val += c1_val * c2_val * weight;

      bli_setiv_real(dest, dest_val, p.left);
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

  bli_setm_scalar(0.0, rm);

  auto &period = ws->get_period_params(_period_index);
  for (lagrange_dist_t dist = 0; dist < ws->states(); dist++) {
    for (lagrange_dist_t i = 0; i < ws->regions(); i++) {
      if (lagrange_bextr(dist, i) != 0) { continue; }
      lagrange_dist_t gain_dist = dist | (1ul << i);

      bli_setijm_real(rm, period.getExtinctionRate(), gain_dist, dist);

      if (dist == 0) { continue; }

      double sum = 0.0;
      for (size_t j = 0; j < ws->regions(); ++j) {
        sum += period.getDispersionRate(j, i) * lagrange_bextr(dist, j);
      }
      bli_setijm_real(rm, sum, dist, gain_dist);
    }
  }

  for (size_t i = 0; i < ws->states(); i++) {
    double sum = 0.0;
    for (size_t j = 0; j < ws->states(); j++) {
      sum += bli_getijm_real(rm, i, j);
    }
    bli_setijm_real(rm, -sum, i, i);
  }

  _last_execution = ws->advance_clock();
  ws->update_rate_matrix_clock(_rate_matrix_index);
}

#if 0
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
#endif

void ExpmOperation::eval(const std::shared_ptr<Workspace> &ws,
                         cntx_t *blis_context, rntm_t *blis_runtime) {
  if ((_rate_matrix_op != nullptr) &&
      (_last_execution > _rate_matrix_op->last_update())) {
    return;
  }

  if (_last_execution == 0) {
    bli_obj_create_conf_to(ws->rate_matrix(_rate_matrix_index), &_A);
  }
  bli_copym(ws->rate_matrix(_rate_matrix_index), &_A);

  // std::cout << "Prob Matrix(t:" << _t << ")" << std::endl;
  // bli_printm("_A", &_A, "%1.10f", "");

  // We place an arbitrary limit on the size of scale exp because if it is too
  // large we run into numerical issues.

  int At_norm = static_cast<int>(bli_inf_normm(&_A) * _t);
  int scale_exp = std::min(30, std::max(0, 1 + At_norm));

  double Ascal = _t / std::pow(2.0, scale_exp);
  auto Ascal_obj = bli_obj_wrap_scalar(Ascal);
  bli_scalm(&Ascal_obj, &_A);

  // q is a magic parameter that controls the number of iterations of the loop
  // higher is more accurate, with each increase of q decreasing error by 4
  // orders of magnitude. Anything above 12 is probably snake oil.
  constexpr int q = 3;
  double c = 0.5;
  double sign = -1.0;
  double signed_c = c * sign;
  // blaze::IdentityMatrix<double, blaze::columnMajor> I(rows);

#if 0
  lagrange_matrix_t X = _A;
  lagrange_matrix_t N = I + c * X;
  lagrange_matrix_t D = I - c * X;
#endif

  if (_last_execution == 0) {
    bli_obj_create_conf_to(&_A, &_X_1);
    bli_obj_create_conf_to(&_A, &_X_2);
    bli_obj_create_conf_to(&_A, &_N);
    bli_obj_create_conf_to(&_A, &_D);
  }

  bli_copym(&_A, &_X_1);

  bli_setm_scalar(0.0, &_X_2);
  bli_setm_scalar(0.0, &_N);
  bli_setm_scalar(0.0, &_D);
  bli_setd(&BLIS_ONE, &_N);
  bli_setd(&BLIS_ONE, &_D);
  // bli_printm("_N post init", &_N, "%1.15f", "");

  obj_t c_obj = bli_obj_wrap_scalar(c);
  obj_t signed_c_obj = bli_obj_wrap_scalar(signed_c);

  bli_axpym_ex(&c_obj, &_X_1, &_N, blis_context, blis_runtime);
  bli_axpym_ex(&signed_c_obj, &_X_1, &_D, blis_context, blis_runtime);

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
    signed_c = c * sign;

    bli_gemm_ex(&BLIS_ONE, &_A, &_X_1, &BLIS_ZERO, &_X_2, blis_context,
                blis_runtime);
    bli_axpym_ex(&c_obj, &_X_2, &_N, blis_context, blis_runtime);
    bli_axpym_ex(&signed_c_obj, &_X_2, &_D, blis_context, blis_runtime);

    i += 1;

    if (i > q) { break; }

    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;
    signed_c = c * sign;

    bli_gemm_ex(&BLIS_ONE, &_A, &_X_2, &BLIS_ZERO, &_X_1, blis_context,
                blis_runtime);
    bli_axpym_ex(&c_obj, &_X_1, &_N, blis_context, blis_runtime);
    bli_axpym_ex(&signed_c_obj, &_X_1, &_D, blis_context, blis_runtime);

    i += 1;
  }

  {
    int n = ws->states();
    // int lwork = n * n;
    int lda = bli_obj_col_stride(&_N);
    int info = 0;
    int *ipiv = (int *)malloc(sizeof(int) * n);

    /*
    double *tau = (double *)malloc(sizeof(double) * n);
    double *workspace = (double *)malloc(sizeof(double) * lwork);

    LAPACK_dgeqrf(&n, &n, static_cast<double *>(bli_obj_buffer(&_D)), &lda, tau,
                  workspace, &lwork, &info);

    LAPACK_dormqr("L", "T", &n, &n, &n,
                  static_cast<double *>(bli_obj_buffer(&_D)), &lda, tau,
                  static_cast<double *>(bli_obj_buffer(&_N)), &lda, workspace,
                  &lwork, &info);

    LAPACK_dtrtrs("U", "N", "N", &n, &n,
                  static_cast<double *>(bli_obj_buffer(&_D)), &lda,
                  static_cast<double *>(bli_obj_buffer(&_N)), &lda, &info);
                 */

    LAPACK_dgesv(&n, &n, static_cast<double *>(bli_obj_buffer(&_D)), &lda, ipiv,
                 static_cast<double *>(bli_obj_buffer(&_N)), &lda, &info);

    // free(tau);
    // free(workspace);
    free(ipiv);
  }

  // bli_printm("final matrix before squaring", &_N, "%1.10f", "");

  auto &r1 = _N;
  auto &r2 = _D;
  for (int i = 0; i < scale_exp; ++i) {
    bli_gemm_ex(&BLIS_ONE, &r1, &r1, &BLIS_ZERO, &r2, blis_context,
                blis_runtime);
    std::swap(r1, r2);
  }

  if (_transposed) {
    bli_obj_set_onlytrans(BLIS_TRANSPOSE, &r1);

    for (size_t i = 1; i < ws->states(); ++i) { bli_setijm_real(&r1, 0, 0, i); }
    bli_setijm_real(&r1, 1, 0, 0);
  }

  // bli_printm("Final Matrix: ", &r1, "%1.10f", "");
  ws->update_prob_matrix(_prob_matrix_index, &r1);

  _last_execution = ws->advance_clock();
}

#if 0
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
#endif

void DispersionOperation::eval(const std::shared_ptr<Workspace> &ws,
                               cntx_t *blis_context, rntm_t *blis_runtime) {
  if (_expm_op != nullptr) {
    if (_expm_op.use_count() > 1) {
      std::lock_guard<std::mutex>(_expm_op->getLock());
      _expm_op->eval(ws, blis_context, blis_runtime);
    } else {
      _expm_op->eval(ws, blis_context, blis_runtime);
    }
  }

  if (ws->last_update_clv(_bot_clv) < ws->last_update_clv(_top_clv)) {
    /* we have already computed stuff, return */
    return;
  }

  /*
  std::cout << "Executing disp op" << std::endl;
  bli_printv("bot clv:", ws->clv(_bot_clv), "%5.5f", "");
  bli_printv("top clv:", ws->clv(_top_clv), "%5.5f", "");
  bli_printm("prob matrix:", ws->prob_matrix(_prob_matrix_index), "%5.5f", "");
  */

  bli_gemv_ex(&BLIS_ONE, ws->prob_matrix(_prob_matrix_index), ws->clv(_bot_clv),
              &BLIS_ZERO, ws->clv(_top_clv), blis_context, blis_runtime);

  /*
  bli_printv("bot clv:", ws->clv(_bot_clv), "%5.5f", "");
  bli_printv("top clv:", ws->clv(_top_clv), "%5.5f", "");
  */

  _last_execution = ws->advance_clock();

  ws->clv_scalar(_top_clv) = ws->clv_scalar(_bot_clv);
  ws->update_clv_clock(_top_clv);
}

#if 0
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
#endif

void SplitOperation::eval(const std::shared_ptr<Workspace> &ws,
                          cntx_t *blis_context, rntm_t *blis_runtime) {
  for (auto &op : _lbranch_ops) {
    if (op.use_count() > 1) {
      std::lock_guard<std::mutex> lock(op->getLock());
      op->eval(ws, blis_context, blis_runtime);
    } else {
      op->eval(ws, blis_context, blis_runtime);
    }
  }
  for (auto &op : _rbranch_ops) {
    if (op.use_count() > 1) {
      std::lock_guard<std::mutex> lock(op->getLock());
      op->eval(ws, blis_context, blis_runtime);
    } else {
      op->eval(ws, blis_context, blis_runtime);
    }
  }

  auto &parent_clv = ws->clv_ref(_parent_clv_index);
  auto &lchild_clv = ws->clv(_lbranch_clv_index);
  auto &rchild_clv = ws->clv(_rbranch_clv_index);

  /*
  bli_printv("parent clv:", parent_clv, "%5.5f", "");
  bli_printv("lchild clv:", lchild_clv, "%5.5f", "");
  bli_printv("rchild clv:", rchild_clv, "%5.5f", "");
  */

  weighted_combine(lchild_clv, rchild_clv, _excl_dists, parent_clv,
                   ws->clv_scalar(_lbranch_clv_index),
                   ws->clv_scalar(_rbranch_clv_index),
                   ws->clv_scalar(_parent_clv_index));

  // bli_printv("parent clv:", parent_clv, "%5.5f", "");
  _last_execution = ws->advance_clock();
  ws->update_clv_clock(_parent_clv_index);
}

#if 0
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
#endif

void ReverseSplitOperation::eval(const std::shared_ptr<Workspace> &ws,
                                 cntx_t *blis_context, rntm_t *blis_runtime) {
  for (auto &op : _branch_ops) {
    if (op.use_count() > 1) {
      std::lock_guard<std::mutex> lock(op->getLock());
      op->eval(ws, blis_context, blis_runtime);
    } else {
      op->eval(ws, blis_context, blis_runtime);
    }
  }

  if (_eval_clvs) {
    auto &ltop_clv = ws->clv(_ltop_clv_index);
    auto &rtop_clv = ws->clv(_rtop_clv_index);

    // bli_printv("ltop_clv", ltop_clv, "%1.20f", "");
    // bli_printv("rtop_clv", rtop_clv, "%1.20f", "");

    reverse_weighted_combine(ltop_clv, rtop_clv, _excl_dists,
                             ws->clv(_bot_clv_index));

    // bli_printv("bot_clv", ws->clv(_bot_clv_index), "%1.20f", "");

    _last_execution = ws->advance_clock();
    ws->update_clv_clock(_bot_clv_index);
  }
}

#if 0
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
#endif

void LLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  double rho = bli_dotv_scalar(ws->clv(_root_clv_index),
                               ws->get_base_frequencies(_prior_index));
  _result = std::log(rho) -
            lagrange_scaling_factor_log * ws->clv_scalar(_root_clv_index);
}

void StateLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  bli_obj_free(_result.get());
  bli_obj_create_conf_to(ws->clv(_lchild_clv_index), _result.get());
  size_t tmp_scalar = 0;

  // bli_printv("lchild clv:", ws->clv(_lchild_clv_index), "%1.15f", "");
  // bli_printv("rchild clv:", ws->clv(_rchild_clv_index), "%1.15f", "");

  weighted_combine(ws->clv(_lchild_clv_index), ws->clv(_rchild_clv_index),
                   /*excl_dists=*/{}, _result.get(),
                   ws->clv_scalar(_lchild_clv_index),
                   ws->clv_scalar(_rchild_clv_index), tmp_scalar);

  tmp_scalar += ws->clv_scalar(_parent_clv_index);

  // bli_printv("parent clv:", ws->clv(_parent_clv_index), "%1.15f", "");

  for (size_t i = 0; i < ws->states(); ++i) {
    double tmp_val = bli_getiv_real(_result.get(), i);
    double parent_val = bli_getiv_real(ws->clv(_parent_clv_index), i);
    bli_setiv_real(_result.get(),
                   std::log(tmp_val * parent_val) -
                       tmp_scalar * lagrange_scaling_factor_log,
                   i);
  }
  // bli_printv("result:", _result.get(), "%1.15f", "");
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
      double lh = bli_getiv_real(parent_clv, dist) *
                  bli_getiv_real(lchild_clv, sp.left) *
                  bli_getiv_real(rchild_clv, sp.right) * weight;
      anc_split.setLikelihood(lh);
      anc_split_vec.push_back(anc_split);
    }
    ret[dist] = anc_split_vec;
  }
  _result = ret;
}
