/* Operation.cpp
 *
 * Created by Ben Bettisworth
 * 2020-10-27
 */

#include "Operation.hpp"

#include <cstdint>
#include <string_view>

#include "Periods.hpp"
#include "logger.hpp"

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "AncSplit.hpp"
#include "Arnoldi.hpp"
#include "Common.hpp"
#include "Utils.hpp"
#include "Workspace.hpp"

namespace lagrange {

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

/* We need this funtion to make ancestral splits.
 * However, we should make sure that we use the other functions whenever
 * possible.
 * */
void generate_splits(Range state,
                     size_t regions,
                     std::vector<RegionSplit> &results) {
  assert(regions > 0);
  results.clear();
  uint64_t valid_region_mask = (1ULL << regions) - 1;
  if (state == 0) { return; }

  if (lagrange_popcount(state) == 1) {
    results.push_back({state, state});
    return;
  }

  for (size_t i = 0; i < regions; i++) {
    uint64_t x = 1ULL << i;
    if ((state & x) == 0) { continue; }

    results.push_back({x, state});
    results.push_back({state, x});

    uint64_t y = (x ^ state) & valid_region_mask;
    results.push_back({x, y});
    results.push_back({y, x});
  }
}

auto split_product(Range left,
                   Range right,
                   const LagrangeConstColVector &c1,
                   const LagrangeConstColVector &c2,
                   const std::function<size_t(Range)> &dist_map) -> double {
  return c1[dist_map(left)] * c2[dist_map(right)];
}

auto fused_join_splits_happy(Range splitting_range,
                             size_t regions,
                             const LagrangeConstColVector &c1,
                             const LagrangeConstColVector &c2) -> double {
  if (splitting_range == 0) { return 0.0; }

  auto set_regions = lagrange_popcount(splitting_range);

  if (set_regions == 1) { return c1[splitting_range] * c2[splitting_range]; }

  double sum = 0.0;

  /* Grab the values that don't change first */
  const double split_c1 = c1[splitting_range];
  const double split_c2 = c2[splitting_range];

  for (size_t i = 0; i < regions; i++) {
    /* Sympatric split index */
    uint64_t sympatric_index = 1ULL << i;

    /* Allopatric split index */
    uint64_t allopatric_index = (sympatric_index ^ splitting_range);

    double valid = lagrange_bextr(splitting_range, i);

    sum += (c1[sympatric_index] * (split_c2 + c2[allopatric_index])
            + c2[sympatric_index] * (split_c1 * c1[allopatric_index]))
           * valid;
  }

  return sum / static_cast<double>(set_regions * 4);
}

auto fused_join_splits(Range splitting_range,
                       size_t regions,
                       const LagrangeConstColVector &c1,
                       const LagrangeConstColVector &c2,
                       const std::function<size_t(Range)> &dist_map) -> double {
  if (splitting_range == 0) { return 0.0; }
  auto set_regions = lagrange_popcount(splitting_range);

  if (set_regions == 1) {
    return c1[dist_map(splitting_range)] * c2[dist_map(splitting_range)];
  }

  double sum = 0.0;

  auto splitting_index = dist_map(splitting_range);
  const double split_c1 = c1[splitting_index];
  const double split_c2 = c2[splitting_index];

  /* Left splits first */
  for (size_t i = 0; i < regions; i++) {
    if (!lagrange_bextr(splitting_range, i)) { continue; }
    /* Sympatric split index */
    Range singleton_range = 1ULL << i;

    /* Allopatric split index */
    Range allopatric_range = (singleton_range ^ splitting_range);

    auto singleton_index = dist_map(singleton_range);
    auto allopatric_index = dist_map(allopatric_range);

    sum += (c1[singleton_index] * (split_c2 + c2[allopatric_index])
            + c2[singleton_index] * (split_c1 + c1[allopatric_index]));
  }
  return sum / static_cast<double>(set_regions * 4);
}

void join_splits(Range splitting_dist,
                 size_t dist_index,
                 size_t regions,
                 const LagrangeConstColVector &c1,
                 const LagrangeConstColVector &c2,
                 LagrangeColVector dest,
                 bool &scale,
                 const std::function<size_t(Range)> &dist_map) {
  double sum = fused_join_splits(splitting_dist, regions, c1, c2, dist_map);
  scale &= sum < lagrange_scale_threshold;

  assert(std::isfinite(sum));

  dest[dist_index] = sum;
}

void join_splits_happy(Range splitting_dist,
                       size_t dist_index,
                       size_t regions,
                       const LagrangeConstColVector &c1,
                       const LagrangeConstColVector &c2,
                       LagrangeColVector dest,
                       bool &scale) {
  double sum = fused_join_splits_happy(splitting_dist, regions, c1, c2);
  scale &= sum < lagrange_scale_threshold;

  dest[dist_index] = sum;
}

constexpr auto weighted_combine_check_happy_path(
    size_t states,
    size_t max_areas,
    const Option<Range> &fixed_dist,
    Range excl_area_mask,
    Range incl_area_mask) {
  return (states == max_areas) && (!fixed_dist.hasValue())
         && (excl_area_mask == 0) && (incl_area_mask == 0);
}

void weighted_combine_happy(const LagrangeConstColVector &c1,
                            const LagrangeConstColVector &c2,
                            size_t states,
                            size_t regions,
                            LagrangeColVector dest,
                            bool &scale) {
  for (size_t i = 0; i < states; i++) {
    join_splits_happy(i, i, regions, c1, c2, dest, scale);
  }
  for (size_t i = 0; i < states; i++) { assert(std::isfinite(dest[i])); }
}

void weighted_combine(const LagrangeConstColVector &c1,
                      const LagrangeConstColVector &c2,
                      size_t states,
                      size_t regions,
                      size_t max_areas,
                      LagrangeColVector dest,
                      size_t c1_scale,
                      size_t c2_scale,
                      size_t &scale_count,
                      const Option<Range> &fixed_dist,
                      Range excl_area_mask,
                      Range incl_area_mask) {
  assert(states != 0);

  bool scale = true;

  if (weighted_combine_check_happy_path(
          states, max_areas, fixed_dist, excl_area_mask, incl_area_mask)) {
    weighted_combine_happy(c1, c2, states, regions, dest, scale);
  } else {
    Range dist = 0;
    size_t index = 0;
    const auto dist_map = invert_dist_map(regions, max_areas);
    const auto dist_map_func = [&dist_map](Range d) -> size_t {
      return dist_map.at(d);
    };

    while (true) {
      auto next_dist_index =
          next_dist(dist, max_areas, index, excl_area_mask, incl_area_mask);
      dist = next_dist_index.first;
      index = next_dist_index.second;

      if (index >= states) { break; }

      if (fixed_dist.hasValue() && fixed_dist.get() != dist) { continue; }

      join_splits(dist, index, regions, c1, c2, dest, scale, dist_map_func);
    }
  }

  scale_count = c1_scale + c2_scale;

  if (scale) {
    for (size_t i = 0; i < states; i++) { dest[i] *= lagrange_scaling_factor; }

    scale_count += 1;
  }
}

void fused_reverse_join_splits(Range splitting_range,
                               size_t regions,
                               const LagrangeConstColVector &c1,
                               const LagrangeConstColVector &c2,
                               LagrangeColVector dest,
                               const std::function<size_t(Range)> &dist_map) {
  if (splitting_range == 0) {
    dest[dist_map(splitting_range)] = 0.0;
    return;
  }
  auto set_regions = lagrange_popcount(splitting_range);

  if (set_regions == 1) {
    dest[dist_map(splitting_range)] +=
        c1[dist_map(splitting_range)] * c2[dist_map(splitting_range)];
    return;
  }

  const double split_c1 = c1[dist_map(splitting_range)];
  const double weight = 1.0 / (4 * set_regions);

  /* Left splits first */
  for (size_t i = 0; i < regions; i++) {
    /* Sympatric split index */
    Range singleton_range = 1ULL << i;

    Range allopatric_range = (singleton_range ^ splitting_range);

    bool valid = lagrange_bextr(splitting_range, i);
    if (!valid) { continue; }

    auto splitting_index = dist_map(splitting_range);
    auto singleton_index = dist_map(singleton_range);
    auto allopatric_index = dist_map(allopatric_range);

    dest[splitting_index] += split_c1 * c2[singleton_index] * weight;
    dest[allopatric_index] += split_c1 * c2[singleton_index] * weight;
    dest[singleton_index] +=
        (c2[allopatric_index] + c2[splitting_index]) * split_c1 * weight;
  }
}

void reverse_weighted_combine(const LagrangeConstColVector &c1,
                              const LagrangeConstColVector &c2,
                              size_t states,
                              size_t regions,
                              size_t max_areas,
                              LagrangeColVector dest,
                              size_t c1_scale,
                              size_t c2_scale,
                              size_t &scale_count,
                              const Option<Range> &fixed_dist,
                              Range excl_area_mask,
                              Range incl_area_mask) {
  assert(states != 0);
  bool scale = true;

  std::vector<RegionSplit> splits;

  if (max_areas == regions && !fixed_dist.hasValue()) {
    const auto identity_func = [](Range d) -> size_t { return d; };
    for (size_t i = 0; i < states; i++) {
      fused_reverse_join_splits(i, regions, c1, c2, dest, identity_func);
      /*
      reverse_join_splits(
          i, states, regions, splits, c1, c2, scale, dest, identity_func);
      */
    }
  } else {
    const auto dist_map = invert_dist_map(regions, max_areas);
    const auto dist_map_func = [&dist_map](Range d) -> size_t {
      return dist_map.at(d);
    };

    Range dist = 0;
    size_t index = 0;

    while (true) {
      auto next_dist_index =
          next_dist(dist, max_areas, index, excl_area_mask, incl_area_mask);
      dist = next_dist_index.first;
      index = next_dist_index.second;

      if (index >= states) { break; }

      if (fixed_dist.hasValue() && fixed_dist.get() != dist) { continue; }

      fused_reverse_join_splits(dist, regions, c1, c2, dest, dist_map_func);
      /*reverse_join_splits(*/
      /*    dist, states, regions, splits, c1, c2, scale, dest,
       * dist_map_func);*/
    }
  }

  for (size_t i = 0; i < states; ++i) {
    assert(std::isfinite(dest[i]));
    scale &= dest[i] < lagrange_scale_threshold;
  }

  scale_count = c1_scale + c2_scale;
  if (scale) {
    for (size_t i = 0; i < states; ++i) { dest[i] *= lagrange_scaling_factor; }

    scale_count += 1;
  }
}

auto make_tabs(size_t tabLevel) -> std::string {
  std::ostringstream tabs;
  tabs << "|";
  for (size_t i = 0; i < tabLevel; ++i) { tabs << " |"; }
  tabs << " ";
  return tabs.str();
}

auto boarder_line(const std::string &tabs, const std::string &corner_char)
    -> std::string {
  std::ostringstream line;

  line << tabs.substr(0, tabs.size() - 2) << corner_char;
  size_t cur_len = line.str().size();
  for (size_t i = 0; i < (80 - cur_len); ++i) { line << "─"; }

  return line.str();
}

auto opening_line(const std::string &tabs) -> std::string {
  return boarder_line(tabs, "┌");
}

auto closing_line(const std::string &tabs) -> std::string {
  return boarder_line(tabs, "└");
}

void MakeRateMatrixOperation::eval(const std::shared_ptr<Workspace> &ws) {
  auto &rm = ws->rateMatrix(_rate_matrix_index);

  for (size_t i = 0; i < ws->matrixSize(); i++) { rm[i] = 0.0; }

  const auto &period = ws->getPeriodParams(_period_index);

  size_t source_index = 0;
  Range source_dist = 0;

  for (source_index = 0, source_dist = 0;
       source_index < ws->restrictedStateCount();
       source_dist =
           next_dist(source_dist, static_cast<uint32_t>(ws->maxAreas())),
      ++source_index) {
    size_t dest_index = 0;
    Range dest_dist = 0;

    for (dest_index = 0, dest_dist = 0; dest_index < ws->restrictedStateCount();
         dest_dist =
             next_dist(dest_dist, static_cast<uint32_t>(ws->maxAreas())),
        ++dest_index) {
      if (lagrange_popcount(source_dist ^ dest_dist) != 1) { continue; }

      /* Source is "gaining" a region, so we add
       */
      if (source_dist < dest_dist && source_dist != 0) {
        double sum = 0.0;
        size_t i = lagrange_fast_log2(source_dist ^ dest_dist);
        for (size_t j = 0; j < ws->regions(); ++j) {
          sum += (lagrange_bextr(source_dist, j) != 0U)
                     ? period.getDispersionRate(i, j)
                     : 0.0;
        }

        rm[ws->computeMatrixIndex(source_index, dest_index)] = sum;

      }
      /* Otherwise, source is loosing a region, so
         we just set the value */
      else if (source_dist != 0) {
        rm[ws->computeMatrixIndex(source_index, dest_index)] =
            period.getExtinctionRate();
      }
    }
  }

  for (size_t i = 0; i < ws->restrictedStateCount(); i++) {
    double sum = 0.0;
    for (size_t j = 0; j < ws->restrictedStateCount(); j++) {
      sum += rm[ws->computeMatrixIndex(i, j)];
    }
    assert(sum >= 0);
    rm[ws->computeMatrixIndex(i, i)] = -sum;
  }

  _last_execution = ws->advanceClock();
  ws->updateRateMatrixAndAdvanceClock(_rate_matrix_index);
}

void MakeRateMatrixOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                          std::ostream &os,
                                          size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "MakeRateMatrixOperation:\n"
     << tabs << "Rate Matrix (index: " << _rate_matrix_index
     << " update: " << ws->lastUpdateRateMatrix(_rate_matrix_index) << "):\n";

  const auto &rm = ws->rateMatrix(_rate_matrix_index);

  for (size_t i = 0; i < ws->restrictedStateCount(); ++i) {
    os << tabs << std::setprecision(10)
       << std::make_tuple(rm + ws->computeMatrixIndex(i, 0),
                          ws->restrictedStateCount())
       << "\n";
  }

  os << tabs << "Period index: " << _period_index << "\n"
     << tabs << ws->getPeriodParams(_period_index).toString() << "\n"
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
  if ((_rate_matrix_op != nullptr)
      && (_last_execution > _rate_matrix_op->lastUpdate())) {
    return;
  }

  if (_last_execution == 0) {
    _A.reset(new LagrangeMatrixBase[ws->matrixSize()]);
    _lapack_work_buffer.reset(
        new LagrangeMatrixBase[ws->restrictedStateCount()]);
  }

#ifdef MKL_ENABLED
  int ret = 0;
#endif  // MKL_ENABLED
  int rows = static_cast<int>(ws->matrixRows());
  int leading_dim = static_cast<int>(ws->leadingDimension());

  assert(rows > 0);
  assert(leading_dim > 0);

  if (_t == 0.0) {
    _N.reset(new LagrangeMatrixBase[ws->matrixSize()]);

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < rows; ++j) {
        _N[ws->computeMatrixIndex(i, j)] = i == j ? 1.0 : 0.0;
      }
    }

    ws->updateProbMatrix(_prob_matrix_index, _N.get());

    _last_execution = ws->advanceClock();
    _execution_count += 1;

    return;
  }

  for (size_t i = 0; i < ws->matrixSize(); i++) {
    _A.get()[i] = ws->rateMatrix(_rate_matrix_index)[i];
  }

  // We place an arbitrary limit on the size of
  // scale exp because if it is too large we run
  // into numerical issues.
#ifdef MKL_ENABLED
  double inf_norm =
      LAPACKE_dlange(CblasRowMajor, 'I', rows, rows, _A.get(), leading_dim);
#else
  double inf_norm = LAPACK_dlange(
      "I", &rows, &rows, _A.get(), &leading_dim, _lapack_work_buffer.get());
#endif
  assert(inf_norm > 0.0);
  int At_norm = static_cast<int>(inf_norm * _t);
  int scale_exp = std::min(30, std::max(0, 1 + At_norm));

  double Ascal = _t / std::pow(2.0, scale_exp);
  assert(Ascal > 0.0);
  cblas_dscal(ws->matrixSize(), Ascal, _A.get(), 1);

  // q is a magic parameter that controls the
  // number of iterations of the loop higher is
  // more accurate, with each increase of q
  // decreasing error by 4 orders of magnitude.
  // Anything above 12 is probably snake oil.
  constexpr int q = 3;
  double c = 0.5;
  double sign = -1.0;

  if (_last_execution == 0) {
    _X_1.reset(new LagrangeMatrixBase[ws->matrixSize()]);
    _X_2.reset(new LagrangeMatrixBase[ws->matrixSize()]);
    _N.reset(new LagrangeMatrixBase[ws->matrixSize()]);
    _D.reset(new LagrangeMatrixBase[ws->matrixSize()]);
  }

  cblas_dcopy(ws->matrixSize(), _A.get(), 1, _X_1.get(), 1);

  for (size_t i = 0; i < ws->matrixSize(); i++) {
    _X_2.get()[i] = 0.0;
    _N.get()[i] = 0.0;
    _D.get()[i] = 0.0;
  }

  for (size_t i = 0; i < static_cast<size_t>(rows); i++) {
    _N.get()[ws->computeMatrixIndex(i, i)] = 1.0;
    _D.get()[ws->computeMatrixIndex(i, i)] = 1.0;
  }

  cblas_daxpy(ws->matrixSize(), c, _X_1.get(), 1, _N.get(), 1);
  cblas_daxpy(ws->matrixSize(), sign * c, _X_1.get(), 1, _D.get(), 1);

  // Using fortran indexing, and we started an
  // iteration ahead to skip some setup.
  // Furthermore, we are going to unroll the loop
  // to allow us to skip some assignments.
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

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                rows,
                rows,
                rows,
                1.0,
                _A.get(),
                leading_dim,
                _X_1.get(),
                leading_dim,
                0.0,
                _X_2.get(),
                leading_dim);
    cblas_daxpy(ws->matrixSize(), c, _X_2.get(), 1, _N.get(), 1);
    cblas_daxpy(ws->matrixSize(), sign * c, _X_2.get(), 1, _D.get(), 1);

    i += 1;

    if (i > q) { break; }

    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                rows,
                rows,
                rows,
                1.0,
                _A.get(),
                leading_dim,
                _X_2.get(),
                leading_dim,
                0.0,
                _X_1.get(),
                leading_dim);
    cblas_daxpy(ws->matrixSize(), c, _X_1.get(), 1, _N.get(), 1);
    cblas_daxpy(ws->matrixSize(), sign * c, _X_1.get(), 1, _D.get(), 1);

    i += 1;
  }

  {
#ifdef MKL_ENABLED
    auto *ipiv = (long long int *)malloc(sizeof(long long int)
                                         * static_cast<size_t>(rows));
#else
    int *ipiv = (int *)malloc(sizeof(int) * static_cast<size_t>(rows));
#endif
    assert(ipiv != nullptr);

#ifdef MKL_ENABLED
    ret = LAPACKE_dgesv(CblasRowMajor,
                        rows,
                        rows,
                        _D.get(),
                        leading_dim,
                        ipiv,
                        _N.get(),
                        leading_dim);
    assert(ret == 0);
#else
    int info = 0;
    LAPACK_dgesv(&rows,
                 &rows,
                 _D.get(),
                 &leading_dim,
                 ipiv,
                 _N.get(),
                 &leading_dim,
                 &info);
    assert(info == 0);
#endif

    free(ipiv);
  }

  auto &r1 = _N;
  auto &r2 = _D;
  for (int i = 0; i < scale_exp; ++i) {
    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                rows,
                rows,
                rows,
                1.0,
                r1.get(),
                leading_dim,
                r1.get(),
                leading_dim,
                0.0,
                r2.get(),
                leading_dim);
    std::swap(r1, r2);
  }

  if (_transposed) {
#ifdef MKL_ENABLED
    mkl_dimatcopy(
        'r', 't', rows, rows, 1.0, r1.get(), leading_dim, leading_dim);
    for (size_t i = 0; i < static_cast<size_t>(rows); i++) {
      r1.get()[ws->computeMatrixIndex(i, 0)] = 0.0;
    }
    r1.get()[0] = 1.0;
#else
    cblas_domatcopy(CblasRowMajor,
                    CblasTrans,
                    rows,
                    rows,
                    1.0,
                    r1.get(),
                    leading_dim,
                    r2.get(),
                    leading_dim);
    for (int i = 0; i < rows; i++) {
      r2.get()[ws->computeMatrixIndex(i, 0)] = 0.0;
    }
    r2.get()[0] = 1.0;
    std::swap(r1, r2);
#endif
  }

  ws->updateProbMatrix(_prob_matrix_index, r1.get());

  _last_execution = ws->advanceClock();
  _execution_count += 1;
}

void ExpmOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                std::ostream &os,
                                size_t tabLevel) const {
  if (_arnoldi_mode) { return; }
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "ExpmOperation:\n"
     << tabs << "Rate Matrix (index: " << _rate_matrix_index << "):\n";

  os << tabs << "Prob Matrix (index: " << _prob_matrix_index
     << " update: " << ws->lastUpdateProbMatrix(_prob_matrix_index) << "):\n";

  const auto &pm = ws->probMatrix(_prob_matrix_index);

  for (size_t i = 0; i < ws->restrictedStateCount(); ++i) {
    os << tabs << std::setprecision(10)
       << std::make_tuple(pm + ws->computeMatrixIndex(i, 0),
                          ws->restrictedStateCount())
       << "\n";
  }

  os << tabs << "t: " << std::setprecision(16) << _t << "\n";
  os << tabs << "_last_execution: " << _last_execution;
  if (_rate_matrix_op != nullptr) {
    os << "\n" << _rate_matrix_op->printStatus(ws, tabLevel + 1) << "\n";
  } else {
    auto &rm = ws->rateMatrix(_rate_matrix_index);
    for (size_t i = 0; i < ws->restrictedStateCount(); ++i) {
      os << tabs << std::setprecision(10)
         << std::make_tuple(rm + ws->computeMatrixIndex(i, 0),
                            ws->restrictedStateCount())
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
  if (evaluated(ws)) { return; }
  if (_child_op != nullptr && _child_op->ready(ws, _last_execution)) {
    _child_op->eval(ws);
  }

  if (_expm_op != nullptr) {
    if (_expm_op->isArnoldiMode()) {
      expm::multiply_arnoldi_chebyshev(ws,
                                       _expm_op->rateMatrix(),
                                       _bot_clv,
                                       _top_clv,
                                       _expm_op->transposed(),
                                       _expm_op->getT());
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
    cblas_dgemv(CblasRowMajor,
                CblasNoTrans,
                ws->restrictedStateCount(),
                ws->restrictedStateCount(),
                1.0,
                ws->probMatrix(_prob_matrix_index),
                ws->leadingDimension(),
                ws->CLV(_bot_clv),
                1,
                0.0,
                ws->CLV(_top_clv),
                1);
  }

  if (_expm_op->isArnoldiMode() && _expm_op->isAdaptive()) {
    const auto &top_clv = ws->CLV(_top_clv);
    for (size_t i = 0; i < ws->CLVSize(); ++i) {
      if (std::isnan(top_clv[i])
          || ((top_clv[i] > 1.0) != (top_clv[i] < 0.0))) {
        fallback();
        eval(ws);
        break;
      }
    }
  }

  _last_execution = ws->advanceClock();

  ws->CLVScalar(_top_clv) = ws->CLVScalar(_bot_clv);
  ws->updateCLVClock(_top_clv);
}

void DispersionOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                      std::ostream &os,
                                      size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);

  os << opening_line(tabs) << "\n";
  os << tabs << "DispersionOperation:\n";
  os << tabs << "Top clv (index: " << _top_clv
     << " update: " << ws->lastUpdateCLV(_top_clv)
     << "): " << std::setprecision(10)
     << std::make_tuple(ws->CLV(_top_clv), ws->restrictedStateCount()) << "\n";
  os << tabs << "Bot clv (index: " << _bot_clv
     << " update: " << ws->lastUpdateCLV(_bot_clv)
     << "): " << std::setprecision(10)
     << std::make_tuple(ws->CLV(_bot_clv), ws->restrictedStateCount()) << "\n";
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

void DispersionOperation::printGraph(std::ostream &os, size_t &index) const {
  auto ptr = (size_t)(static_cast<const void *>(this));
  os << std::format(R"("{}" [label = "Dispersion {}"];)", ptr, index) << "\n";
  os << std::format(R"("{}" -> clv{};)", ptr, _top_clv) << "\n";
  os << std::format(R"(clv{} -> {};)", _bot_clv, ptr) << "\n";

  index++;

  if (_child_op != nullptr) { _child_op->printGraph(os, index); }
}

static void eval_branch_op(std::shared_ptr<DispersionOperation> &op,
                           const std::shared_ptr<Workspace> &ws) {
  if (!op) { return; }
  if (op.use_count() > 1) {
    std::lock_guard<std::mutex> lock(op->getLock());
    op->eval(ws);
  } else {
    op->eval(ws);
  }
}

void SplitOperation::eval(const std::shared_ptr<Workspace> &ws) {
  eval_branch_op(_lbranch_op, ws);
  eval_branch_op(_rbranch_op, ws);

  const auto &parent_clv = ws->CLV(_parent_clv_index);
  const auto &lchild_clv = ws->CLV(_lbranch_clv_index);
  const auto &rchild_clv = ws->CLV(_rbranch_clv_index);

  weighted_combine(lchild_clv,
                   rchild_clv,
                   ws->restrictedStateCount(),
                   ws->regions(),
                   ws->maxAreas(),
                   parent_clv,
                   ws->CLVScalar(_lbranch_clv_index),
                   ws->CLVScalar(_rbranch_clv_index),
                   ws->CLVScalar(_parent_clv_index),
                   _fixed_dist,
                   _excl_area_mask.get(0),
                   _incl_area_mask.get(0));

  _last_execution = ws->advanceClock();
  ws->updateCLVClock(_parent_clv_index);
}

void SplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                 std::ostream &os,
                                 size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "SplitOperation:\n";

  os << tabs << "Lbranch clv (index: " << _lbranch_clv_index
     << ", scalar: " << ws->CLVScalar(_lbranch_clv_index)
     << ", update: " << ws->lastUpdateCLV(_lbranch_clv_index)
     << "): " << std::setprecision(10) << ws->CLVSizeTuple(_lbranch_clv_index)
     << "\n";
  os << tabs << "Rbranch clv (index: " << _rbranch_clv_index
     << ", scalar: " << ws->CLVScalar(_rbranch_clv_index)
     << ", update: " << ws->lastUpdateCLV(_rbranch_clv_index)
     << "): " << std::setprecision(10) << ws->CLVSizeTuple(_rbranch_clv_index)
     << "\n";
  os << tabs << "Parent clv (index: " << _parent_clv_index
     << ", scalar: " << ws->CLVScalar(_parent_clv_index)
     << ", update: " << ws->lastUpdateCLV(_parent_clv_index)
     << "):  " << std::setprecision(10) << ws->CLVSizeTuple(_parent_clv_index)
     << "\n";
  os << tabs << "Last Executed: " << _last_execution << "\n";

  os << tabs << "Left Branch ops:\n";
  if (_lbranch_op != nullptr) {
    os << _lbranch_op->printStatus(ws, tabLevel + 1);
  }
  os << "\n" << tabs << "Right Branch ops:\n";
  if (_rbranch_op != nullptr) {
    os << _rbranch_op->printStatus(ws, tabLevel + 1);
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

void SplitOperation::printGraph(std::ostream &os, size_t &index) const {
  os << std::format(R"({} [label = "Split Operation {}"];)", index, index)
     << "\n";
  os << std::format("{} -> clv{};\n", index, _parent_clv_index);
  os << std::format("clv{} -> {};\n", _rbranch_clv_index, index);
  os << std::format("clv{} -> {};\n", _lbranch_clv_index, index);

  index++;

  _rbranch_op->printGraph(os, index);
  _lbranch_op->printGraph(os, index);
}

void ReverseSplitOperation::eval(const std::shared_ptr<Workspace> &ws) {
  eval_branch_op(_branch_op, ws);

  if (_eval_clvs) {
    /* TODO: This is a trick to avoid double evaluation of the same operation.
     * For some reason, there are identical operations. The solution is to just
     * check if the op has been done, and then skip it. But The better solution
     * is to only produce one operations. */
    if (_last_execution < ws->lastUpdateCLV(_bot_clv_index)) { return; }

    const auto &ltop_clv = ws->CLV(_ltop_clv_index);
    const auto &rtop_clv = ws->CLV(_rtop_clv_index);

    reverse_weighted_combine(ltop_clv,
                             rtop_clv,
                             ws->restrictedStateCount(),
                             ws->regions(),
                             ws->maxAreas(),
                             ws->CLV(_bot_clv_index),
                             ws->CLVScalar(_ltop_clv_index),
                             ws->CLVScalar(_rtop_clv_index),
                             ws->CLVScalar(_bot_clv_index),
                             _fixed_dist,
                             _excl_area_mask.get(0),
                             _incl_area_mask.get(0));

    bool valid = false;
    for (size_t i = 0; i < ws->CLVSize(); ++i) {
      auto val = ws->CLV(_bot_clv_index)[i];
      assert(!std::isnan(val));
      if (std::isfinite(val) && val != 0.0) { valid = true; }
    }
    if (!valid) {
      LOG_ERROR(
          "Reverse split operation produced and invalid result. Bot clv index: "
          "{}",
          _bot_clv_index);
    }
    assert(valid);

    _last_execution = ws->advanceClock();
    ws->updateCLVClock(_bot_clv_index);
  }
}

void ReverseSplitOperation::printStatus(const std::shared_ptr<Workspace> &ws,
                                        std::ostream &os,
                                        size_t tabLevel) const {
  std::string tabs = make_tabs(tabLevel);
  os << opening_line(tabs) << "\n";
  os << tabs << "ReverseSplitOperation:\n";

  os << tabs << "Bot clv (index: " << _bot_clv_index
     << ", update: " << ws->lastUpdateCLV(_bot_clv_index)
     << "): " << std::setprecision(10) << ws->CLVSizeTuple(_bot_clv_index)
     << "\n";
  os << tabs << "Ltop clv (index: " << _ltop_clv_index
     << ", update: " << ws->lastUpdateCLV(_ltop_clv_index)
     << "): " << std::setprecision(10) << ws->CLVSizeTuple(_ltop_clv_index)
     << "\n";
  os << tabs << "Rtop clv (index: " << _rtop_clv_index
     << ", update: " << ws->lastUpdateCLV(_rtop_clv_index)
     << "): " << std::setprecision(10) << ws->CLVSizeTuple(_rtop_clv_index)
     << "\n";

  /*
  if (!_excl_dists.empty()) {
    os << tabs << "Excluded dists:\n";
    for (unsigned long _excl_dist : _excl_dists) {
      os << tabs << _excl_dist << "\n";
    }
  }
  */

  os << tabs << "Eval CLVS: " << _eval_clvs << "\n";

  if (_branch_op != nullptr) {
    os << tabs << "Branch ops:\n";
    if (_branch_op) { os << _branch_op->printStatus(ws, tabLevel + 1); }
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

void ReverseSplitOperation::printGraph(std::ostream &os, size_t &index) const {
  os << std::format(
      R"({} [label = "Reverse Split Operation {}"];)", index, index)
     << "\n";

  os << std::format("{} -> clv{};\n", index, _bot_clv_index);
  os << std::format("clv{} -> {};\n", _ltop_clv_index, index);
  os << std::format("clv{} -> {};\n", _rtop_clv_index, index);

  index++;

  if (_branch_op != nullptr) { _branch_op->printGraph(os, index); }
}

void LLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  double rho = cblas_ddot(ws->CLVSize(),
                          ws->CLV(_root_clv_index),
                          1,
                          ws->getBaseFrequencies(_prior_index),
                          1);
  assert(rho > 0.0);
  _result = std::log(rho)
            - lagrange_scaling_factor_log * ws->CLVScalar(_root_clv_index);

  _last_execution = ws->advanceClock();
}

auto LLHGoal::ready(const std::shared_ptr<Workspace> &ws) const -> bool {
  return ws->lastUpdateCLV(_root_clv_index) > _last_execution;
}

void LLHGoal::printGraph(std::ostream &os, size_t &index) const {
  os << std::format(R"({} [label = "LLHGoal {}"];)", index, index) << "\n";
  os << std::format("clv{} -> {};\n", _root_clv_index, index);
  os << std::format("bf{} -> {};\n", _prior_index, index);

  index++;
}

void StateLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  if (_last_execution == 0) {
    _result.reset(new LagrangeMatrixBase[ws->restrictedStateCount()]);
    _states = ws->restrictedStateCount();
    for (Range i = 0; i < _states; ++i) { _result[i] = 0.0; }
  }

  size_t tmp_scalar = 0;

  weighted_combine(ws->CLV(_lchild_clv_index),
                   ws->CLV(_rchild_clv_index),
                   _states,
                   ws->regions(),
                   ws->maxAreas(),
                   _result.get(),
                   ws->CLVScalar(_lchild_clv_index),
                   ws->CLVScalar(_rchild_clv_index),
                   tmp_scalar,
                   _fixed_dist,
                   _excl_area_mask,
                   _incl_area_mask);

  tmp_scalar += ws->CLVScalar(_parent_clv_index);
  auto *result = _result.get();

  for (size_t i = 0; i < ws->restrictedStateCount(); ++i) {
    double tmp_val = result[i];
    double parent_val = ws->CLV(_parent_clv_index)[i];
    result[i] = std::log(tmp_val * parent_val)
                - tmp_scalar * lagrange_scaling_factor_log;
    assert(!std::isnan(result[i]));
  }
  _last_execution = ws->advanceClock();
}

auto StateLHGoal::ready(const std::shared_ptr<Workspace> &ws) const -> bool {
  return ws->lastUpdateCLV(_lchild_clv_index) > _last_execution
         && ws->lastUpdateCLV(_rchild_clv_index) > _last_execution
         && ws->lastUpdateCLV(_parent_clv_index) > _last_execution;
}

void StateLHGoal::printGraph(std::ostream &os, size_t &index) const {
  os << std::format(R"({} [label = "State Goal {}"];)", index, index) << "\n";

  os << std::format("clv{} -> {};\n", _parent_clv_index, index);
  os << std::format("clv{} -> {};\n", _lchild_clv_index, index);
  os << std::format("clv{} -> {};\n", _rchild_clv_index, index);

  index++;
}

void SplitLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  std::unordered_map<Range, std::vector<AncSplit>> ret;

  const auto &parent_clv = ws->CLV(_parent_clv_index);
  const auto &lchild_clv = ws->CLV(_lchild_clv_index);
  const auto &rchild_clv = ws->CLV(_rchild_clv_index);

  std::vector<RegionSplit> splits;

  for (Range index = 0; index < ws->restrictedStateCount(); index++) {
    if (_fixed_dist.hasValue() && _fixed_dist.get() != index) { continue; }

    if (!check_excl_dist(index, _excl_area_mask)) { continue; }
    if (!check_incl_dist(index, _incl_area_mask)) { continue; }

    std::vector<AncSplit> anc_split_vec;
    generate_splits(index, ws->regions(), splits);
    double weight = 1.0 / splits.size();
    for (auto sp : splits) {
      AncSplit anc_split(index, sp.left, sp.right, weight);
      double lh = parent_clv[index] * lchild_clv[sp.left] * rchild_clv[sp.right]
                  * weight;
      double loglh = std::log(lh);
      assert(!std::isnan(loglh));
      anc_split.setLikelihood(std::log(lh));
      anc_split_vec.push_back(anc_split);
    }
    ret[index] = anc_split_vec;
  }
  _result = ret;
}

auto SplitLHGoal::ready(const std::shared_ptr<Workspace> &ws) const -> bool {
  return ws->lastUpdateCLV(_lchild_clv_index) > _last_execution
         && ws->lastUpdateCLV(_rchild_clv_index) > _last_execution
         && ws->lastUpdateCLV(_parent_clv_index) > _last_execution;
}
}  // namespace lagrange
