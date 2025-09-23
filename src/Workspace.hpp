/* Workspace.h
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#ifndef LAGRANGE_WORKSPACE_H
#define LAGRANGE_WORKSPACE_H

#include <cassert>
#include <cstddef>
#include <expected>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>

#include "Common.hpp"
#include "Periods.hpp"
#include "Utils.hpp"

namespace lagrange {

enum class SetCLVStatus { definite, ambiguous };
enum class SetCLVError { failed };
SetCLVStatus &operator&=(SetCLVStatus &lhs, const SetCLVStatus &rhs);

using SetCLVResult = std::expected<SetCLVStatus, SetCLVError>;

struct NodeReservation {
  NodeReservation() :
      _top_clv{std::numeric_limits<size_t>::max()},
      _bot1_clv{std::numeric_limits<size_t>::max()},
      _bot2_clv{std::numeric_limits<size_t>::max()},
      _top_rclv{std::numeric_limits<size_t>::max()},
      _bot1_rclv{std::numeric_limits<size_t>::max()},
      _bot2_rclv{std::numeric_limits<size_t>::max()} {}

  size_t _top_clv;
  size_t _bot1_clv;
  size_t _bot2_clv;

  size_t _top_rclv;
  size_t _bot1_rclv;
  size_t _bot2_rclv;
};

struct MatrixReservation {
  LagrangeMatrix _matrix = nullptr;
  ClockTick _last_update = 0;
  OpID _op_id{};
};

struct CLVReservation {
  LagrangeColVector _clv = nullptr;
  ClockTick _last_update = 0;
};

class Workspace {
 public:
  Workspace(size_t taxa_count,
            size_t inner_count,
            size_t regions,
            size_t max_areas) :
      _taxa_count{taxa_count},
      _inner_count{inner_count},
      _regions{regions},
      _states{1ULL << regions},
      _max_areas{max_areas},
      _next_free_clv{0},
      _rate_matrix{1},
      _prob_matrix{1},
      _base_frequencies_count{1},
      _base_frequencies{nullptr},
      _clv_scalars{nullptr},
      _node_reservations{nodeCount()},
      _periods{1},
      _current_clock{0},
      _reserved{false} {
    if (taxa_count == 0) {
      throw std::runtime_error{"We cannot make a workspace with zero taxa"};
    }
    if (_max_areas > _regions) {
      throw std::runtime_error{
          "Max areas cannot be larger than the number of regions"};
    }
    if (_max_areas == _regions) {
      _restricted_state_count = _states;
      _leading_dim = _states;
    } else {
      _restricted_state_count =
          lagrange_compute_restricted_state_count(_regions, _max_areas);
      _leading_dim = _restricted_state_count;
    }
  }

  Workspace(size_t taxa_count, size_t regions, size_t max_areas) :
      Workspace(taxa_count, taxa_count - 1, regions, max_areas) {}

  ~Workspace();

  auto CLVScalar(size_t i) -> size_t & {
    assert(_reserved);
    return _clv_scalars[i];
  }

  auto rateMatrix(size_t i) -> LagrangeMatrix & {
    assert(_reserved);
    if (i >= _rate_matrix.size()) {
      throw std::runtime_error{"Rate matrix access out of range"};
    }
    return _rate_matrix[i]._matrix;
  }

  void updateRateMatrix(size_t index, const LagrangeConstMatrix &A) {
    assert(_reserved);
    if (index >= _rate_matrix.size()) {
      throw std::runtime_error{"Rate matrix access out of range when updating"};
    }

    for (size_t i = 0; i < matrixSize(); i++) {
      _rate_matrix[index]._matrix[i] = A[i];
    }

    _rate_matrix[index]._last_update = advanceClock();
  }

  void updateRateMatrixAndAdvanceClock(size_t i) {
    assert(_reserved);
    _rate_matrix[i]._last_update = advanceClock();
  }

  auto probMatrix(size_t i) -> LagrangeMatrix & {
    assert(_reserved);
    if (i >= _prob_matrix.size()) {
      throw std::runtime_error{"Prob matrix access out of range"};
    }
    return _prob_matrix[i]._matrix;
  }

  auto updateProbMatrix(size_t index, const LagrangeConstMatrix &A)
      -> ClockTick {
    if (index >= _prob_matrix.size()) {
      throw std::runtime_error{"Prob matrix access out of range when updating"};
    }

    for (size_t i = 0; i < matrixSize(); i++) {
      _prob_matrix[index]._matrix[i] = A[i];
    }

    _expm_count += 1;

    return _prob_matrix[index]._last_update = advanceClock();
  }

  auto lastUpdateProbMatrix(size_t i) -> ClockTick {
    assert(_reserved);
    return _prob_matrix[i]._last_update;
  }

  auto lastUpdateRateMatrix(size_t i) -> ClockTick {
    assert(_reserved);
    return _rate_matrix[i]._last_update;
  }

  auto CLV(size_t i) -> const LagrangeColVector & {
    assert(_reserved);
    if (i >= CLVCount()) {
      throw std::runtime_error{"CLV access out of range"};
    }
    return _clvs[i]._clv;
  }

  void updateCLV(const LagrangeConstColVector &clv, size_t index) {
    if (index >= CLVCount()) {
      throw std::runtime_error{"CLV access out of range"};
    }

    for (size_t i = 0; i < restrictedStateCount(); i++) {
      _clvs[index]._clv[i] = clv[i];
    }

    _clvs[index]._last_update = advanceClock();
  }

  auto CLVSizeTuple(size_t index) -> std::tuple<LagrangeColVector, size_t> {
    return std::make_tuple(CLV(index), CLVSize());
  }

  void updateCLVClock(size_t index) {
    _clvs[index]._last_update = advanceClock();
  }

  auto lastUpdateCLV(size_t index) -> ClockTick {
    return _clvs[index]._last_update;
  }

  [[nodiscard]] auto states() const -> size_t { return _states; }

  [[nodiscard]] auto regions() const -> size_t { return _regions; }

  [[nodiscard]] auto probMatrixCount() const -> size_t {
    return _prob_matrix.size();
  }

  [[nodiscard]] auto rateMatrixCount() const -> size_t {
    return _rate_matrix.size();
  }

  [[nodiscard]] auto periodCount() const -> size_t { return _periods.size(); }

  [[nodiscard]] auto CLVCount() const -> size_t { return _next_free_clv; }

  [[nodiscard]] auto matrixSize() const -> size_t {
    return leadingDimension() * restrictedStateCount();
  }

  [[nodiscard]] auto nodeCount() const -> size_t {
    return _inner_count + _taxa_count;
  }

  auto suggestProbMatrixIndex() -> size_t {
    size_t suggested_index = _prob_matrix.size();
    _prob_matrix.emplace_back();
    return suggested_index;
  }

  auto reserveRateMatrixIndex(size_t index) -> size_t {
    if (_rate_matrix.size() <= index) { _rate_matrix.resize(index + 1); }
    return index;
  }

  static auto suggestFreqVectorIndex() -> size_t { return 0; }

  auto registerGenericCLV() -> size_t { return registerCLV(); }

  void registerTopCLV(size_t node_id);

  void registerChildrenCLV(size_t node_id);

  SetCLVResult setTipCLVAmbigious(size_t index, Range clv_index);
  SetCLVResult setTipCLV(size_t index, Range dist);

  void registerTopCLVReverse(size_t node_id) {
    _node_reservations[node_id]._top_rclv = registerCLV();
  }

  void registerBot1CLVReverse(size_t node_id) {
    _node_reservations[node_id]._bot1_rclv = registerCLV();
  }

  void registerBot2CLVReverse(size_t node_id) {
    _node_reservations[node_id]._bot2_rclv = registerCLV();
  }

  auto getTopCLV(size_t node_id) -> size_t {
    return _node_reservations[node_id]._top_clv;
  }

  auto getTopCLVReverse(size_t node_id) -> size_t {
    return _node_reservations[node_id]._top_rclv;
  }

  auto getLeftChildCLV(size_t node_id) -> size_t {
    return _node_reservations[node_id]._bot1_clv;
  }

  auto getRightChildCLV(size_t node_id) -> size_t {
    return _node_reservations[node_id]._bot2_clv;
  }

  auto getBot1CLVReverse(size_t node_id) -> size_t {
    return _node_reservations[node_id]._bot1_rclv;
  }

  auto getBot2CLVReverse(size_t node_id) -> size_t {
    return _node_reservations[node_id]._bot2_rclv;
  }

  auto getBaseFrequencies(size_t index) -> LagrangeColVector & {
    return _base_frequencies[index];
  }

  void setBaseFrequenciesByDist(size_t index, Range dist) {
    for (size_t i = 0; i < CLVSize(); ++i) {
      _base_frequencies[index][i] = 0.0;
    }

    size_t tmp_index = 0;
    Range tmp_dist = 0;
    while (true) {
      tmp_dist = next_dist(tmp_dist, maxAreas());
      tmp_index += 1;
      if (dist == tmp_dist) { break; }
    }

    _base_frequencies[index][tmp_index] = 1.0;
  }

  void reserve();

  auto advanceClock() -> ClockTick { return _current_clock++; }

  auto readClock() -> ClockTick { return _current_clock; }

  [[nodiscard]] auto getPeriodParams() const
      -> const std::vector<PeriodParams> & {
    return _periods;
  }

  [[nodiscard]] auto getPeriodParams() -> std::vector<PeriodParams> & {
    return _periods;
  }

  [[nodiscard]] auto getPeriodParams(size_t period_index) -> PeriodParams & {
    return _periods[period_index];
  }

  [[nodiscard]] auto getPeriodParams(size_t period_index) const
      -> const PeriodParams & {
    return _periods[period_index];
  }

  void initPeriodParams(const std::vector<PeriodParams> &p) { _periods = p; }

  void setPeriodParams(size_t period_index, const PeriodParams &p);

  void setPeriodParamsCount(size_t periods) { _periods.resize(periods); }

  void setPeriodAdjustmentMatrix(size_t period_index,
                                 std::shared_ptr<double[]>);

  [[nodiscard]] auto reserved() const -> bool {
    return _base_frequencies != nullptr && !_clvs.empty();
  }

  [[nodiscard]] auto reportNodeVectors(size_t node_id) const -> std::string;

  [[nodiscard]] auto computeMatrixIndex(size_t i, size_t j) const -> size_t {
    return i * leadingDimension() + j;
  }

  void setReversePrior(size_t index) {
    assert(_reserved);
    if (index >= _clvs.size()) {
      throw std::runtime_error{
          "Attempting to access a clv that doesn't exist when setting a "
          "reverse prior"};
    }

    for (size_t i = 0; i < restrictedStateCount(); i++) {
      _clvs[index]._clv[i] = 1.0;
    }

    _clvs[index]._last_update = std::numeric_limits<ClockTick>::max();
  }

  [[nodiscard]] auto leadingDimension() const -> size_t { return _leading_dim; }

  [[nodiscard]] auto restrictedStateCount() const -> size_t {
    return _restricted_state_count;
  }

  [[nodiscard]] auto maxAreas() const -> size_t { return _max_areas; }

  [[nodiscard]] auto matrixRows() const -> size_t {
    return restrictedStateCount();
  }

  [[nodiscard]] auto matrixCols() const -> size_t { return matrixRows(); }

  [[nodiscard]] auto CLVSize() const -> size_t {
    return restrictedStateCount();
  }

  [[nodiscard]] auto EXPMExecutionCount() const -> size_t {
    return _expm_count;
  }

  [[nodiscard]] auto getOptDimsSize() const -> size_t {
    size_t count = 0;
    for (auto p : _periods) { count += p.paramCount(); }
    return count;
  }

  auto makeParameterVector(double rate, double dist) const
      -> std::vector<double> {
    std::vector<double> ret(getOptDimsSize());
    size_t i = 0;
    for (auto p : _periods) {
      ret[i++] = rate;
      ret[i++] = rate;
      if (p.hasAdjustmentMatrix()) { ret[i++] = dist; }
    }
    return ret;
  }

  [[nodiscard]] auto getDefaultParams() const -> std::vector<double> {
    constexpr double rate_default = 0.01;
    constexpr double penalty_default = -1.0;
    return makeParameterVector(rate_default, penalty_default);
  }

  [[nodiscard]] auto getParamMins() const -> std::vector<double> {
    constexpr double rate_min = 1e-7;
    constexpr double penalty_min = -10;
    return makeParameterVector(rate_min, penalty_min);
  }

  [[nodiscard]] auto getParamMaxs() const -> std::vector<double> {
    constexpr double rate_max = 1e2;
    constexpr double penalty_max = 0.0;
    return makeParameterVector(rate_max, penalty_max);
  }

  void setAdjustmentMatrix(size_t index, const std::shared_ptr<double[]> adj) {
    _periods[index].adjustment_matrix = adj;
  }

  template <typename R>
  void setAdjustmentMatrices(const R &mats)
    requires(requires(R r) {
      requires std::ranges::range<R>;
      { *r.begin() } -> std::same_as<std::shared_ptr<double[]>>;
      { r.size() };
    })
  {
    lagrange_assert(_periods.size() == mats.size(),
                    "Setting adjusment matrices failed");
#ifdef __cpp_lib_ranges_zip
    for (auto [a, b] : std::views::zip(_periods, mats)) {
      a.adjustment_matrix = b;
    }
#else
    for (size_t i = 0; i < _periods.size(); ++i) {
      _periods[i].adjustment_matrix = mats[i];
    }
#endif
  }

  template <typename R>
  void setPeriodRangeCounts(const R &ranges)
    requires(requires(R r) {
      requires std::ranges::range<R>;
      { *r.begin() } -> std::same_as<size_t>;
      { r.size() };
    })
  {
    lagrange_assert(_periods.size() == ranges.size(),
                    "Setting adjusment matrices failed");
#ifdef __cpp_lib_ranges_zip
    for (auto [a, b] : std::views::zip(_periods, ranges)) { a.regions = b; }
#else
    for (size_t i = 0; i < _periods.size(); ++i) {
      _periods[i].regions = ranges[i];
    }
#endif
  }

  inline size_t expectedSplitResultSize() const {
    auto n = regions();
    size_t split_count = n * (n + 1) * (2 * n + 1) / 6;
    return split_count * CLVSize() * sizeof(double) * _inner_count;
  }

  inline size_t expectedStateResultSize() const {
    return CLVSize() * sizeof(double) * _inner_count;
  }

 private:
  auto registerCLV() -> size_t { return _next_free_clv++; }

  size_t _taxa_count;
  size_t _inner_count;
  size_t _regions;
  size_t _states;
  size_t _max_areas;
  size_t _restricted_state_count;
  size_t _next_free_clv;

  size_t _leading_dim;

  std::vector<MatrixReservation> _rate_matrix;
  std::vector<MatrixReservation> _prob_matrix;

  size_t _base_frequencies_count;
  LagrangeColVector *_base_frequencies;

  std::vector<CLVReservation> _clvs;
  size_t *_clv_scalars;

  std::vector<NodeReservation> _node_reservations;

  std::vector<PeriodParams> _periods;

  Clock _current_clock;
  bool _reserved;

  size_t _expm_count = 0;
};

}  // namespace lagrange
#endif
