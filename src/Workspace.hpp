/* Workspace.h
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#ifndef LAGRANGE_WORKSPACE_H
#define LAGRANGE_WORKSPACE_H

#include <cassert>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>

#include "Common.hpp"
#include "Utils.hpp"

namespace lagrange {

struct MatrixReservation {
  LagrangeMatrix _matrix = nullptr;
  ClockTick _last_update = 0;
  OpID _op_id{};
};

struct CLVReservation {
  LagrangeColVector _clv = nullptr;
  ClockTick _last_update = 0;
  uint64_t _state_mask;

  size_t vector_size() const { return lagrange_popcount(_state_mask); }
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
      _leading_dim{_states},
      _rate_matrix_reserved_count{0},
      _prob_matrix_reserved_count{0},
      _rate_matrix{},
      _prob_matrix{},
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
    _restricted_state_count =
        lagrange_compute_restricted_state_count(_regions, _max_areas);
  }

  Workspace(size_t taxa_count, size_t regions, size_t max_areas) :
      Workspace(taxa_count, taxa_count - 1, regions, max_areas) {}

  ~Workspace();

  inline auto CLVScalar(size_t i) -> size_t & {
    assert(_reserved);
    return _clv_scalars[i];
  }

  inline auto rateMatrix(size_t i) -> LagrangeMatrix & {
    assert(_reserved);
    if (i >= _rate_matrix.size()) {
      throw std::runtime_error{"Rate matrix access out of range"};
    }
    return _rate_matrix[i]._matrix;
  }

  inline void updateRateMatrix(size_t index, const LagrangeConstMatrix &A) {
    assert(_reserved);
    if (index >= _rate_matrix.size()) {
      throw std::runtime_error{"Rate matrix access out of range when updating"};
    }

    for (size_t i = 0; i < matrixSize(); i++) {
      _rate_matrix[index]._matrix[i] = A[i];
    }

    _rate_matrix[index]._last_update = advanceClock();
  }

  inline void updateRateMatrixClock(size_t i) {
    assert(_reserved);
    _rate_matrix[i]._last_update = advanceClock();
  }

  inline auto probMatrix(size_t i) -> LagrangeMatrix & {
    assert(_reserved);
    if (i >= _prob_matrix.size()) {
      throw std::runtime_error{"Prob matrix access out of range"};
    }
    return _prob_matrix[i]._matrix;
  }

  inline auto updateProbMatrix(size_t index,
                               const LagrangeConstMatrix &A) -> ClockTick {
    if (index >= _prob_matrix.size()) {
      throw std::runtime_error{"Prob matrix access out of range when updating"};
    }

    for (size_t i = 0; i < matrixSize(); i++) {
      _prob_matrix[index]._matrix[i] = A[i];
    }

    _expm_count += 1;

    return _prob_matrix[index]._last_update = advanceClock();
  }

  inline auto lastUpdateProbMatrix(size_t i) -> ClockTick {
    assert(_reserved);
    return _prob_matrix[i]._last_update;
  }

  inline auto lastUpdateRateMatrix(size_t i) -> ClockTick {
    assert(_reserved);
    return _rate_matrix[i]._last_update;
  }

  inline auto CLV(size_t i) -> const LagrangeColVector & {
    assert(_reserved);
    if (i >= CLVCount()) {
      throw std::runtime_error{"CLV access out of range"};
    }
    return _clvs[i]._clv;
  }

  inline void updateCLV(const LagrangeConstColVector &clv, size_t index) {
    if (index >= CLVCount()) {
      throw std::runtime_error{"CLV access out of range"};
    }

    for (size_t i = 0; i < restrictedStateCount(); i++) {
      _clvs[index]._clv[i] = clv[i];
    }

    _clvs[index]._last_update = advanceClock();
  }

  inline auto CLVSizeTuple(size_t index)
      -> std::tuple<LagrangeColVector, size_t> {
    return std::make_tuple(CLV(index), CLVSize());
  }

  inline void updateCLVClock(size_t index) {
    _clvs[index]._last_update = advanceClock();
  }

  inline auto lastUpdateCLV(size_t index) -> ClockTick {
    return _clvs[index]._last_update;
  }

  inline auto states() const -> size_t { return _states; }

  inline auto regions() const -> size_t { return _regions; }

  inline auto probMatrixCount() const -> size_t { return _prob_matrix.size(); }

  inline auto rateMatrixCount() const -> size_t { return _rate_matrix.size(); }

  inline auto CLVCount() const -> size_t { return _next_free_clv; }

  inline auto matrixSize() const -> size_t {
    return leadingDimension() * restrictedStateCount();
  }

  inline auto nodeCount() const -> size_t { return _inner_count + _taxa_count; }

  inline auto suggestProbMatrixIndex() -> size_t {
    size_t suggested_index = _prob_matrix_reserved_count;
    _prob_matrix_reserved_count += 1;
    return suggested_index;
  }

  inline auto suggestRateMatrixIndex() -> size_t {
    size_t suggested_index = _rate_matrix_reserved_count;
    _rate_matrix_reserved_count += 1;
    return suggested_index;
  }

  inline auto reserveRateMatrixIndex(size_t index) -> size_t {
    _rate_matrix_reserved_count =
        std::max(index + 1, _rate_matrix_reserved_count);
    return index;
  }

  static inline auto suggestFreqVectorIndex() -> size_t { return 0; }

  inline auto registerGenericCLV() -> size_t { return registerCLV(); }

  inline void registerTopCLV(size_t node_id) {
    _node_reservations[node_id]._top_clv = registerCLV();
  }

  inline void registerChildrenCLV(size_t node_id) {
    _node_reservations[node_id]._bot1_clv = registerCLV();
    _node_reservations[node_id]._bot2_clv = registerCLV();
  }

  void setTipCLV(size_t index, Range dist);

  inline void setNodeRegionMask(size_t node_id, RangeMask rm) {
    _node_reservations[node_id]._child_region_mask = rm;
  }

  inline void registerTopCLVReverse(size_t node_id) {
    _node_reservations[node_id]._top_rclv = registerCLV();
  }

  inline void registerBot1CLVReverse(size_t node_id) {
    _node_reservations[node_id]._bot1_rclv = registerCLV();
  }

  inline void registerBot2CLVReverse(size_t node_id) {
    _node_reservations[node_id]._bot2_rclv = registerCLV();
  }

  inline auto getTopCLV(size_t node_id) -> size_t {
    return _node_reservations[node_id]._top_clv;
  }

  inline auto getTopCLVReverse(size_t node_id) -> size_t {
    return _node_reservations[node_id]._top_rclv;
  }

  inline auto getLeftChildCLV(size_t node_id) -> size_t {
    return _node_reservations[node_id]._bot1_clv;
  }

  inline auto getRightChildCLV(size_t node_id) -> size_t {
    return _node_reservations[node_id]._bot2_clv;
  }

  inline auto getBot1CLVReverse(size_t node_id) -> size_t {
    return _node_reservations[node_id]._bot1_rclv;
  }

  inline auto getBot2CLVReverse(size_t node_id) -> size_t {
    return _node_reservations[node_id]._bot2_rclv;
  }

  inline auto getBaseFrequencies(size_t index) -> LagrangeColVector & {
    return _base_frequencies[index];
  }

  inline void setBaseFrequenciesByDist(size_t index, Range dist) {
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
  void clean();
  void stage();

  inline auto advanceClock() -> ClockTick { return _current_clock++; }

  inline auto readClock() -> ClockTick { return _current_clock; }

  auto getPeriodParams() const -> const std::vector<PeriodParams> & {
    return _periods;
  }

  auto getPeriodParams(size_t period_index) const -> const PeriodParams & {
    return _periods[period_index];
  }

  void setPeriodParams(size_t period_index, double d, double e);

  void setPeriodParamsCount(size_t periods) { _periods.resize(periods); }

  inline auto reserved() const -> bool {
    return !(_base_frequencies == nullptr || _clvs.empty());
  }

  auto reportNodeVectors(size_t node_id) const -> std::string;

  inline auto computeMatrixIndex(size_t i, size_t j) const -> size_t {
    return i * leadingDimension() + j;
  }

  inline void setReversePrior(size_t index) {
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

  inline auto leadingDimension() const -> size_t { return _leading_dim; }

  inline auto restrictedStateCount() const -> size_t {
    return _restricted_state_count;
  }

  inline auto maxAreas() const -> size_t { return _max_areas; }

  inline auto matrixRows() const -> size_t { return restrictedStateCount(); }

  inline auto matrixCols() const -> size_t { return matrixRows(); }

  inline auto CLVSize() const -> size_t { return restrictedStateCount(); }

  auto EXPMExecutionCount() const -> size_t { return _expm_count; }

  auto getAllRangeMask() const -> RangeMask { return (1ul << regions()) - 1; }

 private:
  inline auto registerCLV() -> size_t { return _next_free_clv++; }

  size_t _taxa_count;
  size_t _inner_count;
  size_t _regions;
  size_t _states;
  size_t _max_areas;
  size_t _restricted_state_count;
  size_t _next_free_clv;

  size_t _leading_dim;

  size_t _rate_matrix_reserved_count;
  size_t _prob_matrix_reserved_count;

  std::vector<MatrixReservation> _rate_matrix;
  std::vector<MatrixReservation> _prob_matrix;

  size_t _base_frequencies_count;
  LagrangeColVector *_base_frequencies;

  std::vector<CLVReservation> _clvs;
  size_t *_clv_scalars;

  std::vector<NodeReservation> _node_reservations;

  std::vector<PeriodParams> _periods;

  Clock _current_clock;
  bool _staged;
  bool _reserved;

  size_t _expm_count = 0;
};

}  // namespace lagrange
#endif
