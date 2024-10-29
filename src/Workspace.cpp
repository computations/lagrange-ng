/* Workspace.cpp
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#include "Workspace.hpp"

#include <limits>

#include "Common.hpp"

namespace lagrange {

Workspace::~Workspace() {
  for (auto &res : _rate_matrix) { delete[] res._matrix; }

  for (auto &res : _prob_matrix) { delete[] res._matrix; }

  if (_base_frequencies != nullptr) {
    for (size_t i = 0; i < _base_frequencies_count; ++i) {
      if (_base_frequencies[i] != nullptr) { delete[] _base_frequencies[i]; }
    }
    delete[] _base_frequencies;
  }

  for (auto &_clv : _clvs) {
    if (_clv._clv == nullptr) { continue; }
    delete[] _clv._clv;
  }
  delete[] _clv_scalars;
}

void Workspace::registerTopCLV(size_t node_id) {
  _node_reservations[node_id]._top_clv = registerCLV();
}

void Workspace::registerChildrenCLV(size_t node_id) {
  _node_reservations[node_id]._bot1_clv = registerCLV();
  _node_reservations[node_id]._bot2_clv = registerCLV();
}

void Workspace::setTipCLV(size_t index, Range clv_index) {
  if (clv_index >= restrictedStateCount()) {
    throw std::runtime_error{
        "Attempted to set a state that is too large for this dataset"};
  }
  _clvs[index]._clv[clv_index] = 1.0;
  _clvs[index]._last_update = std::numeric_limits<size_t>::max();
}

/**
 * Prep for allocation by setting up the various reservations
*/
void Workspace::stage() { 
  _clvs.resize(CLVCount()); 
  _rate_matrix.resize(_rate_matrix_reserved_count);
  _prob_matrix.resize(_prob_matrix_reserved_count);
  _staged = true;
}

void Workspace::clean() {
  if (_base_frequencies != nullptr) {
    for (size_t i = 0; i < _base_frequencies_count; ++i) {
      delete[] _base_frequencies[i];
    }
    delete[] _base_frequencies;
  }

  for (auto &c : _clvs) { delete[] c._clv; }

  delete[] _clv_scalars;

  for (auto &rm : _rate_matrix) { delete[] rm._matrix; }

  for (auto &pm : _prob_matrix) { delete[] pm._matrix; }
  _reserved = false;
}

void Workspace::reserve() {
  assert(!_reserved);
  stage();
  _base_frequencies = new LagrangeColVector[_base_frequencies_count];

  for (size_t i = 0; i < _base_frequencies_count; ++i) {
    _base_frequencies[i] = new LagrangeMatrixBase[CLVSize()];
    for (size_t j = 0; j < restrictedStateCount(); j++) {
      _base_frequencies[i][j] = 1.0;
    }
  }

  for (auto &c : _clvs) {
    c._clv = new LagrangeMatrixBase[CLVSize()];

    for (size_t j = 0; j < restrictedStateCount(); j++) { c._clv[j] = 0.0; }
  }

  _clv_scalars = new size_t[CLVCount()];

  for (size_t i = 0; i < CLVCount(); i++) { _clv_scalars[i] = 0; }

  for (auto &rm : _rate_matrix) {
    rm._matrix = new LagrangeMatrixBase[matrixSize()];
  }

  for (auto &pm : _prob_matrix) {
    pm._matrix = new LagrangeMatrixBase[matrixSize()];
  }

  _reserved = true;
}

void Workspace::setPeriodParams(size_t period_index, double d, double e) {
  _periods[period_index] = {d, e};
}

#if 0
std::string Workspace::report_node_vecs(size_t node_id) const {
  auto reservations = _node_reservations[node_id];
  std::stringstream oss;
  constexpr size_t sentinel_value = std::numeric_limits<size_t>::max();
  if (reservations._top_clv != sentinel_value) {
    oss << "top clv (index: " << reservations._top_clv << ")\n";
    oss << blaze::trans(_clvs[reservations._top_clv]._clv);
  }
  if (reservations._bot1_clv != sentinel_value) {
    oss << "bot1 clv (index: " << reservations._bot1_clv << ")\n";
    oss << blaze::trans(_clvs[reservations._bot1_clv]._clv);
  }
  if (reservations._bot2_clv != sentinel_value) {
    oss << "bot2 clv (index: " << reservations._bot2_clv << ")\n";
    oss << blaze::trans(_clvs[reservations._bot2_clv]._clv);
  }
  if (reservations._top_rclv != sentinel_value) {
    oss << "top rclv (index: " << reservations._top_rclv << ")\n";
    oss << blaze::trans(_clvs[reservations._top_rclv]._clv);
  }
  if (reservations._bot1_rclv != sentinel_value) {
    oss << "bot1 rclv (index: " << reservations._bot1_rclv << ")\n";
    oss << blaze::trans(_clvs[reservations._bot1_rclv]._clv);
  }
  if (reservations._bot2_rclv != sentinel_value) {
    oss << "bot2 rclv (index: " << reservations._bot2_rclv << ")\n";
    oss << blaze::trans(_clvs[reservations._bot2_rclv]._clv);
  }
  return oss.str();
}
#endif
}  // namespace lagrange
