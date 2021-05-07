/* Workspace.cpp
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#include <limits>

#include "Common.h"
#include "Utils.h"
#include "Workspace.h"

Workspace::~Workspace() {
  for (auto &res : _rate_matrix) { delete[] res._matrix; }

  for (auto &res : _prob_matrix) { delete[] res._matrix; }

  if (_base_frequencies != nullptr) {
    for (size_t i = 0; i < _base_frequencies_count; ++i) {
      if (_base_frequencies[i] != nullptr) { delete[] _base_frequencies[i]; }
    }
    delete[] _base_frequencies;
  }

  for (size_t i = 0; i < _clvs.size(); ++i) {
    if (_clvs[i]._clv == nullptr) { continue; }
    delete[] _clvs[i]._clv;
  }
  if (_clv_scalars != nullptr) delete[] _clv_scalars;
}

void Workspace::register_top_clv(size_t node_id) {
  _node_reservations[node_id]._top_clv = register_clv();
}

void Workspace::register_children_clv(size_t node_id) {
  _node_reservations[node_id]._bot1_clv = register_clv();
  _node_reservations[node_id]._bot2_clv = register_clv();
}

void Workspace::set_tip_clv(size_t index, lagrange_dist_t dist) {
  _clvs[index]._clv[dist] = 1.0;
  _clvs[index]._last_update = std::numeric_limits<size_t>::max();
}

void Workspace::reserve() {
  if (_base_frequencies != nullptr) {
    for (size_t i = 0; i < _base_frequencies_count; ++i) {
      delete[] _base_frequencies[i];
    }
    delete[] _base_frequencies;
  }

  _base_frequencies = new lagrange_col_vector_t[_base_frequencies_count];

  for (size_t i = 0; i < _base_frequencies_count; ++i) {
    _base_frequencies[i] = new lagrange_matrix_base_t[states()];
    for (size_t j = 0; j < states(); j++) { _base_frequencies[i][j] = 1.0; }
  }

  _clvs.resize(clv_count());

  for (size_t i = 0; i < _clvs.size(); ++i) {
    if (_clvs[i]._clv != nullptr) { delete[] _clvs[i]._clv; }

    _clvs[i]._clv = new lagrange_matrix_base_t[states()];

    for (size_t j = 0; j < states(); j++) { _clvs[i]._clv[j] = 0.0; }
  }

  if (_clv_scalars != nullptr) { delete[] _clv_scalars; }

  _clv_scalars = new size_t[clv_count()];

  for (size_t i = 0; i < clv_count(); i++) { _clv_scalars[i] = 0; }

  for (size_t i = 0; i < _rate_matrix.size(); i++) {
    if (_rate_matrix[i]._matrix != nullptr) {
      delete[] _rate_matrix[i]._matrix;
    }
    _rate_matrix[i]._matrix = new lagrange_matrix_base_t[matrix_size()];
  }

  for (size_t i = 0; i < _prob_matrix.size(); i++) {
    if (_prob_matrix[i]._matrix != nullptr) {
      delete[] _prob_matrix[i]._matrix;
    }
    _prob_matrix[i]._matrix = new lagrange_matrix_base_t[matrix_size()];
  }

  _reserved = true;
}

void Workspace::set_period_params(size_t period_index, double d, double e) {
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
