/* Workspace.cpp
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#include <cstddef>

#include "Common.h"
#include "Utils.h"
#include "Workspace.h"

Workspace::~Workspace() {
  for (auto &res : _rate_matrix) {
    delete res._matrix;
  }
  for (auto &res : _prob_matrix) {
    delete res._matrix;
  }
  if (_base_frequencies != nullptr) {
    for (size_t i = 0; i < _base_frequencies_count; ++i) {
      if (_base_frequencies[i] != nullptr) {
        bli_obj_free(_base_frequencies[i]);
      }
    }
    delete[] _base_frequencies;
  }
  if (_clvs != nullptr) {
    for (size_t i = 0; i < clv_count(); ++i) {
      if (_clvs[i] != nullptr) {
        bli_obj_free(_clvs[i]);
      }
    }
    delete[] _clvs;
  }
  if (_clv_scalars != nullptr) delete[] _clv_scalars;
  delete[] _lapack_workspace_buffer;
  delete[] _lapack_piv_buffer;
}

void Workspace::register_top_clv(size_t node_id) {
  _node_reservations[node_id]._top_clv = register_clv();
}

void Workspace::register_children_clv(size_t node_id) {
  _node_reservations[node_id]._bot1_clv = register_clv();
  _node_reservations[node_id]._bot2_clv = register_clv();
}

void Workspace::set_tip_clv(size_t index, lagrange_dist_t dist) {
  bli_setiv_real(_clvs[index * _clv_stride], 1.0, dist);
}

void Workspace::reserve() {
  if (_base_frequencies != nullptr) {
    delete[] _base_frequencies;
  }

  _base_frequencies = new lagrange_col_vector_t[_base_frequencies_count];

  for (size_t i = 0; i < _base_frequencies_count; ++i) {
    bli_obj_create(BLIS_DOUBLE, states(), 1, 0, 0, _base_frequencies[i]);
    bli_setv_real(_base_frequencies[i], 1.0);
  }

  if (_clvs != nullptr) {
    delete[] _clvs;
  }

  _clvs = new lagrange_col_vector_t[clv_count()];

  for (size_t i = 0; i < clv_count(); ++i) {
    bli_obj_create(BLIS_DOUBLE, states(), 1, 0, 0, _clvs[i]);
    bli_setv_real(_clvs[i], 0.0);
  }

  if (_clv_scalars != nullptr) {
    delete[] _clv_scalars;
  }

  _clv_scalars = new size_t[clv_count()];

  for (size_t i = 0; i < clv_count(); i++) {
    _clv_scalars[i] = 0;
  }

  for (size_t i = 0; i < _rate_matrix.size(); i++) {
    if (_rate_matrix[i]._matrix != nullptr) {
      delete _rate_matrix[i]._matrix;
    }

    bli_obj_create(BLIS_DOUBLE, states(), states(), 0, 0,
                   _rate_matrix[i]._matrix);
  }
  for (size_t i = 0; i < _rate_matrix.size(); i++) {
    if (_prob_matrix[i]._matrix != nullptr) {
      delete _prob_matrix[i]._matrix;
    }
    bli_obj_create(BLIS_DOUBLE, states(), states(), 0, 0,
                   _prob_matrix[i]._matrix);
  }

  _lapack_workspace_buffer = new double[states()];
  _lapack_piv_buffer = new int[states()];

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
    oss << blaze::trans(_clvs[reservations._top_clv]);
  }
  if (reservations._bot1_clv != sentinel_value) {
    oss << "bot1 clv (index: " << reservations._bot1_clv << ")\n";
    oss << blaze::trans(_clvs[reservations._bot1_clv]);
  }
  if (reservations._bot2_clv != sentinel_value) {
    oss << "bot2 clv (index: " << reservations._bot2_clv << ")\n";
    oss << blaze::trans(_clvs[reservations._bot2_clv]);
  }
  if (reservations._top_rclv != sentinel_value) {
    oss << "top rclv (index: " << reservations._top_rclv << ")\n";
    oss << blaze::trans(_clvs[reservations._top_rclv]);
  }
  if (reservations._bot1_rclv != sentinel_value) {
    oss << "bot1 rclv (index: " << reservations._bot1_rclv << ")\n";
    oss << blaze::trans(_clvs[reservations._bot1_rclv]);
  }
  if (reservations._bot2_rclv != sentinel_value) {
    oss << "bot2 rclv (index: " << reservations._bot2_rclv << ")\n";
    oss << blaze::trans(_clvs[reservations._bot2_rclv]);
  }
  return oss.str();
}
#endif
