/* Workspace.h
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#ifndef LAGRANGE_WORKSPACE_H
#define LAGRANGE_WORKSPACE_H

#include <cstddef>
#include <limits>
#include <new>
#include <stdexcept>

#include "Common.h"
#include "RateModel.h"

struct node_reservation_t {
  node_reservation_t()
      : _top_clv{std::numeric_limits<size_t>::max()},
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

struct period_t {
  double dispersion_rate;
  double extinction_rate;
};

class Workspace {
 public:
  Workspace(size_t taxa_count, size_t inner_count, size_t regions,
            size_t rate_matrix_count, size_t prob_matrix_count)
      : _taxa_count{taxa_count},
        _inner_count{inner_count},
        _regions{regions},
        _states{1ull << regions},
        _next_free_clv{0},
        _rate_matrix_count{rate_matrix_count},
        _rate_matrix{nullptr},
        _prob_matrix_count{prob_matrix_count},
        _prob_matrix{nullptr},
        _base_frequencies_count{1},
        _base_frequencies{nullptr},
        _clv_stride{1},
        _clvs{nullptr},
        _node_reservations{node_count()},
        _periods{1},
        _current_clock{0} {
    if (taxa_count == 0) {
      throw std::runtime_error{"We cannot make a workspace with zero taxa"};
    }
  }

  Workspace(size_t taxa_count, size_t regions)
      : Workspace(taxa_count, taxa_count - 1, regions, 1, 1) {}

  ~Workspace() {
    if (_rate_matrix != nullptr) delete[] _rate_matrix;
    if (_prob_matrix != nullptr) delete[] _prob_matrix;
    if (_base_frequencies != nullptr) delete[] _base_frequencies;
    if (_clvs != nullptr) delete[] _clvs;
  }

  inline lagrange_col_vector_t &clv(size_t i) {
    if (i >= clv_count()) {
      throw std::runtime_error{"CLV access out of range"};
    }
    return _clvs[(i * _clv_stride)];
  }

  inline lagrange_matrix_t &rate_matrix(size_t i) {
    if (i >= _rate_matrix_count) {
      throw std::runtime_error{"Rate matrix access out of range"};
    }
    return _rate_matrix[i];
  }

  inline lagrange_matrix_t &prob_matrix(size_t i) {
    if (i >= _prob_matrix_count) {
      throw std::runtime_error{"Prob matrix access out of range"};
    }
    return _prob_matrix[i];
  }

  inline size_t states() const { return _states; }
  inline size_t regions() const { return _regions; }
  inline size_t prob_matrix_count() const { return _prob_matrix_count; }
  inline size_t rate_matrix_count() const { return _rate_matrix_count; }
  inline size_t clv_count() const { return _next_free_clv; }
  inline size_t matrix_size() const { return _states * _states; }
  inline size_t node_count() const { return _inner_count + _taxa_count; }

  inline size_t suggest_prob_matrix_index() const { return 0; }
  inline size_t suggest_rate_matrix_index() const { return 0; }

  inline size_t register_generic_clv() { return register_clv(); }

  void register_top_clv(size_t node_id) {
    _node_reservations[node_id]._top_clv = register_clv();
  }

  void register_children_clv(size_t node_id) {
    _node_reservations[node_id]._bot1_clv = register_clv();
    _node_reservations[node_id]._bot2_clv = register_clv();
  }

  void set_tip_clv(size_t index, lagrange_dist_t dist) {
    _clvs[index * _clv_stride][dist] = 1.0;
  }

  inline void register_top_clv_reverse(size_t node_id) {
    _node_reservations[node_id]._top_rclv = register_clv();
  }

  inline void register_bot1_clv_reverse(size_t node_id) {
    _node_reservations[node_id]._bot1_rclv = register_clv();
  }

  inline void register_bot2_clv_reverse(size_t node_id) {
    _node_reservations[node_id]._bot2_rclv = register_clv();
  }

  inline size_t get_top_clv(size_t node_id) {
    return _node_reservations[node_id]._top_clv;
  }

  inline size_t get_top_clv_reverse(size_t node_id) {
    return _node_reservations[node_id]._top_rclv;
  }

  inline size_t get_lchild_clv(size_t node_id) {
    return _node_reservations[node_id]._bot1_clv;
  }

  inline size_t get_rchild_clv(size_t node_id) {
    return _node_reservations[node_id]._bot2_clv;
  }

  inline size_t get_bot1_clv_reverse(size_t node_id) {
    return _node_reservations[node_id]._bot1_clv;
  }

  inline size_t get_bot2_clv_reverse(size_t node_id) {
    return _node_reservations[node_id]._bot2_clv;
  }

  inline lagrange_col_vector_t &get_base_frequencies(size_t index) {
    return _base_frequencies[index];
  }

  void reserve() {
    _rate_matrix = new lagrange_matrix_t[_rate_matrix_count];
    _prob_matrix = new lagrange_matrix_t[_prob_matrix_count];
    _base_frequencies = new lagrange_col_vector_t[_base_frequencies_count];
    _clvs = new lagrange_col_vector_t[clv_count()];
    for (size_t i = 0; i < _rate_matrix_count; i++) {
      _rate_matrix[i] = lagrange_matrix_t(_states, _states);
    }
    for (size_t i = 0; i < _prob_matrix_count; i++) {
      _prob_matrix[i] = lagrange_matrix_t(_states, _states);
    }
    for (size_t i = 0; i < _base_frequencies_count; i++) {
      _base_frequencies[i] = lagrange_col_vector_t(_states, 1.0 / _states);
    }
    for (size_t i = 0; i < clv_count(); i++) {
      _clvs[i * _clv_stride] = lagrange_col_vector_t(_states);
      _clvs[i * _clv_stride] = 0.0;
    }
  }

  lagrange_clock_t advance_clock() { return _current_clock++; }

  period_t get_period_params(size_t period_index) const {
    return _periods[period_index];
  }

  void set_period_params(size_t period_index, double d, double e) {
    _periods[period_index] = {d, e};
  }

 private:
  inline size_t register_clv() { return _next_free_clv++; }

  size_t _taxa_count;
  size_t _inner_count;
  size_t _regions;
  size_t _states;
  size_t _next_free_clv;

  size_t _rate_matrix_count;
  lagrange_matrix_t *_rate_matrix;

  size_t _prob_matrix_count;
  lagrange_matrix_t *_prob_matrix;

  size_t _base_frequencies_count;
  lagrange_col_vector_t *_base_frequencies;

  size_t _clv_stride;
  lagrange_col_vector_t *_clvs;

  std::vector<node_reservation_t> _node_reservations;

  std::vector<period_t> _periods;

  lagrange_clock_t _current_clock;
};

#endif
