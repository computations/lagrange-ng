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
#include <sstream>
#include <stdexcept>
#include <string>

#include "Common.h"

struct matrix_reservation_t {
  lagrange_matrix_t *_matrix = nullptr;
  lagrange_clock_t _last_update;
  lagrange_op_id_t _op_id;
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
        _rate_matrix{rate_matrix_count},
        _prob_matrix{prob_matrix_count},
        _base_frequencies_count{1},
        _base_frequencies{nullptr},
        _clv_stride{1},
        _clvs{nullptr},
        _node_reservations{node_count()},
        _periods{1},
        _current_clock{0},
        _reserved{false} {
    if (taxa_count == 0) {
      throw std::runtime_error{"We cannot make a workspace with zero taxa"};
    }
  }

  Workspace(size_t taxa_count, size_t regions)
      : Workspace(taxa_count, taxa_count - 1, regions, 1, 1) {}

  ~Workspace();

  inline lagrange_col_vector_t &clv(size_t i) {
    if (i >= clv_count()) {
      throw std::runtime_error{"CLV access out of range"};
    }
    return _clvs[(i * _clv_stride)];
  }

  inline lagrange_matrix_t &rate_matrix(size_t i) {
    if (i >= _rate_matrix.size()) {
      throw std::runtime_error{"Rate matrix access out of range"};
    }
    return *_rate_matrix[i]._matrix;
  }

  inline lagrange_matrix_t &prob_matrix(size_t i) {
    if (i >= _prob_matrix.size()) {
      throw std::runtime_error{"Prob matrix access out of range"};
    }
    return *_prob_matrix[i]._matrix;
  }

  inline size_t states() const { return _states; }
  inline size_t regions() const { return _regions; }
  inline size_t prob_matrix_count() const { return _prob_matrix.size(); }
  inline size_t rate_matrix_count() const { return _rate_matrix.size(); }
  inline size_t clv_count() const { return _next_free_clv; }
  inline size_t matrix_size() const { return _states * _states; }
  inline size_t node_count() const { return _inner_count + _taxa_count; }

  inline size_t suggest_prob_matrix_index() const { return 0; }
  inline size_t suggest_rate_matrix_index() const { return 0; }
  inline size_t suggest_freq_vector_index() const { return 0; }

  inline size_t register_generic_clv() { return register_clv(); }

  void register_top_clv(size_t node_id);

  void register_children_clv(size_t node_id);

  void set_tip_clv(size_t index, lagrange_dist_t dist);

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
    return _node_reservations[node_id]._bot1_rclv;
  }

  inline size_t get_bot2_clv_reverse(size_t node_id) {
    return _node_reservations[node_id]._bot2_rclv;
  }

  inline lagrange_col_vector_t &get_base_frequencies(size_t index) {
    return _base_frequencies[index];
  }

  void reserve();

  inline lagrange_clock_t advance_clock() { return _current_clock++; }

  const period_t &get_period_params(size_t period_index) const {
    return _periods[period_index];
  }

  void set_period_params(size_t period_index, double d, double e);

  inline bool reserved() const {
    return !(_base_frequencies == nullptr || _clvs == nullptr);
  }

  std::string report_node_vecs(size_t node_id) const;

 private:
  inline size_t register_clv() { return _next_free_clv++; }

  size_t _taxa_count;
  size_t _inner_count;
  size_t _regions;
  size_t _states;
  size_t _next_free_clv;

  std::vector<matrix_reservation_t> _rate_matrix;
  std::vector<matrix_reservation_t> _prob_matrix;

  size_t _base_frequencies_count;
  lagrange_col_vector_t *_base_frequencies;

  size_t _clv_stride;
  lagrange_col_vector_t *_clvs;

  std::vector<node_reservation_t> _node_reservations;

  std::vector<period_t> _periods;

  lagrange_clock_t _current_clock;
  bool _reserved;
};

#endif
