/* Workspace.h
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#ifndef LAGRANGE_WORKSPACE_H
#define LAGRANGE_WORKSPACE_H

#include <cstddef>
#include <limits>
#include <mutex>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include "Common.h"
#include "Utils.h"

struct matrix_reservation_t {
  lagrange_matrix_t _matrix = nullptr;
  lagrange_clock_tick_t _last_update = 0;
  lagrange_op_id_t _op_id;
  std::unique_ptr<std::mutex> _matrix_lock;
};

struct clv_reservation_t {
  lagrange_col_vector_t _clv = nullptr;
  lagrange_clock_tick_t _last_update = 0;
};

class Workspace {
 public:
  Workspace(size_t taxa_count, size_t inner_count, size_t regions,
            size_t max_areas)
      : _taxa_count{taxa_count},
        _inner_count{inner_count},
        _regions{regions},
        _states{1ull << regions},
        _max_areas{max_areas},
        _next_free_clv{0},
        _leading_dim{_states},
        _rate_matrix{1},
        _prob_matrix{1},
        _base_frequencies_count{1},
        _base_frequencies{nullptr},
        _clvs{},
        _clv_scalars{nullptr},
        _node_reservations{node_count()},
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

  Workspace(size_t taxa_count, size_t regions, size_t max_areas)
      : Workspace(taxa_count, taxa_count - 1, regions, max_areas) {}

  ~Workspace();

  inline size_t &clv_scalar(size_t i) { return _clv_scalars[i]; }

  inline lagrange_matrix_t &rate_matrix(size_t i) {
    if (i >= _rate_matrix.size()) {
      throw std::runtime_error{"Rate matrix access out of range"};
    }
    return _rate_matrix[i]._matrix;
  }

  inline void update_rate_matrix(size_t index, const lagrange_matrix_t &A) {
    if (index >= _rate_matrix.size()) {
      throw std::runtime_error{"Rate matrix access out of range when updating"};
    }

    for (size_t i = 0; i < matrix_size(); i++) {
      _rate_matrix[index]._matrix[i] = A[i];
    }

    _rate_matrix[index]._last_update = advance_clock();
  }

  inline void update_rate_matrix_clock(size_t i) {
    _rate_matrix[i]._last_update = advance_clock();
  }

  inline lagrange_matrix_t &prob_matrix(size_t i) {
    if (i >= _prob_matrix.size()) {
      throw std::runtime_error{"Prob matrix access out of range"};
    }
    return _prob_matrix[i]._matrix;
  }

  inline lagrange_clock_tick_t update_prob_matrix(size_t index,
                                                  const lagrange_matrix_t &A) {
    if (index >= _prob_matrix.size()) {
      throw std::runtime_error{"Prob matrix access out of range when updating"};
    }

    for (size_t i = 0; i < matrix_size(); i++) {
      _prob_matrix[index]._matrix[i] = A[i];
    }

    return _prob_matrix[index]._last_update = advance_clock();
  }

  inline lagrange_clock_tick_t last_update_prob_matrix(size_t i) {
    return _prob_matrix[i]._last_update;
  }

  inline lagrange_clock_tick_t last_update_rate_matrix(size_t i) {
    return _rate_matrix[i]._last_update;
  }

  inline const lagrange_col_vector_t &clv(size_t i) {
    if (i >= clv_count()) {
      throw std::runtime_error{"CLV access out of range"};
    }
    return _clvs[i]._clv;
  }

  inline void update_clv(const lagrange_col_vector_t &clv, size_t index) {
    if (index >= clv_count()) {
      throw std::runtime_error{"CLV access out of range"};
    }

    for (size_t i = 0; i < restricted_state_count(); i++) {
      _clvs[index]._clv[i] = clv[i];
    }

    _clvs[index]._last_update = advance_clock();
  }

  inline std::tuple<lagrange_col_vector_t, size_t> clv_size_tuple(
      size_t index) {
    return std::make_tuple(clv(index), clv_size());
  }

  inline void update_clv_clock(size_t index) {
    _clvs[index]._last_update = advance_clock();
  }

  inline lagrange_clock_tick_t last_update_clv(size_t index) {
    return _clvs[index]._last_update;
  }

  inline size_t states() const { return _states; }
  inline size_t regions() const { return _regions; }
  inline size_t prob_matrix_count() const { return _prob_matrix.size(); }
  inline size_t rate_matrix_count() const { return _rate_matrix.size(); }
  inline size_t clv_count() const { return _next_free_clv; }
  inline size_t matrix_size() const {
    return leading_dimension() * restricted_state_count();
  }
  inline size_t node_count() const { return _inner_count + _taxa_count; }

  inline size_t suggest_prob_matrix_index() {
    size_t suggested_index = _prob_matrix.size();
    _prob_matrix.emplace_back();
    return suggested_index;
  }

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

  inline lagrange_clock_tick_t advance_clock() { return _current_clock++; }
  inline lagrange_clock_tick_t read_clock() { return _current_clock; }

  const period_t &get_period_params(size_t period_index) const {
    return _periods[period_index];
  }

  void set_period_params(size_t period_index, double d, double e);

  inline bool reserved() const {
    return !(_base_frequencies == nullptr || _clvs.size() == 0);
  }

  std::string report_node_vecs(size_t node_id) const;

  inline size_t compute_matrix_index(size_t i, size_t j) {
    return i * leading_dimension() + j;
  }

  inline void set_reverse_prior(size_t index) {
    if (index >= _clvs.size()) {
      throw std::runtime_error{
          "Attempting to access a clv that doesn't exist when setting a "
          "reverse prior"};
    }

    for (size_t i = 0; i < restricted_state_count(); i++) {
      _clvs[index]._clv[i] = 1.0;
    }

    _clvs[index]._last_update =
        std::numeric_limits<lagrange_clock_tick_t>::max();
  }

  inline size_t leading_dimension() const { return _leading_dim; }

  inline size_t restricted_state_count() const {
    return _restricted_state_count;
  }

  inline size_t max_areas() const { return _max_areas; }

  inline size_t matrix_rows() const { return restricted_state_count(); }

  inline size_t clv_size() const { return restricted_state_count(); }

  inline std::mutex &get_rate_matrix_lock(size_t index) {
    return *_rate_matrix[index]._matrix_lock;
  }

  inline std::mutex &get_prob_matrix_lock(size_t index) {
    return *_prob_matrix[index]._matrix_lock;
  }

  bool timepoint_is_stale(lagrange_clock_tick_t tick) {
    return tick < _current_clock;
  }

 private:
  inline size_t register_clv() { return _next_free_clv++; }

  size_t _taxa_count;
  size_t _inner_count;
  size_t _regions;
  size_t _states;
  size_t _max_areas;
  size_t _restricted_state_count;
  size_t _next_free_clv;

  size_t _leading_dim;

  std::vector<matrix_reservation_t> _rate_matrix;
  std::vector<matrix_reservation_t> _prob_matrix;

  size_t _base_frequencies_count;
  lagrange_col_vector_t *_base_frequencies;

  std::vector<clv_reservation_t> _clvs;
  size_t *_clv_scalars;

  std::vector<node_reservation_t> _node_reservations;

  std::vector<period_t> _periods;

  lagrange_clock_t _current_clock;
  bool _reserved;
};

#endif
