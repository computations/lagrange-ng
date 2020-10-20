#ifndef LAGRANGE_WORKSPACE_H
#define LAGRANGE_WORKSPACE_H

#include <blaze/Math.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>

#include <cstddef>
#include <new>
#include <stdexcept>

#include "Common.h"
#include "RateModel.h"

class Workspace {
public:
  Workspace(size_t taxa_count, size_t inner_count, size_t regions,
            size_t rate_matrix_count, size_t prob_matrix_count)
      : _taxa_count{taxa_count}, _inner_count{inner_count}, _regions{regions},
        _states{1ull << regions}, _rate_matrix_count{rate_matrix_count},
        _prob_matrix_count{prob_matrix_count}, _clv_stride{1} {
    if (taxa_count == 0) {
      throw std::bad_alloc();
    }
    _rate_matrix = new lagrange_matrix_t[_rate_matrix_count];
    _prob_matrix = new lagrange_matrix_t[_prob_matrix_count];
    _clvs = new lagrange_col_vector_t[clv_count()];
    for (size_t i = 0; i < _rate_matrix_count; i++) {
      _rate_matrix[i] = lagrange_matrix_t(_states, _states);
    }
    for (size_t i = 0; i < _prob_matrix_count; i++) {
      _prob_matrix[i] = lagrange_matrix_t(_states, _states);
    }
    for (size_t i = 0; i < clv_count(); i++) {
      _clvs[i * _clv_stride] = lagrange_col_vector_t(_states);
    }
  }

  Workspace(size_t taxa_count, size_t regions)
      : Workspace(taxa_count, taxa_count - 1, regions, 1, 1) {}

  ~Workspace() {
    delete[] _rate_matrix;
    delete[] _prob_matrix;
    delete[] _clvs;
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
  inline size_t clv_count() const { return 3 * (_inner_count) + _taxa_count; }
  inline size_t matrix_size() const { return _states * _states; }

private:
  size_t _taxa_count;
  size_t _inner_count;
  size_t _regions;
  size_t _states;

  size_t _rate_matrix_count;
  lagrange_matrix_t *_rate_matrix;

  size_t _prob_matrix_count;
  lagrange_matrix_t *_prob_matrix;

  size_t _clv_stride;
  lagrange_col_vector_t *_clvs;
};

#endif
