/*
 * RateMatrix.h
 *
 *  Created on: Aug 14, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef RATEMODEL_H_
#define RATEMODEL_H_

#include <map>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

#include "AncSplit.h"
#include "Common.h"
#include "Utils.h"
#include <blaze/Math.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>

typedef blaze::DynamicMatrix<double, blaze::columnMajor> lagrange_matrix_t;
typedef blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>
    lagrange_complex_matrix_t;
typedef blaze::DynamicVector<double, blaze::columnVector> lagrange_col_vector_t;
typedef blaze::DynamicVector<std::complex<double>, blaze::columnVector>
    lagrange_complex_col_vector_t;

namespace std {
template <> struct hash<std::vector<int>> {
public:
  size_t operator()(const std::vector<int> &vec) const;
};
} // namespace std

class RateModel {
private:
  bool _global_ext;
  unsigned int _area_count;
  lagrange_dist_t _valid_dist_mask;
  int _thread_count;
  vector<string> _labels;
  vector<double> _periods;
  vector<lagrange_dist_t> _dists;
  unordered_map<lagrange_dist_t, int> _dists_int_map;
  unordered_map<int, lagrange_dist_t> _int_dists_map;
  unordered_map<lagrange_dist_t, vector<lagrange_region_split_t>> _iter_dists;
  vector<vector<vector<double>>> _dispersal_params;
  vector<vector<vector<double>>> _dispersal_params_mask;
  vector<vector<double>> _extinction_params;
  inline size_t getDistCount() const { return (1ul << _area_count); }

  vector<lagrange_matrix_t> _rate_matrix;
  vector<lagrange_matrix_t> _rate_matrix_transposed;

  vector<int> _active_zone_counts;
  vector<vector<int>> _ia_s;
  vector<vector<int>> _ja_s;
  vector<vector<double>> _a_s;
  void iter_all_dist_splits();

  lagrange_matrix_t compute_matrix_exponential_ss(lagrange_matrix_t A) const;
  lagrange_matrix_t
  compute_matrix_exponential_eigen(const lagrange_matrix_t &A) const;

  size_t _expm_count;

public:
  RateModel(int na, bool ge, vector<double> pers, bool);
  void set_nthreads(int nthreads);
  int get_nthreads();
  size_t get_expm_count();
  void setup_dists();
  void setup_dists(vector<lagrange_dist_t>, bool);
  void setup_Dmask();
  void setup_D_provided(double d, vector<vector<vector<double>>> &D_mask_in);
  void set_Dmask_cell(int period, int area, int area2, double prob, bool sym);
  void setup_D(double d);
  void setup_E(double e);
  void set_Qdiag(int period);
  void setup_Q();
  vector<vector<double>> setup_fortran_P(int period, double t,
                                         bool store_p_matrices);
  vector<vector<double>> setup_sparse_full_P(int period, double t);
  vector<double> setup_sparse_single_column_P(int period, double t, int column);
  vector<vector<double>> setup_pthread_sparse_P(int period, double t,
                                                vector<int> &columns);
  string Q_repr(int period);
  string P_repr(int period);
  vector<lagrange_dist_t> enumerate_dists();
  vector<lagrange_region_split_t> iter_dist_splits(lagrange_dist_t dist);
  const vector<lagrange_dist_t> &getDists();
  size_t getDistsSize() const;
  const unordered_map<lagrange_dist_t, int> &get_dists_int_map();
  const unordered_map<int, lagrange_dist_t> &get_int_dists_map();
  const vector<lagrange_region_split_t> &
  get_iter_dist_splits(lagrange_dist_t dist) const;
  void remove_dist(vector<int> dist);

  bool _sparse;

  size_t get_num_areas();
  int get_num_periods();

  vector<AncSplit> iter_ancsplits(lagrange_dist_t dist);
  void iter_ancsplits_just_int(lagrange_dist_t &dist, vector<int> &leftdists,
                               vector<int> &rightdists, double &weight);

  vector<int> get_columns_for_sparse(vector<double> &);
  vector<int> get_columns_for_sparse(vector<Superdouble> &);

  /*
   testing storing once optimization has occured map of period and map of bl and
   p matrix map<period,map<branch length,p matrix>>
   */
  unordered_map<int, map<double, vector<vector<double>>>> stored_p_matrices;

  /*
   * get things from stmap
   */
  vector<lagrange_matrix_t> &get_Q();
  // this should be used for getting the eigenvectors and eigenvalues
  bool get_eigenvec_eigenval_from_Q(lagrange_complex_matrix_t &eigenvalues,
                                    lagrange_complex_matrix_t &eigenvectors,
                                    int period);
};

#endif /* RATEMATRIX_H_ */
