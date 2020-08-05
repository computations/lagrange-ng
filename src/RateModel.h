/*
 * RateMatrix.h
 *
 *  Created on: Aug 14, 2009
 *      Author: smitty
 */

#ifndef RATEMODEL_H_
#define RATEMODEL_H_

#include <map>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

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
  int _area_count;
  int _thread_count;
  vector<string> _labels;
  vector<double> _periods;
  vector<vector<int>> _dists;
  unordered_map<vector<int>, vector<vector<vector<int>>>> _iter_dists;
  unordered_map<vector<int>, int> _dists_int_map;
  unordered_map<int, vector<int>> _int_dists_map;
  vector<vector<vector<double>>> _dispersal_params;
  vector<vector<vector<double>>> _dispersal_params_mask;
  vector<vector<double>> _extinction_params;

  vector<lagrange_matrix_t> _rate_matrix;
  vector<lagrange_matrix_t> _rate_matrix_transposed;

  vector<int> _active_zone_counts;
  vector<vector<int>> _ia_s;
  vector<vector<int>> _ja_s;
  vector<vector<double>> _a_s;
  void iter_all_dist_splits();

  lagrange_matrix_t compute_matrix_exponential_ss(lagrange_matrix_t A) const;

public:
  RateModel(int na, bool ge, vector<double> pers, bool);
  void set_nthreads(int nthreads);
  int get_nthreads();
  void setup_dists();
  void setup_dists(vector<vector<int>>, bool);
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
  vector<vector<int>> enumerate_dists();
  vector<vector<vector<int>>> iter_dist_splits(vector<int> &dist);
  const vector<vector<int>> &getDists();
  size_t getDistsSize() const;
  const unordered_map<vector<int>, int> &get_dists_int_map();
  const unordered_map<int, vector<int>> &get_int_dists_map();
  const vector<vector<vector<int>>> &
  get_iter_dist_splits(const vector<int> &dist) const;
  void remove_dist(vector<int> dist);
  bool _sparse;
  int get_num_areas();
  int get_num_periods();
  /*
   testing storing once optimization has occured
   map of period and map of bl and p matrix
   map<period,map<branch length,p matrix>>
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
