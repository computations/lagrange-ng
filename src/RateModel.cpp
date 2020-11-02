/*
 * RateMatrix.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include "AncSplit.h"
#include "Common.h"
#define VERBOSE false

#include <blaze/math/IdentityMatrix.h>
#include <pthread.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>

#include "RateMatrixUtils.h"
#include "RateModel.h"
#include "Utils.h"
using namespace std;

/*
  all of these are for sparse matrix calculation and are used for the fortran
  methods

  WILL PROBABLY CHANGE WHEN MOVED TO C++ CLASSES FOR MATEXP
 */
int get_size_for_coo(const lagrange_matrix_t &inmatrix) {
  int count = 0;
  for (size_t i = 0; i < inmatrix.columns(); ++i) {
    for (auto it = inmatrix.begin(i); it != inmatrix.end(i); ++it) {
      if (*it != 0) {
        count += 1;
      }
    }
  }
  return count;
}

void convert_matrix_to_coo_for_fortran(const lagrange_matrix_t &inmatrix,
                                       int *ia, int *ja, double *a) {
  int count = 0;
  for (unsigned int i = 0; i < inmatrix.rows(); i++) {
    for (unsigned int j = 0; j < inmatrix.columns(); j++) {
      if (inmatrix(i, j) != 0.0) {
        ia[count] = i + 1;
        ja[count] = j + 1;
        a[count] = inmatrix(i, j);
        count += 1;
      }
    }
  }
}

void convert_matrix_to_coo_for_fortran_vector(const lagrange_matrix_t &inmatrix,
                                              vector<int> &ia, vector<int> &ja,
                                              vector<double> &a) {
  int count = 0;
  for (size_t i = 0; i < inmatrix.rows(); ++i) {
    for (size_t j = 0; j < inmatrix.columns(); ++j) {
      if (inmatrix(i, j) != 0.0) {
        ia[count] = i + 1;
        ja[count] = j + 1;
        a[count] = inmatrix(i, j);
        count += 1;
      }
    }
  }
}

void convert_matrix_to_single_row_for_fortran(vector<vector<double>> &inmatrix,
                                              double t, double *H) {
  int count = 0;
  for (unsigned int i = 0; i < inmatrix.size(); i++) {
    for (unsigned int j = 0; j < inmatrix[i].size(); j++) {
      H[i + (j * inmatrix.size())] = inmatrix[i][j] * t;
      count += 1;
    }
  }
}

RateModel::RateModel(int na, bool ge, vector<double> pers, bool sp)
    : _global_ext(ge),
      _area_count(na),
      _valid_dist_mask{(1ull << _area_count) - 1},
      _thread_count(0),
      _periods(pers),
      _expm_count{0},
      _sparse(sp) {}

void RateModel::set_nthreads(int nthreads) { _thread_count = nthreads; }

int RateModel::get_nthreads() { return _thread_count; }

void RateModel::setup_dists() {
  map<int, lagrange_dist_t> a = iterate_all_bv(_area_count);
  if (_global_ext) {
    _dists.push_back(0);
  }
  for (lagrange_dist_t i = 1; i < (1ul << _area_count); i++) {
    _dists.push_back(i);
  }
  /*
   calculate the distribution map
   */
  for (unsigned int i = 0; i < _dists.size(); i++) {
    _dists_int_map[_dists[i]] = i;
  }
  for (unsigned int i = 0; i < _dists.size(); i++) {
    _int_dists_map[i] = _dists[i];
  }
  /*
   * precalculate the iterdists
   */
  iter_all_dist_splits();
}

/*
 * need to make a generator function for setting distributions
 */
void RateModel::setup_dists(vector<lagrange_dist_t> indists, bool include) {
  if (include == true) {
    _dists = indists;
    if (_dists[0] > 0) {
      _dists.push_back(0);
    }
  } else {  // exclude is sent
    _dists.push_back(0);

    for (lagrange_dist_t i = 1; i < getDistCount(); ++i) {
      if (std::find(indists.begin(), indists.end(), i) == indists.end()) {
        _dists.push_back(i);
      }
    }
  }
  /*
   calculate the distribution map
   */
  for (unsigned int i = 0; i < _dists.size(); i++) {
    _dists_int_map[_dists[i]] = i;
  }
  for (unsigned int i = 0; i < _dists.size(); i++) {
    _int_dists_map[i] = _dists[i];
  }
  /*
  precalculate the iterdists
   */
  iter_all_dist_splits();
}

/*
 * just give Dmask a bunch of ones
 * specify particular ones in the Dmask_cell
 */
void RateModel::setup_Dmask() {
  vector<double> cols(_area_count, 1);
  vector<vector<double>> rows(_area_count, cols);
  _dispersal_params_mask =
      vector<vector<vector<double>>>(_periods.size(), rows);
}

void RateModel::set_Dmask_cell(int period, int area, int area2, double prob,
                               bool sym) {
  _dispersal_params_mask[period][area][area2] = prob;
  if (sym) _dispersal_params_mask[period][area2][area] = prob;
}

void RateModel::setup_D(double d) {
  vector<double> cols(_area_count, 1 * d);
  vector<vector<double>> rows(_area_count, cols);
  _dispersal_params = vector<vector<vector<double>>>(_periods.size(), rows);
  for (unsigned int i = 0; i < _dispersal_params.size(); i++) {
    for (unsigned int j = 0; j < _dispersal_params[i].size(); j++) {
      _dispersal_params[i][j][j] = 0.0;
      for (unsigned int k = 0; k < _dispersal_params[i][j].size(); k++) {
        _dispersal_params[i][j][k] =
            _dispersal_params[i][j][k] * _dispersal_params_mask[i][j][k];
      }
    }
  }
  if (VERBOSE) {
    cout << "D" << endl;
    for (unsigned int i = 0; i < _dispersal_params.size(); i++) {
      for (unsigned int j = 0; j < _dispersal_params[i].size(); j++) {
        for (unsigned int k = 0; k < _dispersal_params[i][j].size(); k++) {
          cout << _dispersal_params[i][j][k] << " ";
        }
        cout << endl;
      }
      cout << endl;
    }
  }
}

/*
 * this is for estimating the D matrix
 * in this case, a setup dmatrix is being sent and
 * the dmask is then applied to it
 */

void RateModel::setup_D_provided(double d,
                                 vector<vector<vector<double>>> &D_mask_in) {
  vector<double> cols(_area_count, 1 * d);
  vector<vector<double>> rows(_area_count, cols);
  _dispersal_params = vector<vector<vector<double>>>(_periods.size(), rows);
  for (unsigned int i = 0; i < _dispersal_params.size(); i++) {
    for (unsigned int j = 0; j < _dispersal_params[i].size(); j++) {
      _dispersal_params[i][j][j] = 0.0;
      for (unsigned int k = 0; k < _dispersal_params[i][j].size(); k++) {
        _dispersal_params[i][j][k] = _dispersal_params[i][j][k] *
                                     _dispersal_params_mask[i][j][k] *
                                     D_mask_in[i][j][k];
      }
    }
  }
  if (VERBOSE) {
    cout << "D" << endl;
    for (unsigned int i = 0; i < _dispersal_params.size(); i++) {
      for (unsigned int j = 0; j < _dispersal_params[i].size(); j++) {
        for (unsigned int k = 0; k < _dispersal_params[i][j].size(); k++) {
          cout << _dispersal_params[i][j][k] << " ";
        }
        cout << endl;
      }
      cout << endl;
    }
  }
}

void RateModel::setup_E(double e) {
  vector<double> cols(_area_count, 1 * e);
  _extinction_params = vector<vector<double>>(_periods.size(), cols);
}

void RateModel::set_Qdiag(int period) {
  for (unsigned int i = 0; i < _dists.size(); i++) {
    double sum = 0.0;
    for (size_t j = 0; j < _dists.size(); j++) {
      sum += _rate_matrix[period](i, j);
    }
    sum -= _rate_matrix[period](i, i);
    sum *= -1;

    _rate_matrix[period](i, i) = sum;
  }
}

void RateModel::setup_Q() {
  _rate_matrix = vector<lagrange_matrix_t>(_periods.size(),
                                           {_dists.size(), _dists.size(), 0.0});
  for (unsigned int period = 0; period < _rate_matrix.size();
       period++) {  // periods
    for (unsigned int dist_i = 0; dist_i < _dists.size(); dist_i++) {  // dists
      int s1 = lagrange_popcount(_dists[dist_i]);
      if (s1 > 0) {
        for (unsigned int dist_j = 0; dist_j < _dists.size();
             dist_j++) {  // dists
          int sxor = lagrange_popcount(_dists[dist_i] ^ _dists[dist_j]);
          if (sxor == 1) {
            int s2 = lagrange_popcount(_dists[dist_j]);
            int dest = __builtin_ctzll((_dists[dist_i] ^ _dists[dist_j]));
            double rate = 0.0;
            if (s1 < s2) {
              for (unsigned int src = 0;
                   src < static_cast<unsigned int>(_area_count); src++) {
                rate += _dispersal_params[period][src][dest] *
                        ((_dists[dist_i] >> src) & 1ul);
                /*
                if (_dists[dist_i] & (1ul << _area_count) != 0) {
                  rate += _dispersal_params[period][src][dest];
                }
                */
              }
            } else {
              rate = _extinction_params[period][dest];
            }
            _rate_matrix[period](dist_i, dist_j) = rate;
          }
        }
      }
    }
    set_Qdiag(period);
  }

  /*
   * sparse needs to be transposed for matrix exponential calculation
   */

  if (_sparse == true) {
    _rate_matrix_transposed = vector<lagrange_matrix_t>(_periods.size());
    for (unsigned int p = 0; p < _rate_matrix_transposed.size(); p++) {
      _rate_matrix_transposed[p] = blaze::trans(_rate_matrix[p]);
    }

    // setting up the coo numbs
    _active_zone_counts = vector<int>(_rate_matrix.size(), 0);
    for (unsigned int p = 0; p < _rate_matrix.size(); p++) {  // periods
      _active_zone_counts[p] = get_size_for_coo(_rate_matrix[p]);
    }

    // setup matrix
    _ia_s.clear();
    _ja_s.clear();
    _a_s.clear();
    for (unsigned int p = 0; p < _rate_matrix.size(); p++) {  // periods
      vector<int> ia = vector<int>(_active_zone_counts[p]);
      vector<int> ja = vector<int>(_active_zone_counts[p]);
      vector<double> a = vector<double>(_active_zone_counts[p]);
      convert_matrix_to_coo_for_fortran_vector(_rate_matrix_transposed[p], ia,
                                               ja, a);
      _ia_s.push_back(ia);
      _ja_s.push_back(ja);
      _a_s.push_back(a);
    }
  }

  if (VERBOSE) {
    cout << "Q" << endl;
    for (unsigned int i = 0; i < _rate_matrix.size(); i++) {
      for (unsigned int j = 0; j < _rate_matrix[i].columns(); j++) {
        for (unsigned int k = 0; k < _rate_matrix[i].rows(); k++) {
          cout << _rate_matrix[i](j, k) << " ";
        }
        cout << endl;
      }
      cout << endl;
    }
  }
}

extern "C" {
void wrapalldmexpv_(int *n, int *m, double *t, double *v, double *w,
                    double *tol, double *anorm, double *wsp, int *lwsp,
                    int *iwsp, int *liwsp, int *itrace, int *iflag, int *ia,
                    int *ja, double *a, int *nz, double *res);
void wrapsingledmexpv_(int *n, int *m, double *t, double *v, double *w,
                       double *tol, double *anorm, double *wsp, int *lwsp,
                       int *iwsp, int *liwsp, int *itrace, int *iflag, int *ia,
                       int *ja, double *a, int *nz, double *res);
void wrapdgpadm_(int *ideg, int *m, double *t, double *H, int *ldh, double *wsp,
                 int *lwsp, int *ipiv, int *iexph, int *ns, int *iflag);
}

vector<vector<double>> RateModel::setup_fortran_P(int period, double t,
                                                  bool store_p_matrices) {
  auto &cur_rate_matrix = _rate_matrix[period];
  size_t row_size = cur_rate_matrix.rows();

  lagrange_matrix_t prob_matrix(cur_rate_matrix);
  prob_matrix *= t;

  _expm_count += 1;
  lagrange_matrix_t b = compute_matrix_exponential_ss(prob_matrix);
  /* we need to normalize the matrix */
  for (size_t i = 0; i < row_size; ++i) {
    double sum = 0.0;
    for (size_t j = 0; j < row_size; ++j) {
      sum += b(i, j);
    }
    for (size_t j = 0; j < row_size; ++j) {
      b(i, j) /= sum;
    }
  }

  vector<vector<double>> ret_vector;
  ret_vector.reserve(row_size);
  for (size_t i = 0; i < row_size; ++i) {
    ret_vector.emplace_back(row_size);
  }
  for (size_t i = 0; i < row_size; ++i) {
    for (size_t j = 0; j < row_size; ++j) {
      ret_vector[i][j] = b(i, j);
    }
  }

  if (store_p_matrices == true) {
    stored_p_matrices[period][t] = ret_vector;
  }

  return ret_vector;
}

/* TODO: Actually compute Horner's method for some particular values of q
 */
lagrange_matrix_t RateModel::compute_matrix_exponential_ss(
    lagrange_matrix_t A) const {
  size_t rows = A.rows();
  int scale_exp = std::max(0, 1 + static_cast<int>(blaze::linfNorm(A)));
  A /= std::pow(2.0, scale_exp);
  // q is a magic parameter that controls the number of iterations of the loop
  // higher is more accurate, with each increase of q decreasing error by 4
  // orders of magnitude. Anything above 12 is probably snake oil.
  constexpr int q = 3;
  double c = 0.5;
  double sign = -1.0;
  blaze::IdentityMatrix<double, blaze::columnMajor> I(rows);
  lagrange_matrix_t X = A;
  lagrange_matrix_t N = I + c * A;
  lagrange_matrix_t D = I - c * A;

  // Using fortran indexing, and we started an iteration ahead to skip some
  // setup
  for (int i = 2; i <= q; ++i) {
    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    X = A * X;
    N += c * X;
    sign *= -1.0;
    D += sign * c * X;
  }
  A = blaze::inv(D) * N;
  for (int i = 0; i < scale_exp; ++i) {
    A *= A;
  }
  return A;
}

lagrange_matrix_t RateModel::compute_matrix_exponential_eigen(
    const lagrange_matrix_t &A) const {
  size_t rows = A.rows();
  lagrange_complex_col_vector_t eigenvalues(rows);
  lagrange_complex_matrix_t eigenvectors(rows, rows);

  blaze::eigen(A, eigenvalues, eigenvectors);
  lagrange_complex_matrix_t inverse_eigenvectors = blaze::inv(eigenvectors);
  eigenvalues = blaze::exp(eigenvalues);

  return blaze::real(
      (eigenvectors % blaze::expand(blaze::trans(eigenvalues), rows)) *
      inverse_eigenvectors);
}

#if 0
/*
 * runs the basic padm fortran expokit full matrix exp
 */
vector<vector<double>> RateModel::setup_fortran_P(int period, double t,
                                                  bool store_p_matrices) {
  /*
  return P, the matrix of dist-to-dist transition probabilities,
  from the model's rate matrix (Q) over a time duration (t)
  */
  int ideg = 6;
  int rate_matrix_size = _rate_matrix[period].size();
  int ldh = rate_matrix_size;
  double tolerance = 1;
  int iflag = 0;
  int lwsp = 4 * rate_matrix_size * rate_matrix_size + 6 + 1;
  double *wsp = new double[lwsp];
  int *ipiv = new int[rate_matrix_size];
  int iexph = 0;
  int ns = 0;
  double *H = new double[rate_matrix_size * rate_matrix_size];
  convert_matrix_to_single_row_for_fortran(_rate_matrix[period], t, H);
  wrapdgpadm_(&ideg, &rate_matrix_size, &tolerance, H, &ldh, wsp, &lwsp, ipiv,
              &iexph, &ns, &iflag);
  vector<vector<double>> prob_matrix(
      _rate_matrix[period].size(), vector<double>(_rate_matrix[period].size()));
  for (int i = 0; i < rate_matrix_size; i++) {
    for (int j = 0; j < rate_matrix_size; j++) {
      prob_matrix[i][j] =
          wsp[iexph + (j - 1) * rate_matrix_size + (i - 1) + rate_matrix_size];
    }
  }
  delete[] wsp;
  delete[] ipiv;
  delete[] H;
  for (unsigned int i = 0; i < prob_matrix.size(); i++) {
    double sum = 0.0;
    for (unsigned int j = 0; j < prob_matrix[i].size(); j++) {
      sum += prob_matrix[i][j];
    }
    for (unsigned int j = 0; j < prob_matrix[i].size(); j++) {
      prob_matrix[i][j] = (prob_matrix[i][j] / sum);
    }
  }

  /*
   if store_p_matrices we will store them
   */
  if (store_p_matrices == true) {
    stored_p_matrices[period][t] = prob_matrix;
  }

  if (VERBOSE) {
    cout << "p " << period << " " << t << endl;
    for (unsigned int i = 0; i < prob_matrix.size(); i++) {
      for (unsigned int j = 0; j < prob_matrix[i].size(); j++) {
        cout << prob_matrix[i][j] << " ";
      }
      cout << endl;
    }
  }
  return prob_matrix;
}
#endif

/*
 * runs the sparse matrix fortran expokit matrix exp
 */
vector<vector<double>> RateModel::setup_sparse_full_P(int period, double t) {
  int n = _rate_matrix[period].columns();
  int m = _rate_matrix[period].columns() - 1;  // tweak
  int nz = get_size_for_coo(_rate_matrix[period]);
  int *ia = new int[nz];
  int *ja = new int[nz];
  double *a = new double[nz];
  convert_matrix_to_coo_for_fortran(_rate_matrix[period], ia, ja, a);
  double *v = new double[n];
  for (int i = 0; i < n; i++) {
    v[i] = 0;
  }
  v[0] = 1;
  double *w = new double[n];
  int ideg = 6;
  double tol = 1;
  int iflag = 0;
  int lwsp =
      n * (m + 1) + n + pow((m + 2.), 2) + 4 * pow((m + 2.), 2) + ideg + 1;
  double *wsp = new double[lwsp];
  int liwsp = m + 2;
  int *iwsp = new int[liwsp];
  double t1 = 1;
  double anorm = 0;
  int itrace = 0;
  double *res = new double[n * n];
  wrapalldmexpv_(&n, &m, &t1, v, w, &tol, &anorm, wsp, &lwsp, iwsp, &liwsp,
                 &itrace, &iflag, ia, ja, a, &nz, res);

  vector<vector<double>> p(_rate_matrix[period].rows(),
                           vector<double>(_rate_matrix[period].rows()));

  double *res_ptr = res;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      p[j][i] = *(res_ptr++);
    }
  }

  // filter out impossible dists
  for (unsigned int i = 0; i < _dists.size(); i++) {
    if (_dists[i] > 0) {
      for (unsigned int j = 0; j < _area_count; j++) {
        if (lagrange_bextr(_dists[i], j) == 0) {  // present
          continue;
        }
        double sum1 =
            calculate_vector_double_sum(_dispersal_params_mask[period][j]);
        double sum2 = 0.0;
        for (unsigned int k = 0; k < _dispersal_params_mask[period].size();
             k++) {
          sum2 += _dispersal_params_mask[period][k][j];
        }
        if (sum1 + sum2 == 0) {
          for (unsigned int k = 0; k < p[period].size(); k++) {
            p[period][k] = p[period][k] * 0.0;
          }
          break;
        }
      }
    }
  }
  delete[] res;
  delete[] iwsp;
  delete[] w;
  delete[] wsp;
  delete[] a;
  delete[] ia;
  delete[] ja;
  if (VERBOSE) {
    cout << "p " << period << " " << t << endl;
    for (unsigned int i = 0; i < p.size(); i++) {
      for (unsigned int j = 0; j < p[i].size(); j++) {
        cout << p[i][j] << " ";
      }
      cout << endl;
    }
  }
  return p;
}

/*
 * for returning single P columns
 */
vector<double> RateModel::setup_sparse_single_column_P(int period, double t,
                                                       int column) {
  int n = _rate_matrix[period].rows();
  int m = _area_count - 1;
  int current_zone_count =
      _active_zone_counts[period];  // get_size_for_coo(Q[period],1);
  int *ia = new int[current_zone_count];
  int *ja = new int[current_zone_count];
  double *a = new double[current_zone_count];
  std::copy(_ia_s[period].begin(), _ia_s[period].end(), ia);
  std::copy(_ja_s[period].begin(), _ja_s[period].end(), ja);
  std::copy(_a_s[period].begin(), _a_s[period].end(), a);
  double *v = new double[n];
  for (int i = 0; i < n; i++) {
    v[i] = 0;
  }
  v[column] = 1;  // only return the one column we want
  double *w = new double[n];
  int ideg = 6;
  double tol = 1;
  int iflag = 0;
  int lwsp =
      n * (m + 1) + n + pow((m + 2.), 2) + 4 * pow((m + 2.), 2) + ideg + 1;
  double *wsp = new double[lwsp];
  int liwsp = m + 2;
  int *iwsp = new int[liwsp];
  double t1 = t;  // use to be 1
  double anorm = 0;
  int itrace = 0;
  double *res = new double[n];  // only needs resulting columns
  wrapsingledmexpv_(&n, &m, &t1, v, w, &tol, &anorm, wsp, &lwsp, iwsp, &liwsp,
                    &itrace, &iflag, ia, ja, a, &current_zone_count, res);

  vector<double> p(_rate_matrix[period].rows());
  int count = 0;
  for (int i = 0; i < n; i++) {
    p[i] = res[count];
    count += 1;
  }
  delete[] v;
  delete[] res;
  delete[] iwsp;
  delete[] w;
  delete[] wsp;
  delete[] a;
  delete[] ia;
  delete[] ja;
  if (VERBOSE) {
    cout << "p " << period << " " << t << " " << column << endl;
    for (unsigned int i = 0; i < p.size(); i++) {
      cout << p[i] << " ";
    }
    cout << endl;
  }
  return p;
}

vector<lagrange_region_split_t> RateModel::iter_dist_splits(
    lagrange_dist_t dist) {
  vector<lagrange_region_split_t> ret;

  if (dist == 0) {
    return ret;
  }

  if (lagrange_popcount(dist) == 1) {
    ret.push_back({dist, dist});
    return ret;
  }

  for (unsigned int i = 0; i < _area_count; i++) {
    if (lagrange_bextr(dist, i) == 0) {
      continue;
    }

    lagrange_dist_t x = 1ull << i;
    auto iter = std::find(_dists.begin(), _dists.end(), x);

    if (iter == _dists.end()) {
      continue;
    }

    ret.push_back({x, dist});
    ret.push_back({dist, x});

    lagrange_dist_t y = (x ^ dist) & _valid_dist_mask;
    iter = std::find(_dists.begin(), _dists.end(), y);
    if (iter == _dists.end()) {
      continue;
    }

    ret.push_back({x, y});

    if (lagrange_popcount(y) > 1) {
      ret.push_back({y, x});
    }
  }
  return ret;
}

void RateModel::iter_all_dist_splits() {
  for (auto d : _dists) {
    _iter_dists[d] = iter_dist_splits(d);
  }
}

const vector<lagrange_dist_t> &RateModel::getDists() { return _dists; }

size_t RateModel::getDistsSize() const { return _dists.size(); }

const vector<lagrange_region_split_t> &RateModel::get_iter_dist_splits(
    lagrange_dist_t dist) const {
  return _iter_dists.at(dist);
}

size_t RateModel::get_num_areas() { return _area_count; }

int RateModel::get_num_periods() { return _periods.size(); }

vector<lagrange_matrix_t> &RateModel::get_Q() { return _rate_matrix; }

inline int signof(double d) { return d >= 0 ? 1 : -1; }

inline double roundto(double in) { return floor(in * (1000) + 0.5) / (1000); }

/*
 * this should be used to caluculate the eigenvalues and eigenvectors
 * as U * Q * U-1 -- eigen decomposition
 */
bool RateModel::get_eigenvec_eigenval_from_Q(lagrange_complex_matrix_t &eigval,
                                             lagrange_complex_matrix_t &eigvec,
                                             int period) {
  size_t row_size = _rate_matrix[period].rows();
  lagrange_matrix_t tQ(row_size, row_size, 0.0);

  for (unsigned int i = 0; i < _rate_matrix[period].rows(); i++) {
    for (unsigned int j = 0; j < _rate_matrix[period].columns(); j++) {
      tQ(i, j) = _rate_matrix[period](i, j);
    }
  }

  lagrange_complex_col_vector_t eigva;
  lagrange_complex_matrix_t eigve;
  blaze::eigen(tQ, eigva, eigve);

  bool isImag = false;
  for (unsigned int i = 0; i < _rate_matrix[period].rows(); i++) {
    for (unsigned int j = 0; j < _rate_matrix[period].columns(); j++) {
      if (i == j) {
        eigval(i, j) = eigva[i];
      } else {
        eigval(i, j) = 0;
      }
      eigvec(i, j) = eigve(i, j);
      if (imag(eigvec(i, j)) > 0 || imag(eigval(i, j))) isImag = true;
    }
  }
  if (VERBOSE) {
    cout << eigva << endl;
    cout << tQ - ((eigvec) * (eigval)*blaze::inv(eigvec)) << endl;
  }
  return isImag;
}

size_t RateModel::get_expm_count() { return _expm_count; }

vector<AncSplit> RateModel::iter_ancsplits(lagrange_dist_t dist) {
  vector<AncSplit> ans;
  auto splits = get_iter_dist_splits(dist);
  if (splits.size() == 0) {
    return ans;
  }
  int nsplits = splits.size();
  double weight = 1.0 / nsplits;
  for (unsigned int i = 0; i < splits.size(); i++) {
    AncSplit an(dist, splits[i].left, splits[i].right, weight);
    ans.push_back(an);
  }
  return ans;
}

void RateModel::iter_ancsplits_just_int(lagrange_dist_t &dist,
                                        vector<int> &leftdists,
                                        vector<int> &rightdists,
                                        double &weight) {
  leftdists.clear();
  rightdists.clear();

  auto splits = get_iter_dist_splits(dist);
  auto distsmap = get_dists_int_map();
  size_t dists_size = splits.size() * splits.size();

  leftdists.reserve(dists_size);
  rightdists.reserve(dists_size);

  if (splits.size() == 0) {
    return;
  }
  int nsplits = splits.size();
  weight = 1.0 / nsplits;
  for (unsigned int i = 0; i < splits.size(); i++) {
    leftdists.push_back(distsmap[splits[i].left]);
    rightdists.push_back(distsmap[splits[i].right]);
  }
}

vector<int> RateModel::get_columns_for_sparse(vector<double> &inc) {
  vector<int> ret(inc.size(), 0);
  for (unsigned int i = 0; i < inc.size(); i++) {
    if (inc[i] < 0.0000000001) {
      continue;
    }
    ret[i] = 1;
    auto dis = getDists().at(i);
    for (unsigned int j = 0; j < inc.size(); j++) {
      auto dis2 = getDists().at(j);
      if (lagrange_popcount(dis ^ dis2) == 1) {
        ret[j] = 1;
      }
    }
  }
  return ret;
}

vector<int> RateModel::get_columns_for_sparse(vector<Superdouble> &inc) {
  vector<int> ret(inc.size(), 0);
  for (unsigned int i = 0; i < inc.size(); i++) {
    if (inc[i] < Superdouble(0.0000000001)) {
      continue;
    }
    ret[i] = 1;
    auto dis = getDists().at(i);
    for (unsigned int j = 0; j < inc.size(); j++) {
      auto dis2 = getDists().at(j);
      if (lagrange_popcount(dis ^ dis2) == 1) {
        ret[j] = 1;
      }
    }
  }
  return ret;
}

const unordered_map<lagrange_dist_t, int> &RateModel::get_dists_int_map() {
  return _dists_int_map;
}

const unordered_map<int, lagrange_dist_t> &RateModel::get_int_dists_map() {
  return _int_dists_map;
}
