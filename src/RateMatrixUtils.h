/*
 * RateMatrixUtils.h
 *
 *  Created on: Aug 13, 2009
 *      Author: Stephen A. Smith
 */

/*
  This is essentially a utility file that stores matrix utility functions
  dealing with the rate matrix.
 */

#ifndef RATEMATRIXUTILS_H_
#define RATEMATRIXUTILS_H_
#include "AncSplit.h"
#include "superdouble.h"
#include <memory>
#include <vector>
using namespace std;

#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif

/*
  most of these utilities calculate simple math on matrices and vectors
  and they should only be used if there is the chance that there will
  be a null vector or matrix. otherwise, c++ numerics library should
  be used for speed.
 */
#ifdef BIGTREE
double calculate_vector_mpfr_class_double_sum(vector<mpfr_class> &in);
#endif
double calculate_vector_double_sum(vector<double> &in);
Superdouble calculate_vector_Superdouble_sum(const vector<Superdouble> &in);
int calculate_vector_int_sum(const vector<int> &in);
int get_vector_int_index_from_multi_vector_int(const vector<int> &in,
                                               const vector<vector<int>> &in2);

/*
  used for creating and enumerate distributions
 */
int calculate_vector_int_sum_xor(vector<int> &in, vector<int> &in2);
int locate_vector_int_single_xor(vector<int> &in, vector<int> &in2);

/*
  simple printing functions
 */
void print_vector_int(vector<int> &in);
void print_vector_double(vector<double> &in);

/*
  used for generating all the distributions with maximium number of areas
  involved would be designated in the config file
 */
vector<vector<int>> generate_dists_from_num_max_areas(int totalnum,
                                                      int numareas);

/*
  used for processing custom rate matrix config files designated in the main
  config file
 */
vector<vector<vector<double>>>
processRateMatrixConfigFile(string filename, int numareas, int nperiods);

#endif
