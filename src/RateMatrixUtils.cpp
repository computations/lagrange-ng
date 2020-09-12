/*
 * RateMatrixUtils.cpp
 *
 *  Created on: Aug 13, 2009
 * Author: smitty
 */

#include "RateMatrixUtils.h"
#include "Utils.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
using namespace std;

#ifdef BIGTREE
double calculate_vector_mpfr_class_double_sum(vector<mpfr_class> &in) {
  double sum = 0;
  for (unsigned int i = 0; i < in.size(); i++) {
    double x = in[i].get_d();
    sum += x;
  }
  return sum;
}
#endif

double calculate_vector_double_sum(vector<double> &in) {
  double sum = 0;
  for (unsigned int i = 0; i < in.size(); i++) {
    sum += in[i];
  }
  return sum;
}

Superdouble calculate_vector_Superdouble_sum(vector<Superdouble> &in) {
  Superdouble sum = 0;
  for (unsigned int i = 0; i < in.size(); i++) {
    sum += in[i];
    // cout << in[i] << " sum:" << sum << endl;
  }
  // cout << "endsum:" << sum << endl;
  return sum;
}

/*
 * only used because sometimes will send a null
 */
int calculate_vector_int_sum(const vector<int> &in) {
  int sum = 0;
  for (unsigned int i = 0; i < in.size(); i++) {
    sum += in.at(i);
  }
  return sum;
}

int calculate_vector_int_sum_xor(vector<int> &in, vector<int> &in2) {
  int sum = 0;
  for (unsigned int i = 0; i < in.size(); i++) {
    if (in[i] != in2[i]) {
      sum += 1;
    }
  }
  return sum;
}

int locate_vector_int_single_xor(vector<int> &in, vector<int> &in2) {
  int location = 0;
  for (unsigned int i = 0; i < in.size(); i++) {
    if (in[i] != in2[i]) {
      location = i;
    }
  }
  return location;
}

int get_vector_int_index_from_multi_vector_int(const vector<int> &in,
                                               const vector<vector<int>> &in2) {
  int ret = 0;
  for (unsigned int i = 0; i < in2.size(); i++) {
    size_t sum = 0;
    if (in2[i].size() != in.size()) {
      throw std::invalid_argument{"Input vectors are not of compatible sizes"};
    }
    for (unsigned int j = 0; j < in2[i].size(); j++) {
      if (in2[i][j] == in[j])
        sum += 1;
    }
    if (sum == in.size()) {
      ret = i;
      return ret;
    }
  }
  string dstring;
  for (unsigned int j = 0; j < in.size(); j++) {
    stringstream ss;
    ss << in.at(j);
    dstring.append(ss.str());
  }
  throw std::runtime_error{std::string("the distribution ") + dstring +
                           +" is not included in the possible distributions"};
}

vector<vector<int>> generate_dists_from_num_max_areas(int totalnumareas,
                                                      int numareas) {
  vector<vector<int>> dists;
  map<int, vector<int>> a = iterate_all_bv(totalnumareas);
  // globalextinction
  vector<int> empt;
  for (unsigned int i = 0; i < a[0].size(); i++) {
    empt.push_back(0);
  }
  dists.push_back(empt);

  map<int, vector<int>>::iterator pos;
  for (pos = a.begin(); pos != a.end(); ++pos) {
    int f = pos->first;
    // if(calculate_vector_int_sum(&a[f]) <= numareas)
    if (accumulate(a[f].begin(), a[f].end(), 0) <= numareas)
      dists.push_back(a[f]);
  }
  return dists;
}

/*
 TODO: change this to store to memory instead of creating them
 */
vector<AncSplit> iter_ancsplits(std::shared_ptr<RateModel> rm,
                                vector<int> &dist) {
  vector<AncSplit> ans;
  auto splits = rm->get_iter_dist_splits(dist);
  auto distsmap = rm->get_dists_int_map();
  if (splits.at(0).size() > 0) {
    int nsplits = splits.at(0).size();
    double weight = 1.0 / nsplits;
    for (unsigned int i = 0; i < splits.at(0).size(); i++) {
      AncSplit an(rm, distsmap[dist], distsmap[splits.at(0)[i]],
                  distsmap[splits.at(1)[i]], weight);
      ans.push_back(an);
    }
  }
  return ans;
}

void iter_ancsplits_just_int(std::shared_ptr<RateModel> rm, vector<int> &dist,
                             vector<int> &leftdists, vector<int> &rightdists,
                             double &weight) {
  leftdists.clear();
  rightdists.clear();
  auto splits = rm->get_iter_dist_splits(dist);
  auto distsmap = rm->get_dists_int_map();
  size_t dists_size = splits[0].size() * splits[0][0].size();
  leftdists.reserve(dists_size);
  rightdists.reserve(dists_size);
  if (splits.at(0).size() > 0) {
    int nsplits = splits.at(0).size();
    weight = 1.0 / nsplits;
    for (unsigned int i = 0; i < splits.at(0).size(); i++) {
      leftdists.push_back(distsmap[splits.at(0)[i]]);
      rightdists.push_back(distsmap[splits.at(1)[i]]);
    }
  }
}

void print_vector_int(vector<int> &in) {
  for (unsigned int i = 0; i < in.size(); i++) {
    cout << in[i] << " ";
  }
  cout << endl;
}

void print_vector_double(vector<double> &in) {
  for (unsigned int i = 0; i < in.size(); i++) {
    cout << in[i] << " ";
  }
  cout << endl;
}

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

vector<vector<vector<double>>>
processRateMatrixConfigFile(string filename, int numareas, int nperiods) {
  vector<double> cols(numareas, 1);
  vector<vector<double>> rows(numareas, cols);
  vector<vector<vector<double>>> ratematrix =
      vector<vector<vector<double>>>(nperiods, rows);
  // read file
  ifstream ifs(filename.c_str());
  string line;
  int period = 0;
  int fromarea = 0;
  while (getline(ifs, line)) {
    if (line.size() > 3) {
      vector<string> tokens;
      string del(" ,\t");
      tokens.clear();
      Tokenize(line, tokens, del);
      for (unsigned int j = 0; j < tokens.size(); j++) {
        TrimSpaces(tokens[j]);
      }
      for (unsigned int j = 0; j < tokens.size(); j++) {
        ratematrix[period][fromarea][j] = atof(tokens[j].c_str());
      }
      if (fromarea < numareas - 1) {
        fromarea += 1;
      } else {
        fromarea = 0;
        period += 1;
      }
    }
  }
  ifs.close();
  return ratematrix;
}

/*
 * need to make this much faster
 */
vector<int> get_columns_for_sparse(vector<double> &inc,
                                   std::shared_ptr<RateModel> rm) {
  vector<int> ret(inc.size(), 0);
  for (unsigned int i = 0; i < inc.size(); i++) {
    if (inc[i] > 0.0000000001) {
      ret[i] = 1;
      vector<int> dis = rm->getDists().at(i);
      for (unsigned int j = 0; j < inc.size(); j++) {
        vector<int> dis2 = rm->getDists().at(j);
        int sum = calculate_vector_int_sum_xor(dis, dis2);
        if (sum == 1) {
          ret[j] = 1;
        }
      }
    }
  }
  return ret;
}

vector<int> get_columns_for_sparse(vector<Superdouble> &inc,
                                   std::shared_ptr<RateModel> rm) {
  vector<int> ret(inc.size(), 0);
  for (unsigned int i = 0; i < inc.size(); i++) {
    if (inc[i] > Superdouble(0.0000000001)) {
      ret[i] = 1;
      vector<int> dis = rm->getDists().at(i);
      for (unsigned int j = 0; j < inc.size(); j++) {
        vector<int> dis2 = rm->getDists().at(j);
        int sum = calculate_vector_int_sum_xor(dis, dis2);
        if (sum == 1) {
          ret[j] = 1;
        }
      }
    }
  }
  return ret;
}
