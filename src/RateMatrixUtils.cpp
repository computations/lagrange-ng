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
