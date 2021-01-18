/*
 * Utils.h
 *
 *  Created on: Mar 10, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <map>
#include <string>
#include <vector>

#include "Common.h"

using namespace std;

struct lagrange_region_split_t {
  lagrange_dist_t left;
  lagrange_dist_t right;
};

void Tokenize(const string &str, vector<string> &tokens,
              const string &delimiters = " ");
void TrimSpaces(string &str);
long comb(int m, int n);
vector<vector<int>> iterate(int M, int N);
vector<int> comb_at_index(int m, int n, int i);
vector<vector<int>> dists_by_maxsize(int nareas, int maxsize);
vector<int> idx2bitvect(vector<int> indices, int M);
vector<vector<int>> iterate_all(int m);
map<int, lagrange_dist_t> iterate_all_bv(int m);
map<int, lagrange_dist_t> iterate_all_bv2(int m);

inline uint64_t lagrange_bextr(lagrange_dist_t a, size_t i) {
  return (a >> i) & 1ull;
}

inline size_t lagrange_popcount(lagrange_dist_t a) {
  return __builtin_popcountll(a);
}

inline size_t lagrange_ctz(lagrange_dist_t a) { return __builtin_ctzll(a); }

inline lagrange_dist_t convert_vector_to_lagrange_dist(
    const vector<int> &vec_dist) {
  lagrange_dist_t ret = 0;
  size_t area_count = vec_dist.size();

  for (size_t i = 0; i < area_count; ++i) {
    ret |= static_cast<lagrange_dist_t>(vec_dist[i]) << i;
  }

  return ret;
}

constexpr inline size_t lagrange_fast_log2(size_t x) {
  return __builtin_clzll(x);
}

#endif /* UTILS_H_ */
