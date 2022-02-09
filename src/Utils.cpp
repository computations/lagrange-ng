/*
 * Utils.cpp
 *
 *  Created on: Mar 10, 2009
 *      Author: Stephen A. Smith
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include "Common.h"
#include "Utils.h"

void Tokenize(const std::string &str, std::vector<std::string> &tokens,
              const std::string &delimiters) {
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void TrimSpaces(std::string &str) {
  // Trim Both leading and trailing spaces
  size_t startpos =
      str.find_first_not_of(" \t\r\n");  // Find the first character position
                                         // after excluding leading blank spaces
  size_t endpos = str.find_last_not_of(
      " \t\r\n");  // Find the first character position from reverse af

  // if all spaces or empty return an empty std::string
  if ((std::string::npos == startpos) || (std::string::npos == endpos)) {
    str = "";
  } else
    str = str.substr(startpos, endpos - startpos + 1);
}

std::string lagrange_convert_dist_string(
    lagrange_dist_t dist, const std::vector<std::string> &names) {
  if (dist == 0) { return {}; }
  std::ostringstream oss;

  size_t states = lagrange_clz(dist);
  bool first = true;
  for (size_t i = 0; i < states; ++i) {
    if (lagrange_bextr(dist, i)) {
      if (!first) { oss << "_"; }
      oss << names[i];
      first = false;
    }
  }
  return oss.str();
}

std::vector<std::string> lagrange_convert_dist_to_list(
    lagrange_dist_t dist, const std::vector<std::string> &names) {
  std::vector<std::string> ret;

  size_t states = lagrange_clz(dist);
  for (size_t i = 0; i < states; ++i) {
    if (lagrange_bextr(dist, i)) { ret.push_back(names[i]); }
  }
  return ret;
}

/* Code for the combinations function take from phylourny with permision from
 * author.*/

constexpr size_t factorial_table_size = 11;

constexpr std::array<double, factorial_table_size> factorial_table = {
    1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800,
};

constexpr inline double factorial(uint64_t i) {
  if (i < factorial_table_size) { return factorial_table.at(i); }
  double f = factorial_table[factorial_table_size - 1];
  for (size_t k = factorial_table_size; k <= i; ++k) {
    f *= static_cast<double>(k);
  }
  return f;
}

constexpr inline double combinations(uint64_t n, uint64_t i) {
  return factorial(n) / (factorial(i) * factorial(n - i));
}

size_t lagrange_compute_restricted_state_count(size_t regions,
                                               size_t max_areas) {
  size_t sum = 0;
  for (size_t i = 0; i <= max_areas; i++) { sum += combinations(regions, i); }
  return sum;
}

size_t compute_skips_power_of_2(size_t k, size_t n) {
  size_t skips = 0;
  for (size_t i = n + 1; i < k; ++i) { skips += combinations(k - 1, i); }
  return skips;
}

size_t compute_skips(size_t i, size_t n) {
  size_t skips = 0;
  while (i != 0 && n != 0) {
    size_t first_index = sizeof(i) * 8 - __builtin_clzll(i | 1);
    skips += compute_skips_power_of_2(first_index, n);
    n -= 1;
    i -= 1 << (first_index - 1);
  }
  skips += i;
  return skips;
}

size_t compute_index_from_dist(lagrange_dist_t i, size_t max_areas) {
  if (lagrange_popcount(i) > max_areas) {
    throw std::lagrange_util_dist_index_conversion_exception{};
  }
  size_t skips = compute_skips(i, max_areas);

  return i - skips;
}
