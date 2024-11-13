/*
 * Utils.cpp
 *
 *  Created on: Mar 10, 2009
 *      Author: Stephen A. Smith
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include "Utils.hpp"

#include <array>

#include "Common.hpp"

namespace lagrange {

auto lagrange_convert_dist_string(
    Range dist, const std::vector<std::string> &names) -> std::string {
  if (dist == 0) { return {}; }
  std::ostringstream oss;

  size_t states = lagrange_clz(dist);
  bool first = true;
  for (size_t i = 0; i < states; ++i) {
    if (lagrange_bextr(dist, i) != 0U) {
      if (!first) { oss << "_"; }
      oss << names[i];
      first = false;
    }
  }
  return oss.str();
}

std::vector<std::string> lagrange_convert_dist_to_list(
    Range dist, const std::vector<std::string> &names) {
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

constexpr std::array<size_t, factorial_table_size> factorial_table = {
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5040,
    40320,
    362880,
    3628800,
};

constexpr inline auto factorial(uint64_t i) -> size_t {
  if (i < factorial_table_size) { return factorial_table.at(i); }
  size_t f = factorial_table[factorial_table_size - 1];
  for (size_t k = factorial_table_size; k <= i; ++k) {
    f *= static_cast<double>(k);
  }
  return f;
}

constexpr inline auto combinations(uint64_t n, uint64_t i) -> size_t {
  return factorial(n) / (factorial(i) * factorial(n - i));
}

auto lagrange_compute_restricted_state_count(RangeSize regions,
                                             RangeSize max_areas) -> RangeSize {
  RangeSize sum = 0;
  for (RangeSize i = 0; i <= max_areas; i++) {
    sum += combinations(regions, i);
  }
  return sum;
}

static auto compute_skips_power_of_2(size_t k, size_t n) -> size_t {
  size_t skips = 0;
  for (size_t i = n + 1; i < k; ++i) { skips += combinations(k - 1, i); }
  return skips;
}

static auto compute_skips(size_t i, size_t n) -> size_t {
  constexpr size_t BITS_IN_BYTE = 8;
  size_t skips = 0;
  while (i != 0 && n != 0) {
    size_t first_index = sizeof(i) * BITS_IN_BYTE - __builtin_clzll(i | 1);
    skips += compute_skips_power_of_2(first_index, n);
    n -= 1;
    i -= 1 << (first_index - 1);
  }
  skips += i;
  return skips;
}

auto compute_index_from_dist(Range i, RangeSize max_areas) -> RangeSize {
  if (lagrange_popcount(i) > max_areas) {
    throw UtilDistIndexConversionException{};
  }
  size_t skips = compute_skips(i, max_areas);

  return i - skips;
}

auto lagrange_parse_size_t(const std::string &str) -> size_t {
  int temp = stoi(str);
  if (temp < 0) {
    throw std::invalid_argument{"This argument should be positive"};
  }
  return static_cast<size_t>(temp);
}

auto convert_dist_string_to_dist(
    const std::string &dist, const std::vector<std::string> &names) -> Range {
  Range ret = 0;
  auto start = dist.begin();
  auto end = dist.begin();
  for (start = dist.begin(), end = dist.begin(); end != dist.end(); end++) {
    if (*end != '_') { continue; }

    std::string region_name(start, end - 1);
    end++;
    start = end;
    for (size_t i = 0; i < names.size(); ++i) {
      if (region_name == names[i]) {
        ret |= 1ull << i;
        break;
      }
    }
  }
  return ret;
}

auto convert_dist_binary_string_to_dist(const std::string &dist) -> Range {
  Range d = 0;
  for (RangeIndex i = 0; i < dist.size(); i++) {
    if (dist[i] == '1') { d |= 1ull << (dist.size() - i - 1); }
  }
  return d;
}

}  // namespace lagrange
