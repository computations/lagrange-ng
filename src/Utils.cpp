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
#include <ranges>

#include "Common.hpp"

namespace lagrange {

auto lagrange_convert_dist_string(Range dist,
                                  const std::vector<std::string> &names)
    -> std::string {
  if (dist == 0) { return {}; }
  std::ostringstream oss;

  size_t states = lagrange_clz(dist);
  auto itr = names.begin();
  bool first = true;
  for (size_t i = 0; i < states; ++i) {
    if (lagrange_bextr(dist, i) != 0U) {
      if (!first) { oss << "_"; }
      oss << *itr;
      first = false;
    }
    itr++;
  }
  return oss.str();
}

auto lagrange_convert_dist_to_list(Range dist,
                                   const std::vector<std::string> &names)
    -> std::vector<std::string> {
  std::vector<std::string> ret;

  size_t states = dist == 0 ? 0 : lagrange_clz(dist);
  auto itr = names.begin();
  for (size_t i = 0; i < states; ++i) {
    if (lagrange_bextr(dist, i) != 0u) { ret.push_back(*itr); }
    itr++;
  }
  return ret;
}

auto lagrange_convert_list_to_dist(const std::vector<std::string> &dist,
                                   const std::vector<std::string> &names)
    -> Range {
  Range ret = 0;
  for (const auto &area : dist) {
#ifdef __cpp_lib_ranges_enumerate
    for (const auto &[i, n] : std::views::enumerate(names)) {
      if (area == n) {
#else
    for (size_t i = 0; i < names.size(); ++i) {
      if (area == names[i]) {
#endif
        ret |= 1ULL << i;
        break;
      }
    }
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

constexpr auto factorial(uint64_t i) -> size_t {
  if (i < factorial_table_size) { return factorial_table.at(i); }
  size_t f = factorial_table[factorial_table_size - 1];
  for (size_t k = factorial_table_size; k <= i; ++k) { f *= k; }
  return f;
}

constexpr auto combinations(uint64_t n, uint64_t i) -> size_t {
  return factorial(n) / (factorial(i) * factorial(n - i));
}

auto lagrange_compute_restricted_state_count(size_t regions, size_t max_areas)
    -> size_t {
  size_t sum = 0;
  for (size_t i = 0; i <= max_areas; i++) { sum += combinations(regions, i); }
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

auto compute_index_from_dist(Range i, size_t max_areas) -> size_t {
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

auto convert_dist_string_to_dist(const std::string &dist,
                                 const std::vector<std::string> &names)
    -> Range {
  Range ret = 0;
  auto start = dist.begin();
  auto end = dist.begin();
  for (start = dist.begin(), end = dist.begin(); end != dist.end(); end++) {
    if (*end != '_') { continue; }

    std::string_view region_name(start, end - 1);
    end++;
    start = end;
    for (size_t i = 0; i < names.size(); ++i) {
      if (region_name == names[i]) {
        ret |= 1ULL << i;
        break;
      }
    }
  }
  return ret;
}

auto convert_dist_binary_string_to_dist(const std::string &dist) -> Range {
  Range d = 0;
  for (size_t i = 0; i < dist.size(); i++) {
    if (dist[i] == '1') { d |= 1ULL << (dist.size() - i - 1); }
  }
  return d;
}

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)

  #include <windows.h>

size_t get_system_memory() {
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  GlobalMemoryStatusEx(&status);
  return status.ullTotalPhys;
}

#elif __unix__

  #include <unistd.h>

size_t get_system_memory() {
  auto page_count = sysconf(_SC_PHYS_PAGES);
  auto page_size = sysconf(_SC_PAGE_SIZE);
  return page_count * page_size;
}
#endif

}  // namespace lagrange
