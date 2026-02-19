#include "Utils.hpp"

#include <stdexcept>
#include <vector>

#include "gtest/gtest.h"

using namespace lagrange;

TEST(Utils, convert_dist_binary_string_to_dist) {
  EXPECT_EQ(convert_dist_binary_string_to_dist("0"), 0b0);
  EXPECT_EQ(convert_dist_binary_string_to_dist("1"), 0b1);
  EXPECT_EQ(convert_dist_binary_string_to_dist("10"), 0b10);
  EXPECT_EQ(convert_dist_binary_string_to_dist("01"), 0b01);
  EXPECT_EQ(convert_dist_binary_string_to_dist("101"), 0b101);
  EXPECT_EQ(convert_dist_binary_string_to_dist("011"), 0b011);
  EXPECT_EQ(convert_dist_binary_string_to_dist("111"), 0b111);
  EXPECT_EQ(convert_dist_binary_string_to_dist("1000"), 0b1000);
}

TEST(Utils, convert_dist_string_to_dist) {
  std::vector<std::string> names = {"A", "B", "C"};
  EXPECT_EQ(convert_dist_string_to_dist("", names), 0);
}

TEST(Utils, lagrange_convert_dist_to_list) {
  std::vector<std::string> names = {"A", "B", "C"};

  auto result = lagrange_convert_dist_to_list(0b001, names);
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result[0], "A");

  result = lagrange_convert_dist_to_list(0b010, names);
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result[0], "B");

  result = lagrange_convert_dist_to_list(0b100, names);
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result[0], "C");

  result = lagrange_convert_dist_to_list(0b011, names);
  EXPECT_EQ(result.size(), 2);
  EXPECT_EQ(result[0], "A");
  EXPECT_EQ(result[1], "B");

  result = lagrange_convert_dist_to_list(0b111, names);
  EXPECT_EQ(result.size(), 3);
  EXPECT_EQ(result[0], "A");
  EXPECT_EQ(result[1], "B");
  EXPECT_EQ(result[2], "C");

  result = lagrange_convert_dist_to_list(0b000, names);
  EXPECT_EQ(result.size(), 0);
}

TEST(Utils, lagrange_convert_dist_string) {
  std::vector<std::string> names = {"A", "B", "C"};

  EXPECT_EQ(lagrange_convert_dist_string(0b001, names), "A");
  EXPECT_EQ(lagrange_convert_dist_string(0b010, names), "B");
  EXPECT_EQ(lagrange_convert_dist_string(0b100, names), "C");
  EXPECT_EQ(lagrange_convert_dist_string(0b011, names), "A_B");
  EXPECT_EQ(lagrange_convert_dist_string(0b101, names), "A_C");
  EXPECT_EQ(lagrange_convert_dist_string(0b110, names), "B_C");
  EXPECT_EQ(lagrange_convert_dist_string(0b111, names), "A_B_C");
  EXPECT_EQ(lagrange_convert_dist_string(0b000, names), "");
}

TEST(Utils, lagrange_convert_list_to_dist) {
  std::vector<std::string> names = {"A", "B", "C"};

  EXPECT_EQ(lagrange_convert_list_to_dist({"A"}, names), 0b001);
  EXPECT_EQ(lagrange_convert_list_to_dist({"B"}, names), 0b010);
  EXPECT_EQ(lagrange_convert_list_to_dist({"C"}, names), 0b100);
  EXPECT_EQ(lagrange_convert_list_to_dist({"A", "B"}, names), 0b011);
  EXPECT_EQ(lagrange_convert_list_to_dist({"A", "C"}, names), 0b101);
  EXPECT_EQ(lagrange_convert_list_to_dist({"B", "C"}, names), 0b110);
  EXPECT_EQ(lagrange_convert_list_to_dist({"A", "B", "C"}, names), 0b111);
  EXPECT_EQ(lagrange_convert_list_to_dist({}, names), 0b000);
}

TEST(Utils, lagrange_compute_restricted_state_count) {
  EXPECT_EQ(lagrange_compute_restricted_state_count(1, 1), 2);
  EXPECT_EQ(lagrange_compute_restricted_state_count(2, 1), 3);
  EXPECT_EQ(lagrange_compute_restricted_state_count(2, 2), 4);
  EXPECT_EQ(lagrange_compute_restricted_state_count(3, 1), 4);
  EXPECT_EQ(lagrange_compute_restricted_state_count(3, 2), 7);
  EXPECT_EQ(lagrange_compute_restricted_state_count(3, 3), 8);
  EXPECT_EQ(lagrange_compute_restricted_state_count(4, 1), 5);
  EXPECT_EQ(lagrange_compute_restricted_state_count(4, 2), 11);
  EXPECT_EQ(lagrange_compute_restricted_state_count(4, 3), 15);
  EXPECT_EQ(lagrange_compute_restricted_state_count(4, 4), 16);
}

TEST(Utils, compute_index_from_dist_basic) {
  EXPECT_EQ(compute_index_from_dist(0b0001, 4), 1);
  EXPECT_EQ(compute_index_from_dist(0b0010, 4), 2);
  EXPECT_EQ(compute_index_from_dist(0b0100, 4), 4);
  EXPECT_EQ(compute_index_from_dist(0b1000, 4), 8);

  EXPECT_THROW(compute_index_from_dist(0b11111, 4),
               UtilDistIndexConversionException);
}

TEST(Utils, compute_index_from_dist_max_areas) {
  EXPECT_THROW(compute_index_from_dist(0b0111, 2),
               UtilDistIndexConversionException);
  EXPECT_THROW(compute_index_from_dist(0b1111, 2),
               UtilDistIndexConversionException);
}

TEST(Utils, lagrange_parse_size_t) {
  EXPECT_EQ(lagrange_parse_size_t("0"), 0);
  EXPECT_EQ(lagrange_parse_size_t("1"), 1);
  EXPECT_EQ(lagrange_parse_size_t("100"), 100);
  EXPECT_EQ(lagrange_parse_size_t("999999"), 999999);

  EXPECT_THROW(lagrange_parse_size_t("-1"), std::invalid_argument);
  EXPECT_THROW(lagrange_parse_size_t("abc"), std::invalid_argument);
}

TEST(Utils, next_dist_basic) {
  EXPECT_EQ(next_dist(0b0000, 1), 0b0001);
  EXPECT_EQ(next_dist(0b0001, 1), 0b0010);
  EXPECT_EQ(next_dist(0b0010, 1), 0b0100);
  EXPECT_EQ(next_dist(0b0100, 1), 0b1000);
  EXPECT_EQ(next_dist(0b1000, 1), 0b10000);
}

TEST(Utils, next_dist_max_2) {
  EXPECT_EQ(next_dist(0b0000, 2), 0b0001);
  EXPECT_EQ(next_dist(0b0001, 2), 0b0010);
  EXPECT_EQ(next_dist(0b0010, 2), 0b0011);
  EXPECT_EQ(next_dist(0b0011, 2), 0b0100);
  EXPECT_EQ(next_dist(0b0100, 2), 0b0101);
  EXPECT_EQ(next_dist(0b0101, 2), 0b0110);
}

TEST(Utils, next_dist_with_exclusions) {
  auto result = next_dist(0b0000, 3, 0, 0b0100, 0);
  EXPECT_EQ(result.first, 0b0001);
  EXPECT_EQ(result.second, 1);

  result = next_dist(0b0001, 3, 1, 0b0010, 0);
  EXPECT_EQ(result.first, 0b0100);
  EXPECT_EQ(result.second, 4);

  result = next_dist(0b0000, 3, 0, 0, 0b0001);
  EXPECT_EQ(result.first, 0b0001);
  EXPECT_EQ(result.second, 1);

  result = next_dist(0b0001, 3, 1, 0, 0b0011);
  EXPECT_EQ(result.first, 0b0011);
  EXPECT_EQ(result.second, 3);
}

TEST(Utils, invert_dist_map_basic) {
  auto map = invert_dist_map(2, 2);
  EXPECT_EQ(map.size(), 4);
  EXPECT_EQ(map[0b00], 0);
  EXPECT_EQ(map[0b01], 1);
  EXPECT_EQ(map[0b10], 2);
  EXPECT_EQ(map[0b11], 3);
}

TEST(Utils, invert_dist_map_3_areas) {
  auto map = invert_dist_map(3, 2);
  EXPECT_EQ(map.size(), 8);
  EXPECT_EQ(map[0b000], 0);
  EXPECT_EQ(map[0b001], 1);
  EXPECT_EQ(map[0b010], 2);
  EXPECT_EQ(map[0b100], 4);
  EXPECT_EQ(map[0b011], 3);
  EXPECT_EQ(map[0b101], 5);
  EXPECT_EQ(map[0b110], 6);
  EXPECT_EQ(map[0b111], std::numeric_limits<size_t>::max());
}

TEST(Utils, lagrange_popcount) {
  EXPECT_EQ(lagrange_popcount(0b0), 0);
  EXPECT_EQ(lagrange_popcount(0b1), 1);
  EXPECT_EQ(lagrange_popcount(0b10), 1);
  EXPECT_EQ(lagrange_popcount(0b11), 2);
  EXPECT_EQ(lagrange_popcount(0b101), 2);
  EXPECT_EQ(lagrange_popcount(0b111), 3);
  EXPECT_EQ(lagrange_popcount(0b11111111), 8);
}

TEST(Utils, lagrange_bextr) {
  EXPECT_EQ(lagrange_bextr(0b1010, 0), 0);
  EXPECT_EQ(lagrange_bextr(0b1010, 1), 1);
  EXPECT_EQ(lagrange_bextr(0b1010, 2), 0);
  EXPECT_EQ(lagrange_bextr(0b1010, 3), 1);
}

TEST(Utils, check_excl_dist) {
  EXPECT_TRUE(check_excl_dist(0b0001, 0b0010));
  EXPECT_TRUE(check_excl_dist(0b0001, 0b1000));
  EXPECT_FALSE(check_excl_dist(0b0001, 0b0001));
  EXPECT_FALSE(check_excl_dist(0b0011, 0b0001));
}

TEST(Utils, check_incl_dist) {
  EXPECT_TRUE(check_incl_dist(0b0011, 0b0001));
  EXPECT_TRUE(check_incl_dist(0b0011, 0b0010));
  EXPECT_TRUE(check_incl_dist(0b0011, 0b0011));
  EXPECT_FALSE(check_incl_dist(0b0001, 0b0011));
  EXPECT_FALSE(check_incl_dist(0b0100, 0b0011));
}

TEST(Utils, convert_vector_to_lagrange_dist) {
  EXPECT_EQ(convert_vector_to_lagrange_dist({0, 0, 0}), 0b000);
  EXPECT_EQ(convert_vector_to_lagrange_dist({1, 0, 0}), 0b001);
  EXPECT_EQ(convert_vector_to_lagrange_dist({0, 1, 0}), 0b010);
  EXPECT_EQ(convert_vector_to_lagrange_dist({0, 0, 1}), 0b100);
  EXPECT_EQ(convert_vector_to_lagrange_dist({1, 1, 0}), 0b011);
  EXPECT_EQ(convert_vector_to_lagrange_dist({1, 1, 1}), 0b111);
}

TEST(Utils, lagrange_fast_log2) {
  EXPECT_EQ(lagrange_fast_log2(1), 1);
  EXPECT_EQ(lagrange_fast_log2(2), 2);
  EXPECT_EQ(lagrange_fast_log2(3), 2);
  EXPECT_EQ(lagrange_fast_log2(4), 3);
  EXPECT_EQ(lagrange_fast_log2(7), 3);
  EXPECT_EQ(lagrange_fast_log2(8), 4);
}
