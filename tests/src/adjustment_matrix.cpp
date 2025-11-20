#include "AdjustmentMatrix.hpp"
#include "gtest/gtest.h"

TEST(AdjustmentMatrix, simple_symmetric) {
  constexpr auto csv_string =
      "from, to, dist\na, b, 1.0\n a, c, 2.0\n b, c, 3.0";
  std::vector<std::string> areanames = {"c", "b", "a"};
  lagrange::AdjustmentMatrix adj(std::istringstream{csv_string}, areanames);
  EXPECT_EQ(adj.type(), lagrange::AdjustmentMatrixType::symmetric);
  auto matrix = adj.to_matrix();

  EXPECT_EQ(matrix[0 * 3 + 0], 0.0);
  EXPECT_EQ(matrix[1 * 3 + 1], 0.0);
  EXPECT_EQ(matrix[2 * 3 + 2], 0.0);

  EXPECT_EQ(matrix[0 * 3 + 1], 1.0);
  EXPECT_EQ(matrix[1 * 3 + 0], 1.0);

  EXPECT_EQ(matrix[0 * 3 + 2], 2.0);
  EXPECT_EQ(matrix[2 * 3 + 0], 2.0);

  EXPECT_EQ(matrix[1 * 3 + 2], 3.0);
  EXPECT_EQ(matrix[2 * 3 + 1], 3.0);
}

TEST(AdjustmentMatrix, simple_nonsymmetric) {
  constexpr auto csv_string =
      "from, to, dist\n"
      "a, b, 1.0\n"
      "c, a, 5.0\n"
      "a, c, 2.0\n"
      "b, c, 3.0\n"
      "b, a, 4.0\n"
      "c, b, 6.0";
  std::vector<std::string> areanames = {"c", "b", "a"};
  lagrange::AdjustmentMatrix adj(std::istringstream{csv_string}, areanames);
  EXPECT_EQ(adj.type(), lagrange::AdjustmentMatrixType::nonsymmetric);
  auto matrix = adj.to_matrix();

  EXPECT_EQ(matrix[0 * 3 + 0], 0.0);
  EXPECT_EQ(matrix[1 * 3 + 1], 0.0);
  EXPECT_EQ(matrix[2 * 3 + 2], 0.0);

  EXPECT_EQ(matrix[0 * 3 + 1], 1.0);
  EXPECT_EQ(matrix[1 * 3 + 0], 4.0);

  EXPECT_EQ(matrix[0 * 3 + 2], 2.0);
  EXPECT_EQ(matrix[2 * 3 + 0], 5.0);

  EXPECT_EQ(matrix[1 * 3 + 2], 3.0);
  EXPECT_EQ(matrix[2 * 3 + 1], 6.0);
}

TEST(AdjustmentMatrix, simple_check_error1) {
  constexpr auto csv_string =
      "from, to, dist\n"
      "a, b, 1.0\n"
      "c, a, 5.0\n"
      "a, c, 2.0\n"
      "b, c, 3.0\n"
      "b, a, 4.0\n"
      "b, c, 6.0";
  std::vector<std::string> areanames = {"a", "b", "c"};
  EXPECT_ANY_THROW(lagrange::AdjustmentMatrix adj(
      std::istringstream{csv_string}, areanames));
}

TEST(AdjustmentMatrix, simple_check_error2) {
  constexpr auto csv_string =
      "from, to, dist\n"
      "a, b, 1.0\n"
      "c, a, 5.0\n"
      "a, c, 2.0\n"
      "b, c, 3.0\n"
      "b, a, 4.0";
  std::vector<std::string> areanames = {"a", "b", "c"};
  EXPECT_ANY_THROW(lagrange::AdjustmentMatrix adj(
      std::istringstream{csv_string}, areanames));
}

TEST(AdjustmentMatrix, simple_check_error3) {
  constexpr auto csv_string =
      "from, to, dist\na, c, 2.0\n b, c, 3.0";
  std::vector<std::string> areanames = {"a", "b", "c"};
  EXPECT_ANY_THROW(lagrange::AdjustmentMatrix adj(
      std::istringstream{csv_string}, areanames));
}
