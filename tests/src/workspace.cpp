#include "Workspace.h"
#include "environment.hpp"
#include "gtest/gtest.h"
#include <stdexcept>

TEST(Workspace, simple0) {
  constexpr size_t regions = 3;
  constexpr size_t states = 1 << regions;
  auto ws = Workspace(10, regions);

  EXPECT_EQ(ws.states(), states);

  EXPECT_EQ(ws.clv_count(), 3 * 9 + 10);
  EXPECT_EQ(ws.clv(0).size(), states);
  EXPECT_THROW(ws.clv(ws.clv_count()), std::runtime_error);

  EXPECT_EQ(ws.rate_matrix_count(), 1);
  EXPECT_EQ(ws.rate_matrix(0).rows(), states);
  EXPECT_EQ(ws.rate_matrix(0).columns(), states);
  EXPECT_THROW(ws.rate_matrix(1), std::runtime_error);

  EXPECT_EQ(ws.rate_matrix_count(), 1);
  EXPECT_EQ(ws.prob_matrix(0).rows(), states);
  EXPECT_EQ(ws.prob_matrix(0).columns(), states);
  EXPECT_THROW(ws.prob_matrix(1), std::runtime_error);
}

TEST(Workspace, minsize) {
  EXPECT_THROW(Workspace(0, 2), std::bad_alloc);
  EXPECT_NO_THROW(Workspace(1, 2));
}
