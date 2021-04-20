#include <stdexcept>

#include "Workspace.h"
#include "environment.hpp"
#include "gtest/gtest.h"

TEST(Workspace, simple0) {
  constexpr size_t regions = 3;
  constexpr size_t states = 1 << regions;
  Workspace ws(10, regions);

  ws.reserve();

  EXPECT_EQ(ws.states(), states);

  EXPECT_EQ(ws.rate_matrix_count(), 1);
  EXPECT_EQ(ws.rate_matrix(0).rows(), states);
  EXPECT_EQ(ws.rate_matrix(0).columns(), states);
  EXPECT_THROW(ws.rate_matrix(1), std::runtime_error);

  EXPECT_EQ(ws.rate_matrix_count(), 1);
  EXPECT_EQ(ws.prob_matrix(0).rows(), states);
  EXPECT_EQ(ws.prob_matrix(0).columns(), states);
  EXPECT_THROW(ws.prob_matrix(1), std::runtime_error);
}

TEST(Workspace, simple1) {
  constexpr size_t regions = 3;
  constexpr size_t states = 1 << regions;
  Workspace ws(10, regions);
  size_t clv_index = ws.register_generic_clv();
  ws.reserve();
  EXPECT_EQ(clv_index, 0);
  EXPECT_EQ(ws.clv(clv_index).size(), states);
}

TEST(Workspace, minsize) {
  EXPECT_ANY_THROW(Workspace(0, 2));
  EXPECT_NO_THROW(Workspace(1, 2));
}

TEST(Workspace, setters) {
  constexpr size_t regions = 3;
  constexpr size_t states = 1 << regions;
  Workspace ws(10, regions);
  size_t clv_index = ws.register_generic_clv();
  ws.reserve();
  ws.get_base_frequencies(0) = 1.0;
  EXPECT_EQ(ws.get_base_frequencies(0)[0], 1.0);
}
