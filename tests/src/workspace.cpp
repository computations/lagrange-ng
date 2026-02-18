#include <stdexcept>

#include "Workspace.hpp"
#include "gtest/gtest.h"

using namespace lagrange;

TEST(Workspace, simple0) {
  constexpr size_t regions = 3;
  constexpr size_t states = 1 << regions;
  Workspace ws(10, regions, regions);

  ws.reserve();

  EXPECT_EQ(ws.states(), states);

  EXPECT_EQ(ws.rateMatrixCount(), 1);
  EXPECT_THROW(ws.rateMatrix(1), std::runtime_error);

  EXPECT_EQ(ws.probMatrixCount(), 1);
  EXPECT_THROW(ws.probMatrix(1), std::runtime_error);
}

TEST(Workspace, simple1) {
  constexpr size_t regions = 3;
  Workspace ws(10, regions, regions);
  size_t clv_index = ws.registerGenericCLV();
  ws.reserve();
  EXPECT_EQ(clv_index, 0);
}

TEST(Workspace, minsize) {
  EXPECT_ANY_THROW(Workspace(0, 2, 2));
  EXPECT_NO_THROW(Workspace(1, 2, 2));
}

TEST(Workspace, maxareas_throws) { EXPECT_ANY_THROW(Workspace(2, 2, 12)); }
