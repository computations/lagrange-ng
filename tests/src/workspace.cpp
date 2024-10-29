#include <stdexcept>

#include "Workspace.hpp"
#include "gtest/gtest.h"

using namespace lagrange;

TEST(Workspace, simple0) {
  constexpr size_t regions = 3;
  constexpr size_t states = 1 << regions;
  Workspace ws(10, regions, regions);

  EXPECT_EQ(ws.suggestProbMatrixIndex(), 0);
  EXPECT_EQ(ws.reserveRateMatrixIndex(0), 0);

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

TEST(Workspace, maxareas) { EXPECT_ANY_THROW(Workspace(2, 2, 12)); }

/*
TEST(Workspace, setters) {
  constexpr size_t regions = 3;
  constexpr size_t states = 1 << regions;
  Workspace ws(10, regions);
  size_t clv_index = ws.register_generic_clv();
  ws.reserve();
  ws.get_base_frequencies(0) = 1.0;
  EXPECT_EQ(ws.get_base_frequencies(0)[0], 1.0);
}
*/
