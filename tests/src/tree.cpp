#include "Tree.hpp"

#include <memory>
#include <string>
#include <unordered_map>

#include "Common.hpp"
#include "Periods.hpp"
#include "TreeReader.hpp"
#include "Workspace.hpp"
#include "gtest/gtest.h"

using namespace lagrange;

class TreeTest : public ::testing::Test {
 protected:
  void SetUp() override {
    _basic_tree_newick = "((a:1.0,b:1.0):1.0, c:2.0);";
    _basic_tree_dist_data = {{"a", 0b01}, {"b", 0b11}, {"c", 0b01}};
    _basic_ws = std::make_shared<Workspace>(/*taxa=*/3,
                                            /*regions=*/2,
                                            /*maxareas=*/2);
    /*
    _arbitrary_rate_matrix = {
        {0.00000, 0.00000, 0.00000, 0.00000},
        {0.08851, -1.71766, 0.00000, 1.62915},
        {0.75966, 0.00000, -0.81937, 0.05971},
        {0.00000, 0.37406, 1.26782, -1.64188},
    };
    */

    _arbitrary_rate_matrix.reset(new LagrangeMatrixBase[4 * 4]{
        0.00000,
        0.00000,
        0.00000,
        0.00000,
        0.08851,
        -1.71766,
        0.00000,
        1.62915,
        0.75966,
        0.00000,
        -0.81937,
        0.05971,
        0.00000,
        0.37406,
        1.26782,
        -1.64188,
    });

    _rate_matrix_op = std::make_shared<MakeRateMatrixOperation>(
        _rate_matrix_index, _basic_period_index);
  }

  std::shared_ptr<Tree> parse_tree(const std::string& newick) {
    return TreeReader().readTree(newick);
  }

  size_t _rate_matrix_index = 0;

  std::string _basic_tree_newick;

  std::unordered_map<std::string, Range> _basic_tree_dist_data;

  std::shared_ptr<Workspace> _basic_ws;

  std::shared_ptr<MakeRateMatrixOperation> _rate_matrix_op;

  std::unique_ptr<LagrangeMatrixBase[]> _arbitrary_rate_matrix;

  PeriodTimes _basic_periods;
  const size_t _basic_period_index = 0;
};

TEST_F(TreeTest, simple0) { parse_tree(_basic_tree_newick); }

TEST_F(TreeTest, generate) {
  auto t = parse_tree(_basic_tree_newick);
  t->setPeriods(_basic_periods);
  PeriodRateMatrixMap rm_map;
  BranchProbMatrixMap pm_map;
  auto ops = t->generateForwardOperations(*_basic_ws, rm_map, pm_map);

  EXPECT_EQ(ops.size(), 2);
}

TEST_F(TreeTest, generateOperationsSimple1) {
  auto t = parse_tree(_basic_tree_newick);
  t->setPeriods(_basic_periods);
  PeriodRateMatrixMap rm_map;
  BranchProbMatrixMap pm_map;
  auto ops = t->generateBackwardOperations(*_basic_ws, rm_map, pm_map);

  EXPECT_EQ(ops.size(), 4);
}
