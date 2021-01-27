#include <iomanip>
#include <memory>
#include <unordered_map>

#include "Context.h"
#include "environment.hpp"
#include "gtest/gtest.h"
#include "tree_reader.h"

class ContextTest : public ::testing::Test {
 protected:
  void SetUp() override {
    _basic_tree_newick = "((a:1.0,b:1.0):1.0, c:2.0);";
    _basic_tree = parse_tree(_basic_tree_newick);
    _basic_tree_data = {{"a", 0b01}, {"b", 0b01}, {"c", 0b11}};
  }

  std::shared_ptr<Tree> parse_tree(const std::string& newick) {
    return TreeReader().readTree(newick);
  }

  std::string _basic_tree_newick;
  std::shared_ptr<Tree> _basic_tree;
  std::unordered_map<std::string, lagrange_dist_t> _basic_tree_data;
  size_t _basic_tree_data_region_count = 2;
};

TEST_F(ContextTest, simple0) {
  Context context(_basic_tree, 2);
  context.registerLHGoal();
  context.registerTipClvs(_basic_tree_data);
  context.init();
}

TEST_F(ContextTest, error0) {
  Context context(_basic_tree, 2);
  EXPECT_THROW(context.registerTipClvs(_basic_tree_data), std::runtime_error);
}

TEST_F(ContextTest, computelh1) {
  Context context(_basic_tree, 2);
  context.registerLHGoal();
  context.init();
  context.registerTipClvs(_basic_tree_data);

  double llh = context.computeLLH();
  constexpr double regression_llh = -1.7596288538749982;

  EXPECT_NEAR(llh, regression_llh, 1e-9);
}

TEST_F(ContextTest, optimizeSimple0) {
  Context context(_basic_tree, 2);
  context.registerLHGoal();
  context.init();
  context.updateRates({10.5, 1.5});
  context.registerTipClvs(_basic_tree_data);

  double initial_llh = context.computeLLH();
  context.optimize();
  double llh = context.computeLLH();

  EXPECT_GT(llh, initial_llh);
}

TEST_F(ContextTest, LHGoal0) {
  Context context(_basic_tree, 2);
  context.registerLHGoal();
  context.init();
  context.updateRates({10.5, 1.5});
  context.registerTipClvs(_basic_tree_data);

  context.computeLHGoal();
}

TEST_F(ContextTest, StateGoal0) {
  Context context(_basic_tree, 2);
  context.registerLHGoal();
  context.registerStateLHGoal();
  context.init();
  context.updateRates({10.5, 1.5});
  context.registerTipClvs(_basic_tree_data);

  context.computeLHGoal();
  auto states = context.computeStateGoal();
  EXPECT_EQ(states.size(), 2);
}
