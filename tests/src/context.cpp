#include "Context.hpp"

#include <memory>
#include <unordered_map>

#include "Common.hpp"
#include "ConfigFile.hpp"
#include "Periods.hpp"
#include "TreeReader.hpp"
#include "WorkerState.hpp"
#include "Workspace.hpp"
#include "gtest/gtest.h"

using namespace lagrange;

class ContextTest : public ::testing::Test {
 protected:
  void SetUp() override {
    _basic_tree_newick = "((a:1.0,b:1.0)n:1.0, c:2.0)r;";
    _basic_tree = parse_tree(_basic_tree_newick);
    _basic_periods = Periods();
    _basic_tree->setPeriods(_basic_periods);
    _basic_tree_data = {{"a", 0b01}, {"b", 0b01}, {"c", 0b11}};
  }

  std::shared_ptr<Tree> parse_tree(const std::string& newick) {
    return TreeReader().readTree(newick);
  }

  std::string _basic_tree_newick;
  std::shared_ptr<Tree> _basic_tree;
  std::unordered_map<std::string, Range> _basic_tree_data;
  Periods _basic_periods;
  size_t _basic_tree_data_region_count = 2;
  WorkerState _worker_state;
};

TEST_F(ContextTest, simple0) {
  Context context(_basic_tree, 2, 2);
  context.registerLHGoal();
  EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
              == SetCLVStatus::definite);
  context.init();
}

TEST_F(ContextTest, error0) {
  Context context(_basic_tree, 2, 2);
  EXPECT_THROW(EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
                           == SetCLVStatus::definite),
               std::runtime_error);
}

TEST_F(ContextTest, computelh1) {
  Context context(_basic_tree, 2, 2);
  context.registerLHGoal();
  context.init();
  EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
              == SetCLVStatus::definite);

  double llh = context.computeLLH(_worker_state);
  constexpr double regression_llh = -1.9960966944483829;

  EXPECT_NEAR(llh, regression_llh, 1e-9);

  llh = context.computeLLH(_worker_state);
  EXPECT_NEAR(llh, regression_llh, 1e-9);
}

TEST_F(ContextTest, optimizeSimple0) {
  Context context(_basic_tree, 2, 2);
  context.registerLHGoal();
  context.init();
  context.updateRates({{10.5, 1.5}});
  context.set_opt_method(OptimizationMethod::BFGS);
  EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
              == SetCLVStatus::definite);

  auto wc = context.makeThreadContext();
  double initial_llh = context.computeLLH(_worker_state);
  context.optimizeAndComputeValues(_worker_state,
                                   wc,
                                   false,
                                   false,
                                   LagrangeOperationMode::OPTIMIZE);
  double llh = context.computeLLH(_worker_state);

  EXPECT_GT(llh, initial_llh);
}

TEST_F(ContextTest, StateGoal0) {
  Context context(_basic_tree, 2, 2);
  context.registerLHGoal();
  context.registerStateLHGoal();
  context.init();
  context.updateRates({{10.5, 1.5}});
  EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
              == SetCLVStatus::definite);

  context.computeLLH(_worker_state);
  auto states = context.computeStateGoal(_worker_state);
  EXPECT_EQ(states.size(), 2);
}
