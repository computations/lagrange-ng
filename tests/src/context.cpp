#include "Context.hpp"

#include <memory>
#include <unordered_map>

#include "Common.hpp"
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
    _basic_periods = PeriodTimes();
    _basic_tree->setPeriods(_basic_periods);
    _basic_tree_data = {{"a", 0b01}, {"b", 0b01}, {"c", 0b11}};
  }

  std::shared_ptr<Tree> parse_tree(const std::string& newick) {
    return TreeReader().readTree(newick);
  }

  std::string _basic_tree_newick;
  std::shared_ptr<Tree> _basic_tree;
  std::unordered_map<std::string, Range> _basic_tree_data;
  PeriodTimes _basic_periods;
  size_t _basic_tree_data_region_count = 2;
  // WorkerState _worker_state;
};

TEST_F(ContextTest, simple0) {
  Context context(_basic_tree, 2, 2);
  context.registerLHGoal();
  EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
              == SetCLVStatus::definite);
  context.init({{.dispersion_rate = 0.1,
                 .extinction_rate = 0.1,
                 .distance_penalty = 1.0,
                 .name = "",
                 .regions = 2}});
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
  context.init({{.dispersion_rate = 0.01,
                 .extinction_rate = 0.01,
                 .distance_penalty = 1.0,
                 .name = "",
                 .regions = 2}});
  EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
              == SetCLVStatus::definite);

  auto worker_context = context.makeThreadContext();
  worker_context.setTotalThreads(1);
  worker_context.initBarrier();

  WorkerState worker_state(0);
  double llh = context.computeLLH(worker_state, worker_context);
  constexpr double regression_llh = -1.9960966944483829;

  EXPECT_NEAR(llh, regression_llh, 1e-9);

  llh = context.computeLLH(worker_state, worker_context);
  EXPECT_NEAR(llh, regression_llh, 1e-9);
}

TEST_F(ContextTest, optimizeSimple0) {
  Context context(_basic_tree, 2, 2);
  context.registerLHGoal();
  context.init({{.dispersion_rate = 10.5,
                 .extinction_rate = 1.5,
                 .distance_penalty = 1,
                 .name = "",
                 .regions = 2}});
  // context.updateRates({{10.5, 1.5}});
  context.set_opt_method(OptimizationMethod::BFGS);

  auto worker_context = context.makeThreadContext();
  worker_context.setTotalThreads(1);
  worker_context.initBarrier();

  WorkerState worker_state(0);

  EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
              == SetCLVStatus::definite);

  double initial_llh = context.computeLLH(worker_state, worker_context);
  context.setRunMode(LagrangeOperationMode::OPTIMIZE);
  context.optimizeAndComputeValues(worker_state, worker_context, false, false);
  double llh = context.computeLLH(worker_state, worker_context);

  EXPECT_GT(llh, initial_llh);
}

TEST_F(ContextTest, evaluateUsesConfiguredRates) {
  auto evaluate = [this](double dispersion, double extinction) {
    auto tree = parse_tree(_basic_tree_newick);
    tree->setPeriods(_basic_periods);
    Context context(tree, 2, 2);
    context.registerLHGoal();
    context.init({{.dispersion_rate = dispersion,
                   .extinction_rate = extinction,
                   .distance_penalty = 1.0,
                   .name = "",
                   .regions = 2}});

    auto worker_context = context.makeThreadContext();
    worker_context.setTotalThreads(1);
    worker_context.initBarrier();

    WorkerState worker_state(0);
    EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
                == SetCLVStatus::definite);

    double llh = context.computeLLH(worker_state, worker_context);
    context.setRunMode(LagrangeOperationMode::EVALUATE);
    context.optimizeAndComputeValues(
        worker_state, worker_context, false, false);

    auto params = context.currentParams();
    EXPECT_DOUBLE_EQ(params[0].dispersion_rate, dispersion);
    EXPECT_DOUBLE_EQ(params[0].extinction_rate, extinction);
    return llh;
  };

  double first_llh = evaluate(0.1, 0.2);
  double second_llh = evaluate(2.0, 3.0);
  EXPECT_NE(first_llh, second_llh);
}

TEST_F(ContextTest, StateGoal0) {
  Context context(_basic_tree, 2, 2);
  context.registerLHGoal();
  context.registerStateLHGoal();
  context.init({{.dispersion_rate = 10.5,
                 .extinction_rate = 1.5,
                 .distance_penalty = 1,
                 .name = "",
                 .regions = 2}});
  // context.updateRates({{10.5, 1.5}});

  auto worker_context = context.makeThreadContext();
  worker_context.setTotalThreads(1);
  worker_context.initBarrier();

  WorkerState worker_state(0);

  EXPECT_TRUE(context.registerTipClvs(_basic_tree_data)
              == SetCLVStatus::definite);

  context.computeLLH(worker_state, worker_context);
  auto states = context.computeAndGetStateGoals(worker_state, worker_context);
  EXPECT_EQ(states.size(), 2);
}
