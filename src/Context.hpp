/*
 * Context.h
 *     Created: 2021-03-23
 *      Author: Ben Bettisworth
 */

#ifndef LAGRANGE_CONTEXT_H
#define LAGRANGE_CONTEXT_H

#include <memory>
#include <nlopt.hpp>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "AncSplit.hpp"
#include "Common.hpp"
#include "ConfigFile.hpp"
#include "Operation.hpp"
#include "Tree.hpp"
#include "WorkerState.hpp"
#include "Workspace.hpp"

namespace lagrange {

class Context {
 public:
  Context(std::shared_ptr<Tree> tree, size_t regions, size_t max_areas) :
      _tree{std::move(tree)},
      _workspace{std::make_shared<Workspace>(
          _tree->getTipCount(), regions, max_areas)} {}

  void registerLHGoal();
  void registerStateLHGoal();
  void registerStateLHGoal(
      const std::unordered_set<std::string>& mrca_keys,
      const std::unordered_map<std::string, std::shared_ptr<MRCAEntry>>&
          mrca_map);
  void registerSplitLHGoal();
  void registerSplitLHGoal(
      const std::unordered_set<std::string>& mrca_keys,
      const std::unordered_map<std::string, std::shared_ptr<MRCAEntry>>&
          mrca_map);

  [[nodiscard]] auto getStateResults() const
      -> std::unordered_map<size_t, std::unique_ptr<LagrangeMatrixBase[]>>;
  [[nodiscard]] auto getSplitResults() const -> SplitReturnList;

  [[nodiscard]] SetCLVStatus registerTipClvs(
      const std::unordered_map<std::string, Range>& dist_data);

  void optimizeAndComputeValues(WorkerState& ts,
                                bool states,
                                bool splits,
                                const LagrangeOperationMode& mode);

  auto computeLLH(WorkerState& ts) -> double;
  auto computeLLH(WorkerState& ts, WorkerContext& tc) -> double;
  auto computeStateGoal(WorkerState& ts)
      -> std::unordered_map<size_t, std::unique_ptr<LagrangeMatrixBase[]>>;
  auto computeSplitGoal(WorkerState& ts) -> SplitReturnList;

  [[nodiscard]] auto toString() const -> std::string;

  void updateRates(const std::vector<PeriodParams>& p);
  void init();

  [[nodiscard]] auto stateCount() const -> size_t {
    return _workspace->restrictedStateCount();
  }

  void setParams(double e, double d) { _workspace->setPeriodParams(0, e, d); }

  [[nodiscard]] auto currentParams() const -> std::vector<PeriodParams>;

  auto makeThreadContext() -> WorkerContext {
    WorkerContext tc{_forward_operations,
                     _reverse_operations,
                     _llh_goal,
                     _state_lh_goal,
                     _split_lh_goal};
    return tc;
  }

  void set_lh_epsilon(double lhe) { _lh_epsilon = lhe; }

  void useArnoldi(bool mode_set = true, bool adaptive = true) const;

  [[nodiscard]] auto getPeriodCount() const -> size_t;

  void set_opt_method(const OptimizationMethod&);

 private:
  void registerForwardOperations();
  void registerBackwardOperations();

  void computeForwardOperations(WorkerState& ts, WorkerContext& tc);
  void computeBackwardOperations(WorkerState& ts, WorkerContext& tc);

  void computeStateGoal(WorkerState& ts, WorkerContext& tc);
  void computeSplitGoal(WorkerState& ts, WorkerContext& tc);

  auto optimize(WorkerState& ts, WorkerContext& tc) -> double;

  void extractRateMatrixOperations();

  void registerGoals(const std::function<void(Node&)>& func);
  void registerGoals(
      const std::unordered_set<std::string>& mrca_keys,
      const std::unordered_map<std::string, std::shared_ptr<MRCAEntry>>&
          mrca_map,
      const std::function<void(Node&)>& func);

  [[nodiscard]] auto computeDerivative() const -> bool;

  auto makeStateGoalCB() -> std::function<void(Node&)>;
  auto makeSplitGoalCB() -> std::function<void(Node&)>;

  nlopt::algorithm _opt_method;

  double _lh_epsilon;

  std::shared_ptr<Tree> _tree;
  std::shared_ptr<Workspace> _workspace;
  std::vector<LLHGoal> _llh_goal;
  std::vector<StateLHGoal> _state_lh_goal;
  std::vector<SplitLHGoal> _split_lh_goal;

  std::vector<std::shared_ptr<SplitOperation>> _forward_operations;
  std::vector<std::shared_ptr<ReverseSplitOperation>> _reverse_operations;
  std::vector<std::shared_ptr<MakeRateMatrixOperation>> _rate_matrix_ops;
};

}  // namespace lagrange
#endif /* end of include guard */
