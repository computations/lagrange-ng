/*
 * Context.h
 *     Created: 2021-03-23
 *      Author: Ben Bettisworth
 */

#ifndef LAGRANGE_CONTEXT_H
#define LAGRANGE_CONTEXT_H

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "AncSplit.hpp"
#include "Common.hpp"
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
          _tree->getTipCount(), regions, max_areas)},
      _rate_matrix_ops{} {}

  void registerLHGoal();
  void registerStateLHGoal();
  void registerStateLHGoal(
      const std::vector<std::string>& mrca_keys,
      const std::unordered_map<std::string, std::shared_ptr<MRCAEntry>>
          mrca_map);
  void registerSplitLHGoal();
  void registerSplitLHGoal(
      const std::vector<std::string>& mrca_keys,
      const std::unordered_map<std::string, std::shared_ptr<MRCAEntry>>
          mrca_map);

  auto getStateResults() const
      -> std::vector<std::unique_ptr<LagrangeMatrixBase[]>>;
  auto getSplitResults() const -> SplitReturnList;

  void registerTipClvs(const std::unordered_map<std::string, Range>& dist_data);

  void optimizeAndComputeValues(WorkerState& ts,
                                bool states,
                                bool splits,
                                const LagrangeOperationMode& mode);

  auto computeLLH(WorkerState& ts) -> double;
  auto computeLLH(WorkerState& ts, WorkerContext& tc) -> double;
  auto computeStateGoal(WorkerState& ts)
      -> std::vector<std::unique_ptr<LagrangeMatrixBase[]>>;
  auto computeSplitGoal(WorkerState& ts) -> SplitReturnList;

  auto toString() const -> std::string;

  void updateRates(const std::vector<PeriodParams>& p);
  void init();

  auto stateCount() const -> size_t {
    return _workspace->restrictedStateCount();
  }

  void setParams(double e, double d) { _workspace->setPeriodParams(0, e, d); }

  auto currentParams() const -> std::vector<PeriodParams>;

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

  size_t getPeriodCount() const;

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
      const std::vector<std::string> mrca_keys,
      const std::unordered_map<std::string, std::shared_ptr<MRCAEntry>>&
          mrca_map,
      const std::function<void(Node&)>& func);

  std::function<void (Node &)> makeStateGoalCB();
  std::function<void (Node &)> makeSplitGoalCB();

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
