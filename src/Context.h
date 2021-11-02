/*
 * Context.h
 *     Created: 2021-03-23
 *      Author: Ben Bettisworth
 */

#ifndef LAGRANGE_CONTEXT_H
#define LAGRANGE_CONTEXT_H

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "AncSplit.h"
#include "Common.h"
#include "Operation.h"
#include "Tree.h"
#include "WorkerState.h"
#include "Workspace.h"

class Context {
 public:
  Context(std::shared_ptr<Tree> tree, size_t regions)
      : _tree{tree},
        _workspace{std::make_shared<Workspace>(_tree->getExternalNodeCount(),
                                               regions)},
        _rate_matrix_op{std::make_shared<MakeRateMatrixOperation>(
            _workspace->suggest_rate_matrix_index())} {}

  void registerLHGoal();
  void registerStateLHGoal();
  void registerSplitLHGoal();

  std::vector<std::unique_ptr<lagrange_matrix_base_t[]>> getStateResults();
  lagrange_split_list_t getSplitResults();

  void registerTipClvs(
      const std::unordered_map<std::string, lagrange_dist_t>& dist_data);

  void optimizeAndComputeValues(WorkerState& ts, bool states, bool splits,
                                bool output);
  double computeLLH(WorkerState& ts);
  double computeLLH(WorkerState& ts, WorkerContext& tc);
  std::vector<std::unique_ptr<lagrange_matrix_base_t[]>> computeStateGoal(
      WorkerState& ts);
  lagrange_split_list_t computeSplitGoal(WorkerState& ts);

  std::string toString() const;

  void updateRates(const period_t& p);
  void init();

  size_t stateCount() const { return _workspace->states(); }

  period_t currentParams() const;

  WorkerContext makeThreadContext() {
    WorkerContext tc{_forward_operations, _reverse_operations, _llh_goal,
                     _state_lh_goal, _split_lh_goal};
    return tc;
  }

 private:
  void registerForwardOperations();
  void registerBackwardOperations();

  void computeForwardOperations(WorkerState& ts, WorkerContext& tc);
  void computeBackwardOperations(WorkerState& ts, WorkerContext& tc);

  void computeStateGoal(WorkerState& ts, WorkerContext& tc);
  void computeSplitGoal(WorkerState& ts, WorkerContext& tc);

  double optimize(WorkerState& ts, WorkerContext& tc);

  std::shared_ptr<Tree> _tree;
  std::shared_ptr<Workspace> _workspace;
  std::vector<LLHGoal> _llh_goal;
  std::vector<StateLHGoal> _state_lh_goal;
  std::vector<SplitLHGoal> _split_lh_goal;

  std::vector<std::shared_ptr<SplitOperation>> _forward_operations;
  std::vector<std::shared_ptr<ReverseSplitOperation>> _reverse_operations;
  std::shared_ptr<MakeRateMatrixOperation> _rate_matrix_op;
};

#endif /* end of include guard */
