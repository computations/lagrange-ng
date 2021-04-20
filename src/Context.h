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
#include <vector>

#include "AncSplit.h"
#include "Common.h"
#include "Operation.h"
#include "ThreadState.h"
#include "Workspace.h"
#include "tree.h"

typedef std::vector<std::unordered_map<lagrange_dist_t, std::vector<AncSplit>>>
    lagrange_split_goal_return_t;

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

  double computeLHGoal();
  std::vector<lagrange_col_vector_t> computeStateGoal();
  lagrange_split_goal_return_t computeSplitGoal();

  void registerTipClvs(
      const std::unordered_map<std::string, lagrange_dist_t>& dist_data);

  double optimize();
  double computeLLH();

  period_derivative_t computeDLLH(double initial_lh);

  std::string toString() const;

  void updateRates(const period_t& p);
  void init();

  period_t currentParams() const;

  std::string treeCLVStatus() const;

 private:
  void registerForwardOperations();
  void registerBackwardOperations();

  void computeForwardOperations();
  void computeBackwardOperations();

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
