#ifndef LAGRANGE_CONTEXT_H
#define LAGRANGE_CONTEXT_H

#include <memory>
#include <vector>

#include "Common.h"
#include "Operation.h"
#include "Workspace.h"
#include "tree.h"

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

  void registerTipClvs(
      const std::unordered_map<std::string, lagrange_dist_t>& dist_data);

  void optimize();
  double computeLH();
  double computeLLH();

  period_derivative_t computeDLLH(double initial_lh);

  std::string toString() const;

  void updateRates(const period_t& p);
  void init();

 private:
  void registerForwardOperations();
  void registerBackwardOperations();

  void computeForwardOperations();
  void computeBackwardOperations();

  std::shared_ptr<Tree> _tree;
  std::shared_ptr<Workspace> _workspace;
  std::vector<LHGoal> _lh_goal;
  std::vector<StateLHGoal> _state_lh_goal;
  std::vector<SplitLHGoal> _split_lh_goal;

  std::vector<SplitOperation> _forward_operations;
  std::vector<ReverseSplitOperation> _reverse_operations;
  std::shared_ptr<MakeRateMatrixOperation> _rate_matrix_op;
};

#endif /* end of include guard */
