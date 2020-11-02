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
  Context(std::shared_ptr<Tree>, size_t regions);

  void registerLHGoal();
  void registerStateLHGoal();
  void registerSplitLHGoal();

 private:
  void registerForwardOperations();
  void registerBackwardOperations();

  Workspace _workspace;
  std::vector<LHGoal> _lh_goal;
  std::vector<StateLHGoal> _state_lh_goal;
  std::vector<SplitLHGoal> _split_lh_goal;

  std::vector<SplitOperation> _forward_operations;
  std::vector<ReverseSplitOperation> _reverse_operations;

  std::shared_ptr<Tree> _tree;
};

#endif /* end of include guard */
