/*
 * tree.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef TREE_H
#define TREE_H

#include <cstddef>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include "Alignment.hpp"
#include "AncSplit.hpp"
#include "Common.hpp"
#include "Fossil.hpp"
#include "MRCA.hpp"
#include "Node.hpp"
#include "Operation.hpp"
#include "Workspace.hpp"

namespace lagrange {

class Tree {
 public:
  Tree();
  explicit Tree(std::shared_ptr<Node> root);

  ~Tree();

  auto generateForwardOperations(Workspace &ws)
      -> std::vector<std::shared_ptr<SplitOperation>>;

  auto generateBackwardOperations(Workspace &ws)
      -> std::vector<std::shared_ptr<ReverseSplitOperation>>;

  auto generateForwardOperations(Workspace &ws,
                                 PeriodRateMatrixMap &period_rm_map,
                                 BranchProbMatrixMap &period_pm_map)
      -> std::vector<std::shared_ptr<SplitOperation>>;
  auto generateBackwardOperations(Workspace &ws,
                                  PeriodRateMatrixMap &period_rm_map,
                                  BranchProbMatrixMap &period_pm_map)
      -> std::vector<std::shared_ptr<ReverseSplitOperation>>;

  auto traversePreorderInternalNodesOnly() const -> std::vector<size_t>;
  auto traversePreorderInternalNodesOnlyNumbers() const -> std::vector<size_t>;

  auto inverseNodeIdMap() const -> std::unordered_map<size_t, size_t>;

  void applyPreorderInternalOnly(const std::function<void(Node &)> &func);

  [[nodiscard]] SetCLVStatus assignTipData(
      Workspace &ws, const std::unordered_map<std::string, Range> &dist_data);

  void assignMCRALabels(const MRCAMap &mrca_map);

  auto getNode(size_t num) -> std::shared_ptr<Node>;

  auto getNodeCount() const -> size_t;
  auto getTipCount() const -> size_t;
  auto getInternalCount() const -> size_t;

  auto getRoot() -> std::shared_ptr<Node>;
  void processRoot();

  auto getMRCA(const std::shared_ptr<MRCAEntry> &mrca) -> std::shared_ptr<Node>;

  auto getParent(const std::shared_ptr<Node> &n) const -> std::shared_ptr<Node>;

  void setHeightTopDown();
  void setHeightBottomUp();
  void setPeriods(const Periods &periods);

  auto getNewick() const -> std::string;
  auto getNewickLambda(
      const std::function<std::string(const Node &)> &newick_lambda) const
      -> std::string;

  auto checkAlignmentConsistency(const Alignment &align) const -> bool;

  void assignFossils(const std::vector<Fossil> &fossils);

  void assignStateResult(std::unique_ptr<LagrangeMatrixBase[]>,
                         const MRCALabel &);
  void assignSplitResult(const SplitReturn &, const MRCALabel &);

 private:
  void processReRoot(const std::shared_ptr<Node> &node);
  static void exchangeInfo(const std::shared_ptr<Node> &node1,
                           const std::shared_ptr<Node> &node2);
  void postOrderProcessRoot(const std::shared_ptr<Node> &node);
  auto getGreatestDistance(std::shared_ptr<Node> inNode) -> double;

  auto findNode(const std::shared_ptr<Node> &n) -> bool;

  std::shared_ptr<Node> _root;
  MRCAMap _mrcas;
  size_t _node_count;
  size_t _tip_count;
};

}  // namespace lagrange

#endif /* TREE_H_ */
