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
#include <utility>
#include <vector>

#include "Alignment.h"
#include "Common.h"
#include "Fossil.h"
#include "Node.h"
#include "Operation.h"

class Tree {
 private:
  std::shared_ptr<Node> _root;
  std::vector<std::shared_ptr<Node>> _nodes;
  std::vector<std::shared_ptr<Node>> _internal_nodes;
  std::vector<std::shared_ptr<Node>> _external_nodes;
  std::unordered_map<std::string, std::shared_ptr<Node>> _mrcas;
  size_t _internal_node_count;
  size_t _external_node_count;

  void processReRoot(const std::shared_ptr<Node> &node);
  static void exchangeInfo(const std::shared_ptr<Node> &node1,
                           const std::shared_ptr<Node> &node2);
  void postOrderProcessRoot(const std::shared_ptr<Node> &node);
  auto getGreatestDistance(std::shared_ptr<Node> inNode) -> double;

  auto findNode(const std::shared_ptr<Node> &n) -> bool;

 public:
  Tree();
  explicit Tree(std::shared_ptr<Node> root);

  void addExternalNode(const std::shared_ptr<Node> &tn);
  void addInternalNode(const std::shared_ptr<Node> &tn);
  void pruneExternalNode(std::shared_ptr<Node> node);

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

  void applyPreorderInternalOnly(const std::function<void(Node &)> &func);

  void assignTipData(
      Workspace &ws,
      const std::unordered_map<std::string, lagrange_dist_t> &dist_data);

  auto getExternalNode(size_t num) -> std::shared_ptr<Node>;
  auto getExternalNode(const std::string &name) -> std::shared_ptr<Node>;
  auto getInternalNode(size_t num) -> std::shared_ptr<Node>;
  auto getInternalNode(const std::string &name) -> std::shared_ptr<Node>;
  auto getNode(size_t num) -> std::shared_ptr<Node>;

  auto getNodeCount() const -> unsigned int;
  auto getExternalNodeCount() const -> unsigned int;
  auto getInternalNodeCount() const -> unsigned int;

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

  bool checkAlignmentConsistency(const Alignment &align) const;

  void assignFossils(const std::vector<Fossil> &fossils);

  ~Tree();
};

#endif /* TREE_H_ */
