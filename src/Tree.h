/*
 * tree.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef TREE_H_
#define TREE_H_

#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Common.h"
#include "Node.h"
#include "Operation.h"

class Tree {
 private:
  std::shared_ptr<Node> _root;
  std::vector<std::shared_ptr<Node>> _nodes;
  std::vector<std::shared_ptr<Node>> _internal_nodes;
  std::vector<std::shared_ptr<Node>> _external_nodes;
  size_t _internal_node_count;
  size_t _external_node_count;

  void processReRoot(std::shared_ptr<Node> node);
  void exchangeInfo(std::shared_ptr<Node> node1, std::shared_ptr<Node> node2);
  void postOrderProcessRoot(std::shared_ptr<Node> node);
  void setHeightFromRootToNode(std::shared_ptr<Node> inNode, double newHeight);
  double getGreatestDistance(std::shared_ptr<Node> inNode);

  bool findNode(std::shared_ptr<Node> n);

 public:
  Tree();
  explicit Tree(std::shared_ptr<Node> root);

  void addExternalNode(std::shared_ptr<Node> tn);
  void addInternalNode(std::shared_ptr<Node> tn);
  void pruneExternalNode(std::shared_ptr<Node> node);

  std::vector<std::shared_ptr<SplitOperation>> generateForwardOperations(
      Workspace &ws, const std::shared_ptr<MakeRateMatrixOperation> &rm);

  std::vector<std::shared_ptr<ReverseSplitOperation>>
  generateBackwardOperations(
      Workspace &ws, const std::shared_ptr<MakeRateMatrixOperation> &rm);

  std::vector<std::shared_ptr<SplitOperation>> generateForwardOperations(
      Workspace &ws, PeriodRateMatrixMap &period_rm_map,
      BranchProbMatrixMap &period_pm_map);
  std::vector<std::shared_ptr<ReverseSplitOperation>>
  generateBackwardOperations(Workspace &ws, PeriodRateMatrixMap &period_rm_map,
                             BranchProbMatrixMap &period_pm_map);

  std::vector<size_t> traversePreorderInternalNodesOnly() const;
  std::vector<size_t> traversePreorderInternalNodesOnlyNumbers() const;

  void assignTipData(
      Workspace &ws,
      const std::unordered_map<std::string, lagrange_dist_t> &dist_data);

  std::shared_ptr<Node> getExternalNode(size_t num);
  std::shared_ptr<Node> getExternalNode(const std::string &name);
  std::shared_ptr<Node> getInternalNode(size_t num);
  std::shared_ptr<Node> getInternalNode(const std::string &name);
  std::shared_ptr<Node> getNode(size_t num);

  unsigned int getNodeCount() const;
  unsigned int getExternalNodeCount() const;
  unsigned int getInternalNodeCount() const;

  std::shared_ptr<Node> getRoot();
  void processRoot();

  std::shared_ptr<Node> getMRCA(const std::vector<std::string> &innodes);
  std::shared_ptr<Node> getMRCA(
      const std::vector<std::shared_ptr<Node>> &innodes);

  std::shared_ptr<Node> getParent(std::shared_ptr<Node> n) const;

  void setHeightFromRootToNodes();
  void setHeightFromTipToNodes();

  std::string getNewick() const;
  std::string getNewickLambda(
      const std::function<std::string(const Node &)> &newick_lambda) const;

  ~Tree();
};

#endif /* TREE_H_ */
