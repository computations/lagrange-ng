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
#include "Operation.h"
#include "node.h"

class Tree {
 private:
  std::shared_ptr<Node> _root;
  vector<std::shared_ptr<Node>> _nodes;
  vector<std::shared_ptr<Node>> _internal_nodes;
  vector<std::shared_ptr<Node>> _external_nodes;
  unsigned int _internal_node_count;
  unsigned int _external_node_count;

  void processReRoot(std::shared_ptr<Node> node);
  void exchangeInfo(std::shared_ptr<Node> node1, std::shared_ptr<Node> node2);
  void postOrderProcessRoot(std::shared_ptr<Node> node);
  void setHeightFromRootToNode(std::shared_ptr<Node> inNode, double newHeight);
  double getGreatestDistance(std::shared_ptr<Node> inNode);

  bool findNode(std::shared_ptr<Node> n);

 public:
  Tree();
  Tree(std::shared_ptr<Node> root);

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

  std::shared_ptr<Node> getExternalNode(int num);
  std::shared_ptr<Node> getExternalNode(string &name);
  std::shared_ptr<Node> getInternalNode(int num);
  std::shared_ptr<Node> getInternalNode(string &name);
  std::shared_ptr<Node> getNode(int num);

  unsigned int getNodeCount() const;
  unsigned int getExternalNodeCount() const;
  unsigned int getInternalNodeCount() const;

  std::shared_ptr<Node> getRoot();
  void setRoot(std::shared_ptr<Node> inroot);
  void unRoot(std::shared_ptr<Node> inroot);
  void reRoot(std::shared_ptr<Node> inroot);
  void tritomyRoot(std::shared_ptr<Node> toberoot);
  void processRoot();

  std::shared_ptr<Node> getMRCA(vector<string> innodes);
  std::shared_ptr<Node> getMRCA(vector<std::shared_ptr<Node>> innodes);

  double getLongestPathRootToTip() const;
  std::shared_ptr<Node> getParent(std::shared_ptr<Node> n) const;

  void setHeightFromRootToNodes();
  void setHeightFromTipToNodes();

  std::string getNewick() const;
  std::string getNewickLambda(
      const std::function<string(const Node &)> &newick_lambda) const;

  void setStateStrings(const std::vector<size_t> &id_map,
                       const std::vector<lagrange_col_vector_t> &dist_lhs,
                       const std::vector<std::string> &names);
  void setSplitStrings(const std::vector<size_t> &id_map,
                       const std::vector<lagrange_col_vector_t> &dist_lhs,
                       const std::vector<std::string> &names);

  ~Tree();
};

#endif /* TREE_H_ */
