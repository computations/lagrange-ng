/*
 * tree.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include <math.h>

#include <exception>
#include <limits>
#include <unordered_map>
#include <utility>

#include "Operation.h"
#include "Tree.h"

Tree::Tree() : Tree(nullptr) {}

Tree::Tree(std::shared_ptr<Node> inroot)
    : _root(std::move(inroot)),
      _internal_node_count(0),
      _external_node_count(0) {
  processRoot();
  _root->assignId();
  for (unsigned int i = 0; i < getNodeCount(); i++) {
    getNode(i)->initExclDistVector();
  }

  setHeightFromTipToNodes();
}

void Tree::addExternalNode(const std::shared_ptr<Node> &tn) {
  _external_nodes.push_back(tn);
  _external_node_count = _external_nodes.size();
  _nodes.push_back(tn);
}

void Tree::addInternalNode(const std::shared_ptr<Node> &tn) {
  _internal_nodes.push_back(tn);
  _internal_node_count = _internal_nodes.size();
  _nodes.push_back(tn);
}

auto Tree::getExternalNode(size_t num) -> std::shared_ptr<Node> {
  return _external_nodes.at(num);
}

/*
 * could precompute this, check for run time differences
 */
auto Tree::getExternalNode(const std::string &name) -> std::shared_ptr<Node> {
  std::shared_ptr<Node> ret = nullptr;
  for (auto &_external_node : _external_nodes) {
    if (_external_node->getName() == name) { ret = _external_node; }
  }
  return ret;
}

auto Tree::getInternalNode(size_t num) -> std::shared_ptr<Node> {
  return _internal_nodes.at(num);
}

/*
 * could precompute this, check for run time differences
 */
auto Tree::getInternalNode(const std::string &name) -> std::shared_ptr<Node> {
  std::shared_ptr<Node> ret = nullptr;
  for (auto &_internal_node : _internal_nodes) {
    if (_internal_node->getName() == name) { ret = _internal_node; }
  }
  return ret;
}

auto Tree::getExternalNodeCount() const -> unsigned int {
  return _external_node_count;
}

auto Tree::getInternalNodeCount() const -> unsigned int {
  return _internal_node_count;
}

auto Tree::getNode(size_t num) -> std::shared_ptr<Node> {
  return _nodes.at(num);
}

auto Tree::getNodeCount() const -> unsigned int { return _nodes.size(); }

auto Tree::getRoot() -> std::shared_ptr<Node> { return _root; }

auto Tree::getMRCA(const std::vector<std::string> &outgroup)
    -> std::shared_ptr<Node> {
  if (outgroup.size() == 1) { return getExternalNode(outgroup[0]); }
  std::shared_ptr<Node> mcra;
  std::vector<std::shared_ptr<Node>> outgroup_nodes;
  outgroup_nodes.reserve(outgroup.size());
  for (const auto &o : outgroup) {
    outgroup_nodes.push_back(getExternalNode(o));
  }
  return getMRCA(outgroup_nodes);
}

auto Tree::getMRCA(const std::vector<std::shared_ptr<Node>> &outgroup)
    -> std::shared_ptr<Node> {
  return getMRCAWithNode(_root, outgroup);
}

void Tree::setHeightFromRootToNodes() {
  setHeightFromRootToNode(_root, _root->getBL());
}

void Tree::setHeightFromRootToNode(const std::shared_ptr<Node> &inNode,
                                   double newHeight) {
  if (inNode != _root) {
    newHeight += inNode->getBL();
    inNode->setHeight(newHeight);
  } else {
    inNode->setHeight(newHeight);
  }
  for (size_t i = 0; i < inNode->getChildCount(); i++) {
    setHeightFromRootToNode(inNode->getChild(i), newHeight);
  }
}

/*
 * only makes sense for ultrametric trees
 */
void Tree::setHeightFromTipToNodes() { _root->setHeightRecursive(); }

/*
 * private
 */
void Tree::processRoot() {
  _nodes.clear();
  _internal_nodes.clear();
  _external_nodes.clear();
  _internal_node_count = 0;
  _external_node_count = 0;
  if (_root == nullptr) { return; }
  postOrderProcessRoot(_root);
}

void Tree::processReRoot(const std::shared_ptr<Node> &node) {
  if (node != _root || node->isExternal()) { return; }
  if (getParent(node) != nullptr) { processReRoot(getParent(node)); }
  // Exchange branch label, length et cetera
  exchangeInfo(getParent(node), node);
  // Rearrange topology
  std::shared_ptr<Node> parent = getParent(node);
  node->addChild(parent);
  parent->removeChild(node);
}

void Tree::exchangeInfo(const std::shared_ptr<Node> &node1,
                        const std::shared_ptr<Node> &node2) {
  std::string swaps;
  double swapd = NAN;
  swaps = node1->getName();
  node1->setName(node2->getName());
  node2->setName(swaps);

  swapd = node1->getBL();
  node1->setBL(node2->getBL());
  node2->setBL(swapd);
}

void Tree::postOrderProcessRoot(const std::shared_ptr<Node> &node) {
  if (node == nullptr) { return; }
  if (node->getChildCount() > 0) {
    for (size_t i = 0; i < node->getChildCount(); i++) {
      postOrderProcessRoot(node->getChild(i));
    }
  }
  if (node->isExternal()) {
    addExternalNode(node);
    node->setNumber(_external_node_count);
  } else {
    addInternalNode(node);
    node->setNumber(_internal_node_count);
  }
}

/*
 * end private
 */

Tree::~Tree() { _root.reset(); }

auto Tree::findNode(std::shared_ptr<Node> n) -> bool {
  return _root->findNode(std::move(n));
}

auto Tree::getParent(const std::shared_ptr<Node> &n) const
    -> std::shared_ptr<Node> {
  return getParentWithNode(_root, n);
}

auto Tree::generateForwardOperations(
    Workspace &ws, const std::shared_ptr<MakeRateMatrixOperation> &rm)
    -> std::vector<std::shared_ptr<SplitOperation>> {
  PeriodRateMatrixMap rm_map;
  BranchProbMatrixMap pm_map;
  rm_map[0] = rm;
  return generateForwardOperations(ws, rm_map, pm_map);
}

auto Tree::generateForwardOperations(Workspace &ws, PeriodRateMatrixMap &rm_map,
                                     BranchProbMatrixMap &pm_map)
    -> std::vector<std::shared_ptr<SplitOperation>> {
  return _root->traverseAndGenerateForwardOperations(ws, rm_map, pm_map).first;
}

auto Tree::generateBackwardOperations(
    Workspace &ws, const std::shared_ptr<MakeRateMatrixOperation> &rm)
    -> std::vector<std::shared_ptr<ReverseSplitOperation>> {
  PeriodRateMatrixMap rm_map;
  BranchProbMatrixMap pm_map;
  rm_map[0] = rm;
  return generateBackwardOperations(ws, rm_map, pm_map);
}

auto Tree::generateBackwardOperations(Workspace &ws,
                                      PeriodRateMatrixMap &rm_map,
                                      BranchProbMatrixMap &pm_map)
    -> std::vector<std::shared_ptr<ReverseSplitOperation>> {
  auto ret = _root->traverseAndGenerateBackwardOperations(ws, rm_map, pm_map);
  (*ret.first.begin())
      ->makeRootOperation(ws.get_bot1_clv_reverse(_root->getId()));
  (*(ret.first.begin() + 1))
      ->makeRootOperation(ws.get_bot2_clv_reverse(_root->getId()));

  return ret.first;
}

auto Tree::traversePreorderInternalNodesOnly() const -> std::vector<size_t> {
  std::vector<size_t> ret;
  ret.reserve(getInternalNodeCount());
  _root->traverseAndGenerateBackwardNodeIdsInternalOnly(ret);
  return ret;
}

auto Tree::traversePreorderInternalNodesOnlyNumbers() const
    -> std::vector<size_t> {
  std::vector<size_t> ret;
  ret.reserve(getInternalNodeCount());
  _root->traverseAndGenerateBackwardNodeNumbersInternalOnly(ret);
  return ret;
}

void Tree::assignTipData(
    Workspace &ws,
    const std::unordered_map<std::string, lagrange_dist_t> &dist_data) {
  _root->assignTipData(ws, dist_data);
}

auto Tree::getNewick() const -> std::string { return _root->getNewick() + ";"; }

auto Tree::getNewickLambda(
    const std::function<std::string(const Node &)> &newick_lambda) const
    -> std::string {
  return _root->getNewickLambda(newick_lambda) + ";";
}
