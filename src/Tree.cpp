/*
 * tree.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include "Tree.hpp"

#include <cmath>
#include <sstream>
#include <unordered_map>
#include <utility>

#include "MRCA.hpp"
#include "Node.hpp"
#include "Operation.hpp"
#include "logger.hpp"

namespace lagrange {

Tree::Tree() : Tree(nullptr) {}

Tree::Tree(std::shared_ptr<Node> inroot) : _root(std::move(inroot)) {
  _root->assignId();
  _node_count = _root->getCount();
  _tip_count = _root->getTipCount();
}

auto Tree::getNodeCount() const -> size_t { return _node_count; }

auto Tree::getTipCount() const -> size_t { return _tip_count; }

auto Tree::getInternalCount() const -> size_t {
  return _node_count - _tip_count;
}

auto Tree::getNode(size_t id) -> std::shared_ptr<Node> {
  return getNodeById(_root, id);
}

auto Tree::getRoot() -> std::shared_ptr<Node> { return _root; }

auto Tree::getMRCA(const std::shared_ptr<MRCAEntry> &mrca)
    -> std::shared_ptr<Node> {
  std::vector<std::shared_ptr<Node>> members;
  getNodesByMRCAEntry(_root, mrca, members);
  return getMRCAWithNodes(_root, members);
}

void Tree::setHeightTopDown() { _root->setHeightReverse(); }

/*
 * only makes sense for ultrametric trees
 */
auto Tree::setHeightBottomUp() -> double { return _root->setHeight(); }

  void Tree::scaleBranchLengths(double scale){
  _root->scaleBranchLength(scale);
}

void Tree::processReRoot(const std::shared_ptr<Node> &node) {
  if (node != _root || node->isTip()) { return; }
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

Tree::~Tree() { _root.reset(); }

auto Tree::findNode(const std::shared_ptr<Node> &n) -> bool {
  return _root->findNode(n);
}

auto Tree::getParent(const std::shared_ptr<Node> &n) const
    -> std::shared_ptr<Node> {
  return getParentWithNode(_root, n);
}

auto Tree::generateForwardOperations(Workspace &ws)
    -> std::vector<std::shared_ptr<SplitOperation>> {
  PeriodRateMatrixMap rm_map;
  BranchProbMatrixMap pm_map;
  return generateForwardOperations(ws, rm_map, pm_map);
}

auto Tree::generateForwardOperations(Workspace &ws,
                                     PeriodRateMatrixMap &rm_map,
                                     BranchProbMatrixMap &pm_map)
    -> std::vector<std::shared_ptr<SplitOperation>> {
  return _root->traverseAndGenerateForwardOperations(ws, rm_map, pm_map).first;
}

auto Tree::generateBackwardOperations(Workspace &ws)
    -> std::vector<std::shared_ptr<ReverseSplitOperation>> {
  PeriodRateMatrixMap rm_map;
  BranchProbMatrixMap pm_map;
  return generateBackwardOperations(ws, rm_map, pm_map);
}

auto Tree::generateBackwardOperations(Workspace &ws,
                                      PeriodRateMatrixMap &rm_map,
                                      BranchProbMatrixMap &pm_map)
    -> std::vector<std::shared_ptr<ReverseSplitOperation>> {
  auto ret =
      _root->traverseAndGenerateBackwardOperations(ws, rm_map, pm_map, true);
  (ret.first[0])->makeRootOperation(ws.getBot1CLVReverse(_root->getId()));
  (ret.first[1])->makeRootOperation(ws.getBot2CLVReverse(_root->getId()));

  return ret.first;
}

auto Tree::traversePreorderInternalNodesOnly() const -> std::vector<size_t> {
  std::vector<size_t> ret;
  ret.reserve(getInternalCount());
  _root->traverseAndGenerateBackwardNodeIdsInternalOnly(ret);
  return ret;
}

auto Tree::traversePreorderInternalNodesOnlyNumbers() const
    -> std::vector<size_t> {
  std::vector<size_t> ret;
  ret.reserve(getInternalCount());
  _root->traverseAndGenerateBackwardNodeNumbersInternalOnly(ret);
  return ret;
}

void Tree::assignTipData(
    Workspace &ws, const std::unordered_map<std::string, Range> &dist_data) {
  _root->assignTipData(ws, dist_data);
}

auto Tree::getNewick() const -> std::string { return _root->getNewick() + ";"; }

auto Tree::getNewickLambda(const std::function<std::string(const Node &)>
                               &newick_lambda) const -> std::string {
  return _root->getNewickLambda(newick_lambda) + ";";
}

bool Tree::checkAlignmentConsistency(const Alignment &align) const {
  auto count = _root->checkAlignmentConsistency(align, 0);
  if (count != align.taxa_count) {
    std::ostringstream oss;
    oss << "Taxa present in alignment that are not presesnt on the tree";
    throw std::runtime_error{oss.str()};
  }
  return true;
}

void Tree::assignFossils(const std::vector<Fossil> &fossils) {
  for (const auto &f : fossils) { _root->assignFossil(f); }
}

void Tree::applyPreorderInternalOnly(const std::function<void(Node &)> &func) {
  _root->applyPreorderInternalOnly(func);
}

void Tree::setPeriods(const Periods &periods) {
  auto period_func = [&](Node &n) { n.setPeriodSegments(periods); };
  _root->applyPreorder(period_func);
}

void Tree::assignMCRALabels(const MRCAMap &mrca_map) {
  for (const auto &kv : mrca_map) {
    auto n = getMRCA(kv.second);
    if (!n) {
      LOG_ERROR(
          "MRCA '%s' not found, please check that the tips exist in the tree",
          kv.first.c_str());
    }
    n->setMRCALabel(kv.first);
  }
}

void Tree::assignStateResult(std::unique_ptr<LagrangeMatrixBase[]> r,
                             const MRCALabel &mrca_label) {
  getNodesByMRCALabel(_root, mrca_label)->assignAncestralState(std::move(r));
}

void Tree::assignSplitResult(const SplitReturn &r,
                             const MRCALabel &mrca_label) {
  getNodesByMRCALabel(_root, mrca_label)->assignAncestralSplit(r);
}
}  // namespace lagrange
