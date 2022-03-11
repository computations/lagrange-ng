/*
 * node.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Common.h"
#include "Node.h"
#include "Operation.h"

Node::Node()
    : _branch_length(0.0), _height(0.0), _number(0), _id(0), _period(0) {}

Node::Node(double bl, size_t innumber, std::string inname)
    : _branch_length(bl),
      _height(0.0),
      _number(innumber),
      _id(0),
      _period(0),
      _label(std::move(inname)) {}

auto Node::isExternal() const -> bool { return _children.empty(); }

auto Node::isInternal() const -> bool { return !_children.empty(); }

auto Node::getNumber() const -> size_t { return _number; }

void Node::setNumber(size_t n) { _number = n; }

auto Node::getBL() const -> double { return _branch_length; }

void Node::setBL(double bl) { _branch_length = bl; }

void Node::setHeight(double he) { _height = he; }

auto Node::hasChild(const std::shared_ptr<Node> &test) -> bool {
  return std::any_of(_children.begin(), _children.end(),
                     [&test](auto n) { return n == test; });
}

auto Node::addChild(const std::shared_ptr<Node> &c) -> bool {
  if (!hasChild(c)) {
    _children.push_back(c);
    return true;
  }
  return false;
}

auto Node::removeChild(const std::shared_ptr<Node> &c) -> bool {
  if (hasChild(c)) {
    for (auto it = _children.begin(); it != _children.end(); it++) {
      if (*it == c) {
        _children.erase(it);
        return true;
      }
    }
  }
  return false;
}

auto Node::getChild(size_t c) const -> std::shared_ptr<Node> {
  return _children.at(c);
}

auto Node::getName() const -> std::string { return _label; }

void Node::setName(const std::string &s) { _label = s; }

void Node::setComment(const std::string &s) { _comment = s; }

auto Node::getNewick() const -> std::string {
  static auto newick_lambda = [](const Node &n) { return n.getName(); };

  return getNewickLambda(newick_lambda);
}

auto Node::getNewickLambda(
    const std::function<std::string(const Node &)> &newick_lambda) const
    -> std::string {
  std::ostringstream newick_oss;
  for (size_t i = 0; i < getChildCount(); i++) {
    if (i == 0) { newick_oss << "("; }

    newick_oss << getChild(i)->getNewickLambda(newick_lambda);

    if (i == getChildCount() - 1) {
      newick_oss << ")";
    } else {
      newick_oss << ",";
    }
  }
  newick_oss << newick_lambda(*this);
  return newick_oss.str();
}

auto Node::getChildCount() const -> size_t { return _children.size(); }

void Node::setSplitString(const std::string &splitstring) {
  _split_string = splitstring;
}

void Node::setStateString(const std::string &statestring) {
  _state_string = statestring;
}

void Node::setStochString(const std::string &stochstring) {
  _stoch_string = stochstring;
}

auto Node::getStateString() const -> std::string { return _state_string; }

auto Node::getSplitString() const -> std::string { return _split_string; }

auto Node::getStochString() const -> std::string { return _stoch_string; }

void Node::initExclDistVector() {
  _excluded_dists = std::make_shared<std::vector<lagrange_dist_t>>();
}

auto Node::getMaxHeightRecursive() const -> double {
  double max_height = 0.0;
  for (const auto &c : _children) {
    max_height = std::max(max_height, c->getMaxHeightRecursive());
  }
  return max_height + _branch_length;
}

auto Node::getMaxHeight() const -> double {
  double max_height = 0.0;
  for (const auto &c : _children) {
    max_height = std::max(max_height, c->getMaxHeightRecursive());
  }
  return max_height;
}

void Node::setHeightRecursive() {
  _height = getMaxHeight();
  for (auto &c : _children) { c->setHeightRecursive(); }
}

auto getMRCAWithNode(const std::shared_ptr<Node> &current,
                     const std::vector<std::shared_ptr<Node>> &leaves)
    -> std::shared_ptr<Node> {
  if (current->_children.empty()) {
    if (std::any_of(current->_children.begin(), current->_children.end(),
                    [&current](auto &n) { return n == current; })) {
      return current;
    }
    return {nullptr};
  }
  std::shared_ptr<Node> mrca = nullptr;
  for (auto &c : current->_children) {
    auto tmp_mrca = getMRCAWithNode(c, leaves);
    if (tmp_mrca != nullptr) {
      // if this is our second match, then we can return with this pointer.
      if (mrca != nullptr) { return current; }
      mrca = tmp_mrca;
    }
  }
  return mrca;
}

auto Node::findNode(const std::shared_ptr<Node> &n) -> bool {
  if (this == n.get()) { return true; }
  return std::any_of(_children.begin(), _children.end(),
                     [&n](auto &c) { return c->findNode(n); });
}

auto getParentWithNode(const std::shared_ptr<Node> &current,
                       const std::shared_ptr<Node> &n)
    -> std::shared_ptr<Node> {
  for (const auto &c : current->_children) {
    if (n == c) { return current; }
    auto ret = getParentWithNode(c, n);
    if (ret != nullptr) { return ret; }
  }
  return {nullptr};
}

auto Node::traverseAndGenerateForwardOperations(
    Workspace &ws, PeriodRateMatrixMap &pm_map,
    BranchProbMatrixMap &bm_map) const
    -> std::pair<std::vector<std::shared_ptr<SplitOperation>>,
                 std::shared_ptr<DispersionOperation>> {
  if (_children.size() != 2 && !_children.empty()) {
    throw std::runtime_error{
        "Tree is not bifircating when generating operations"};
  }

  if (_children.empty()) {
    ws.register_top_clv(_id);
    return {{}, generateDispersionOperations(ws, pm_map, bm_map)};
  }

  std::vector<std::shared_ptr<SplitOperation>> split_ops;

  auto lchild =
      _children[0]->traverseAndGenerateForwardOperations(ws, pm_map, bm_map);
  auto rchild =
      _children[1]->traverseAndGenerateForwardOperations(ws, pm_map, bm_map);

  split_ops.reserve(lchild.first.size() + rchild.first.size() + 1);
  split_ops.insert(split_ops.end(), lchild.first.begin(), lchild.first.end());
  split_ops.insert(split_ops.end(), rchild.first.begin(), rchild.first.end());

  ws.register_top_clv(_id);
  ws.register_children_clv(_id);
  lchild.second->terminate_top(ws.get_lchild_clv(_id));
  rchild.second->terminate_top(ws.get_rchild_clv(_id));

  split_ops.push_back(std::make_shared<SplitOperation>(
      ws.get_top_clv(_id), lchild.second, rchild.second));
  return {split_ops, generateDispersionOperations(ws, pm_map, bm_map)};
}

auto Node::traverseAndGenerateBackwardOperations(
    Workspace &ws, PeriodRateMatrixMap &rm_map,
    BranchProbMatrixMap &pm_map) const
    -> std::pair<std::vector<std::shared_ptr<ReverseSplitOperation>>,
                 std::shared_ptr<DispersionOperation>> {
  if (_children.size() != 2 && !_children.empty()) {
    throw std::runtime_error{
        "Tree is not bifircating when generating operations"};
  }

  if (_children.empty()) { return {{}, {}}; }

  std::vector<std::shared_ptr<ReverseSplitOperation>> rsplit_ops;

  ws.register_top_clv_reverse(_id);

  auto disp_ops = generateDispersionOperationsReverse(ws, rm_map, pm_map);

  ws.register_bot1_clv_reverse(_id);
  rsplit_ops.push_back(std::make_shared<ReverseSplitOperation>(
      ws.get_bot1_clv_reverse(_id), ws.get_rchild_clv(_id), disp_ops));

  ws.register_bot2_clv_reverse(_id);
  rsplit_ops.push_back(std::make_shared<ReverseSplitOperation>(
      ws.get_bot2_clv_reverse(_id), ws.get_lchild_clv(_id), disp_ops));

  if (_children[0]->isInternal()) {
    auto child_trav =
        _children[0]->traverseAndGenerateBackwardOperations(ws, rm_map, pm_map);
    auto child_disp_op = child_trav.second;

    child_disp_op->terminate_bot(ws.get_bot1_clv_reverse(_id));
    rsplit_ops.insert(rsplit_ops.end(), child_trav.first.begin(),
                      child_trav.first.end());
  }

  if (_children[1]->isInternal()) {
    auto child_trav =
        _children[1]->traverseAndGenerateBackwardOperations(ws, rm_map, pm_map);
    auto child_disp_op = child_trav.second;

    child_disp_op->terminate_bot(ws.get_bot2_clv_reverse(_id));
    rsplit_ops.insert(rsplit_ops.end(), child_trav.first.begin(),
                      child_trav.first.end());
  }
  return {rsplit_ops, disp_ops};
}

auto Node::generateDispersionOperations(Workspace &ws,
                                        PeriodRateMatrixMap &rm_map,
                                        BranchProbMatrixMap &pm_map) const
    -> std::shared_ptr<DispersionOperation> {
  return std::make_shared<DispersionOperation>(
      ws.get_top_clv(_id), getProbMatrixOperation(ws, rm_map, pm_map));
}

auto Node::generateDispersionOperationsReverse(
    Workspace &ws, PeriodRateMatrixMap &rm_map,
    BranchProbMatrixMap &pm_map) const -> std::shared_ptr<DispersionOperation> {
  return std::make_shared<DispersionOperation>(
      ws.get_top_clv_reverse(_id), std::numeric_limits<size_t>::max(),
      getProbMatrixOperation(ws, rm_map, pm_map, true));
}

void Node::traverseAndGenerateBackwardNodeIdsInternalOnly(
    std::vector<size_t> &ret) const {
  ret.push_back(_id);
  for (const auto &c : _children) {
    if (c->isInternal()) {
      c->traverseAndGenerateBackwardNodeIdsInternalOnly(ret);
    }
  }
}

void Node::traverseAndGeneratePostorderNodeIdsInternalOnly(
    std::vector<size_t> &ret) const {
  for (const auto &c : _children) {
    if (c->isInternal()) {
      c->traverseAndGeneratePostorderNodeIdsInternalOnly(ret);
    }
  }
  ret.push_back(_id);
}

void Node::traverseAndGenerateBackwardNodeNumbersInternalOnly(
    std::vector<size_t> &ret) const {
  ret.push_back(_number);
  for (const auto &c : _children) {
    if (c->isInternal()) {
      c->traverseAndGenerateBackwardNodeNumbersInternalOnly(ret);
    }
  }
}

void Node::assignTipData(
    Workspace &ws,
    const std::unordered_map<std::string, size_t> &distrib_data) const {
  if (_children.empty()) {
    ws.set_tip_clv(ws.get_top_clv(_id), distrib_data.at(_label));
  } else {
    for (const auto &c : _children) { c->assignTipData(ws, distrib_data); }
  }
}

void Node::assignIdRecursive(size_t &id) {
  _id = id++;
  for (auto &c : _children) { c->assignIdRecursive(id); }
}

void Node::assignId() {
  size_t id = 0;
  assignIdRecursive(id);
}

auto Node::getId() const -> size_t { return _id; }

void Node::setSplitStringRecursive(
    const std::vector<size_t> &id_map,
    const std::vector<lagrange_col_vector_t> &dist_lhs, size_t states,
    const std::vector<std::string> &names) {
  if (isExternal()) { return; }

  lagrange_dist_t best_dist =
      lagrange_compute_best_dist(dist_lhs[id_map[getNumber()]], states);
  _split_string = lagrange_convert_dist_string(best_dist, names);
  for (auto &c : _children) {
    c->setSplitStringRecursive(id_map, dist_lhs, states, names);
  }
}

void Node::setStateStringRecursive(
    const std::vector<size_t> &id_map,
    const std::vector<lagrange_col_vector_t> &dist_lhs, size_t states,
    const std::vector<std::string> &names) {
  if (isExternal()) { return; }

  lagrange_dist_t best_dist =
      lagrange_compute_best_dist(dist_lhs[id_map[getNumber()]], states);
  _state_string = lagrange_convert_dist_string(best_dist, names);
  for (auto &c : _children) {
    c->setStateStringRecursive(id_map, dist_lhs, states, names);
  }
}
