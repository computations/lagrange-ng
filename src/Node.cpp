/*
 * node.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "Common.h"
#include "Node.h"
#include "Operation.h"

Node::Node()
    : _branch_length(0.0),
      _height(0.0),
      _number(0),
      _period(0),
      _label(""),
      _comment(""),
      _children{} {}

Node::Node(double bl, int innumber, const std::string &inname)
    : _branch_length(bl),
      _height(0.0),
      _number(innumber),
      _label(inname),
      _comment(""),
      _children{} {}

std::vector<std::shared_ptr<Node>> Node::getChildren() { return _children; }

bool Node::isExternal() const { return _children.size() == 0; }

bool Node::isInternal() const { return _children.size() != 0; }

int Node::getNumber() const { return _number; }

void Node::setNumber(int n) { _number = n; }

double Node::getBL() { return _branch_length; }

void Node::setBL(double bl) { _branch_length = bl; }

double Node::getHeight() { return _height; }

void Node::setHeight(double he) { _height = he; }

bool Node::hasChild(std::shared_ptr<Node> test) {
  for (unsigned int i = 0; i < _children.size(); i++) {
    if (_children.at(i) == test) { return true; }
  }
  return false;
}

bool Node::addChild(std::shared_ptr<Node> c) {
  if (hasChild(c) == false) {
    _children.push_back(c);
    return true;
  }
  return false;
}

bool Node::removeChild(std::shared_ptr<Node> c) {
  if (hasChild(c) == true) {
    for (auto it = _children.begin(); it != _children.end(); it++) {
      if (*it == c) {
        _children.erase(it);
        return true;
      }
    }
  }
  return false;
}

std::shared_ptr<Node> Node::getChild(int c) const { return _children.at(c); }

std::string Node::getName() const { return _label; }

void Node::setName(const std::string &s) { _label = s; }

std::string Node::getComment() const { return _comment; }

void Node::setComment(const std::string &s) { _comment = s; }

std::string Node::getNewick() const {
  static auto newick_lambda = [](const Node &n) { return n.getName(); };

  return getNewickLambda(newick_lambda);
}

std::string Node::getNewickLambda(
    const std::function<std::string(const Node &)> &newick_lambda) const {
  std::ostringstream newick_oss;
  for (int i = 0; i < getChildCount(); i++) {
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

int Node::getChildCount() const { return _children.size(); }

void Node::setSplitString(const std::string &splitstring) {
  _split_string = splitstring;
}

void Node::setStateString(const std::string &statestring) {
  _state_string = statestring;
}

void Node::setStochString(const std::string &stochstring) {
  _stoch_string = stochstring;
}

std::string Node::getStateString() const { return _state_string; }

std::string Node::getSplitString() const { return _split_string; }

std::string Node::getStochString() const { return _stoch_string; }

void Node::initExclDistVector() {
  _excluded_dists = std::make_shared<std::vector<lagrange_dist_t>>();
}

std::shared_ptr<std::vector<lagrange_dist_t>> &Node::getExclDistVector() {
  return _excluded_dists;
}

double Node::getMaxHeightRecursive() const {
  double max_height = 0.0;
  for (auto &c : _children) {
    max_height = std::max(max_height, c->getMaxHeightRecursive());
  }
  return max_height + _branch_length;
}

double Node::getMaxHeight() const {
  double max_height = 0.0;
  for (auto &c : _children) {
    max_height = std::max(max_height, c->getMaxHeightRecursive());
  }
  return max_height;
}

void Node::setHeightRecursive() {
  _height = getMaxHeight();
  for (auto &c : _children) { c->setHeightRecursive(); }
}

void Node::pruneNode(std::shared_ptr<Node> n) {
  for (auto it = _children.begin(); it != _children.end(); ++it) {
    if (*it == n) {
      _children.erase(it);
      if (_children.size() == 1) {
        auto child = _children[0];
        double tmp_branch_length = _branch_length + child->_branch_length;
        (*this) = *child;
        _branch_length = tmp_branch_length;
        setHeightRecursive();
      }
      return;
    }
  }
  for (auto &c : _children) { c->pruneNode(n); }
}

std::shared_ptr<Node> getMRCAWithNode(
    const std::shared_ptr<Node> &current,
    const std::vector<std::shared_ptr<Node>> &leaves) {
  if (current->_children.size() == 0) {
    for (auto &n : leaves) {
      if (n == current) { return current; }
    }
    return {nullptr};
  } else {
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
}

bool Node::findNode(std::shared_ptr<Node> n) {
  if (this == n.get()) { return true; }
  for (auto c : _children) {
    if (c->findNode(n)) { return true; }
  }
  return false;
}

std::shared_ptr<Node> getParentWithNode(const std::shared_ptr<Node> &current,
                                        const std::shared_ptr<Node> &n) {
  for (auto c : current->_children) {
    if (n == c) { return current; }
    auto ret = getParentWithNode(c, n);
    if (ret != nullptr) { return ret; }
  }
  return {nullptr};
}

void Node::traverseAndGenerateForwardOperations(
    Workspace &ws, PeriodRateMatrixMap &rm_map, BranchProbMatrixMap &pm_map,
    std::vector<Operation> &out_vector) const {
  if (_children.size() != 2 && _children.size() != 0) {
    throw std::runtime_error{
        "Tree is not bifircating when generating operations"};
  }

  /* Steps:
   * 0. Check that there are children.
   * 1. Generate Operations from Children
   * 2. Generate Expm Ops
   * 3. Generate Dispersion Operations from children
   * 4. Make split operation
   * 5. Bind into an Operation and add to vector
   */

  /* Step 0 */
  if (_children.size() == 0) {
    ws.register_top_clv(_id);
    return;
  }

  /* Step 1 */
  _children[0]->traverseAndGenerateForwardOperations(ws, rm_map, pm_map,
                                                     out_vector);
  _children[1]->traverseAndGenerateForwardOperations(ws, rm_map, pm_map,
                                                     out_vector);

  /* Steps 2 & 3 */
  std::vector<ExpmOperation> expm_ops;
  /* First the left child */

  ws.register_children_clv(_id);
  ws.register_top_clv(_id);

  std::vector<DispersionOperation> left_disp_ops;
  _children[0]->generateDispersionOperations(ws, rm_map, pm_map, left_disp_ops,
                                             expm_ops);
  left_disp_ops.back().terminate_top(getLeftCLVIndex(ws));

  /* Then then right child */
  std::vector<DispersionOperation> right_disp_ops;
  _children[1]->generateDispersionOperations(ws, rm_map, pm_map, right_disp_ops,
                                             expm_ops);
  right_disp_ops.back().terminate_top(getRightCLVIndex(ws));

  /* Step 4 */
  SplitOperation sp_op(getLeftCLVIndex(ws), getRightCLVIndex(ws),
                       getTopCLVIndex(ws));

  /* Step 5 */
  out_vector.emplace_back(sp_op, left_disp_ops, right_disp_ops, expm_ops);
}

void Node::traverseAndGenerateBackwardOperations(
    Workspace &ws, PeriodRateMatrixMap &rm_map, BranchProbMatrixMap &pm_map,
    std::vector<ReverseOperation> &out_vector) const {
  if (_children.size() != 2 && _children.size() != 0) {
    throw std::runtime_error{
        "Tree is not bifircating when generating operations"};
  }

  /* The reverse operation is just weird, and we need to generate 2 for every
   * node we visit. One of the two children will be the "parent" node. This will
   * generate the 2 operations that we said would be made.
   *
   * Steps:
   * 0. Check if there are children.
   * 1. Generate Expm Operations
   * 2. Generate Disperial Operations
   * 3. Generate Reverse Split Operation
   * 4. Terminate Dispersal Operations
   * 5. Generate Reverse Operations
   */

  /* Step 1 */
  if (_children.size() == 0) { return; }

  /* Steps 2 and 3 */
  std::vector<DispersionOperation> disp_ops;
  std::vector<ExpmOperation> expm_ops;

  ws.register_top_clv_reverse(_id);

  generateDispersionOperationsReverse(ws, rm_map, pm_map, disp_ops, expm_ops);

  /* Step 4 */
  ws.register_bot1_clv_reverse(_id);
  size_t dispersion_end_clv_index = disp_ops.back().top_clv_index();
  ReverseSplitOperation rsplit_op1(ws.get_bot1_clv_reverse(_id),
                                   dispersion_end_clv_index,
                                   ws.get_rchild_clv(_id));

  ws.register_bot2_clv_reverse(_id);
  ReverseSplitOperation rsplit_op2(ws.get_bot2_clv_reverse(_id),
                                   dispersion_end_clv_index,
                                   ws.get_lchild_clv(_id));

  out_vector.emplace_back(rsplit_op1, disp_ops, expm_ops);
  out_vector.emplace_back(rsplit_op2, disp_ops, expm_ops);

  if (_children[0]->isInternal()) {
    _children[0]->traverseAndGenerateBackwardOperations(ws, rm_map, pm_map,
                                                        out_vector);
    if (out_vector.back().is_bot_terminated()) {
      out_vector.back().terminate_bot(ws.get_bot1_clv_reverse(_id));
    }
  }

  if (_children[1]->isInternal()) {
    _children[1]->traverseAndGenerateBackwardOperations(ws, rm_map, pm_map,
                                                        out_vector);
    if (out_vector.back().is_bot_terminated()) {
      out_vector.back().terminate_bot(ws.get_bot1_clv_reverse(_id));
    }
  }
}

void Node::generateDispersionOperations(
    Workspace &ws, PeriodRateMatrixMap &rm_map, BranchProbMatrixMap &pm_map,
    std::vector<DispersionOperation> &disp_ops,
    std::vector<ExpmOperation> &expm_ops) const {
  auto pm = getProbMatrixOperation(ws, rm_map, pm_map);
  expm_ops.emplace_back(pm);
  disp_ops.emplace_back(std::numeric_limits<size_t>::max(), ws.get_top_clv(_id),
                        pm.get_index());
}

void Node::generateDispersionOperationsReverse(
    Workspace &ws, PeriodRateMatrixMap &rm_map, BranchProbMatrixMap &pm_map,
    std::vector<DispersionOperation> &disp_ops,
    std::vector<ExpmOperation> &expm_ops) const {
  auto pm = getProbMatrixOperation(ws, rm_map, pm_map, true);
  expm_ops.push_back(pm);
  DispersionOperation disp_op(ws.get_top_clv_reverse(_id),
                              std::numeric_limits<size_t>::max(),
                              pm.get_index());
  disp_ops.push_back(disp_op);
}

void Node::traverseAndGenerateBackwardNodeIds(std::vector<size_t> &ret) const {
  ret.push_back(_id);
  for (auto &c : _children) {
    c->traverseAndGenerateBackwardNodeIdsInternalOnly(ret);
  }
}

void Node::traverseAndGenerateBackwardNodeIdsInternalOnly(
    std::vector<size_t> &ret) const {
  ret.push_back(_id);
  for (auto &c : _children) {
    if (c->isInternal()) {
      c->traverseAndGenerateBackwardNodeIdsInternalOnly(ret);
    }
  }
}

void Node::traverseAndGeneratePostorderNodeIdsInternalOnly(
    std::vector<size_t> &ret) const {
  for (auto &c : _children) {
    if (c->isInternal()) {
      c->traverseAndGeneratePostorderNodeIdsInternalOnly(ret);
    }
  }
  ret.push_back(_id);
}

void Node::traverseAndGenerateBackwardNodeNumbersInternalOnly(
    std::vector<size_t> &ret) const {
  ret.push_back(_number);
  for (auto &c : _children) {
    if (c->isInternal()) {
      c->traverseAndGenerateBackwardNodeNumbersInternalOnly(ret);
    }
  }
}

void Node::assignTipData(
    Workspace &ws,
    const std::unordered_map<std::string, size_t> &distrib_data) const {
  if (_children.size() == 0) {
    ws.set_tip_clv(ws.get_top_clv(_id), distrib_data.at(_label));
  } else {
    for (auto &c : _children) { c->assignTipData(ws, distrib_data); }
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

size_t Node::getId() const { return _id; }

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
