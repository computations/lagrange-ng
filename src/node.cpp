/*
 * node.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <cstring>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <vector>
#include <unordered_map>

using namespace std;

#include "Common.h"
#include "node.h"

Node::Node()
    : _branch_length(0.0),
      _height(0.0),
      _number(0),
      _label(""),
      _comment(""),
      _children{} {}

Node::Node(double bl, int innumber, const string &inname)
    : _branch_length(bl),
      _height(0.0),
      _number(innumber),
      _label(inname),
      _comment(""),
      _children{} {}

vector<std::shared_ptr<Node>> Node::getChildren() { return _children; }

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
    if (_children.at(i) == test) {
      return true;
    }
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

string Node::getName() { return _label; }

void Node::setName(const string &s) { _label = s; }

string Node::getComment() { return _comment; }

void Node::setComment(const string &s) { _comment = s; }

string Node::getNewick(bool branch_lengths) const {
  std::ostringstream newick_oss;
  for (int i = 0; i < getChildCount(); i++) {
    if (i == 0) {
      newick_oss << "(";
    }
    newick_oss << getChild(i)->getNewick(branch_lengths);
    if (branch_lengths == true) {
      newick_oss << ":" << getChild(i)->getBL();
    }
    if (i == getChildCount() - 1) {
      newick_oss << ")";
    } else {
      newick_oss << ",";
    }
  }
  if (_label.size() > 0) {
    newick_oss << _label;
  }
  return newick_oss.str();
}

string Node::getNewick(
    bool branch_lengths,
    const std::function<string(const Node &)> &node_lambda) const {
  std::ostringstream newick_oss;
  for (int i = 0; i < getChildCount(); i++) {
    if (i == 0) {
      newick_oss << "(";
    }
    newick_oss << getChild(i)->getNewick(branch_lengths, node_lambda);
    if (branch_lengths == true) {
      newick_oss << ":" << getChild(i)->getBL();
    }
    if (i == getChildCount() - 1) {
      newick_oss << ")";
    } else {
      newick_oss << ",";
    }
  }
  if (isInternal() == true) {
    newick_oss << node_lambda(*this);
  } else {  // EXTERNAL
    if (_label.size() > 0) newick_oss << _label;
  }
  return newick_oss.str();
}

string Node::getNewickLambda(
    const std::function<string(const Node &)> &length_lambda) const {
  std::ostringstream newick_oss;
  for (int i = 0; i < getChildCount(); i++) {
    if (i == 0) {
      newick_oss << "(";
    }

    newick_oss << getChild(i)->getNewickLambda(length_lambda);
    newick_oss << length_lambda(*this);

    if (i == getChildCount() - 1) {
      newick_oss << ")";
    } else {
      newick_oss << ",";
    }
  }
  if (_label.size() > 0) {
    newick_oss << _label;
  }
  return newick_oss.str();
}

int Node::getChildCount() const { return _children.size(); }

void Node::setConditionalVector(const vector<Superdouble> &v) {
  _conditionals = v;
}

void Node::setReverseBits(const vector<Superdouble> &v) { _reverse_bits = v; }

const vector<Superdouble> &Node::getConditionalVector() const {
  return _conditionals;
}

const vector<Superdouble> &Node::getReverseBits() const {
  return _reverse_bits;
}

void Node::setSplitString(const string &splitstring) {
  _split_string = splitstring;
}

void Node::setStateString(const string &statestring) {
  _state_string = statestring;
}

void Node::setStochString(const string &stochstring) {
  _stoch_string = stochstring;
}

string Node::getStateString() const { return _state_string; }

string Node::getSplitString() const { return _split_string; }

string Node::getStochString() const { return _stoch_string; }

void Node::initSegVector() { _branch_segments = vector<BranchSegment>(); }

vector<BranchSegment> &Node::getSegVector() { return _branch_segments; }

void Node::initExclDistVector() {
  _excluded_dists = std::make_shared<vector<lagrange_dist_t>>();
}

std::shared_ptr<vector<lagrange_dist_t>> &Node::getExclDistVector() {
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
  for (auto &c : _children) {
    c->setHeightRecursive();
  }
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
  for (auto &c : _children) {
    c->pruneNode(n);
  }
}

std::shared_ptr<Node> getMRCAWithNode(
    const std::shared_ptr<Node> &current,
    const std::vector<std::shared_ptr<Node>> &leaves) {
  if (current->_children.size() == 0) {
    for (auto &n : leaves) {
      if (n == current) {
        return current;
      }
    }
    return {nullptr};
  } else {
    std::shared_ptr<Node> mrca = nullptr;
    for (auto &c : current->_children) {
      auto tmp_mrca = getMRCAWithNode(c, leaves);
      if (tmp_mrca != nullptr) {
        // if this is our second match, then we can return with this pointer.
        if (mrca != nullptr) {
          return current;
        }
        mrca = tmp_mrca;
      }
    }
    return mrca;
  }
}

bool Node::findNode(std::shared_ptr<Node> n) {
  if (this == n.get()) {
    return true;
  }
  for (auto c : _children) {
    if (c->findNode(n)) {
      return true;
    }
  }
  return false;
}

std::shared_ptr<Node> getParentWithNode(const std::shared_ptr<Node> &current,
                                        const std::shared_ptr<Node> &n) {
  for (auto c : current->_children) {
    if (n == c) {
      return current;
    }
    auto ret = getParentWithNode(c, n);
    if (ret != nullptr) {
      return ret;
    }
  }
  return {nullptr};
}

std::pair<std::vector<SplitOperation>, std::shared_ptr<DispersionOperation>>
Node::traverseAndGenerateForwardOperations(
    Workspace &ws,
    const std::shared_ptr<MakeRateMatrixOperation> &rm_op) const {
  if (_children.size() != 2 && _children.size() != 0) {
    throw std::runtime_error{
        "Tree is not bifircating when generating operations"};
  }

  if (_children.size() == 0) {
    ws.register_top_clv(_id);
    return {{}, generateDispersionOperations(ws, rm_op)};
  }

  std::vector<SplitOperation> split_ops;

  auto lchild = _children[0]->traverseAndGenerateForwardOperations(ws, rm_op);
  auto rchild = _children[1]->traverseAndGenerateForwardOperations(ws, rm_op);

  split_ops.reserve(lchild.first.size() + rchild.first.size() + 1);
  split_ops.insert(split_ops.end(), lchild.first.begin(), lchild.first.end());
  split_ops.insert(split_ops.end(), rchild.first.begin(), rchild.first.end());

  ws.register_top_clv(_id);
  ws.register_children_clv(_id);
  lchild.second->terminate_top(ws.get_lchild_clv(_id));
  rchild.second->terminate_top(ws.get_rchild_clv(_id));

  split_ops.emplace_back(ws.get_top_clv(_id), lchild.second, rchild.second);
  return {split_ops, generateDispersionOperations(ws, rm_op)};
}

std::pair<std::vector<ReverseSplitOperation>,
          std::shared_ptr<DispersionOperation>>
Node::traverseAndGenerateBackwardOperations(
    Workspace &ws,
    const std::shared_ptr<MakeRateMatrixOperation> &rm_op) const {
  if (_children.size() != 2 && _children.size() != 0) {
    throw std::runtime_error{
        "Tree is not bifircating when generating operations"};
  }

  if (_children.size() == 0) {
    return {{}, {}};
  }

  std::vector<ReverseSplitOperation> rsplit_ops;

  ws.register_top_clv_reverse(_id);

  auto disp_ops = generateDispersionOperationsReverse(ws, rm_op);

  ws.register_bot1_clv_reverse(_id);
  rsplit_ops.emplace_back(ws.get_bot1_clv_reverse(_id), ws.get_rchild_clv(_id),
                          disp_ops);

  ws.register_bot2_clv_reverse(_id);
  rsplit_ops.emplace_back(ws.get_bot2_clv_reverse(_id), ws.get_lchild_clv(_id),
                          disp_ops);

  if (_children[0]->isInternal()) {
    auto child_trav =
        _children[0]->traverseAndGenerateBackwardOperations(ws, rm_op);
    auto child_disp_op = child_trav.second;

    child_disp_op->terminate_bot(ws.get_bot1_clv_reverse(_id));
    rsplit_ops.insert(rsplit_ops.end(), child_trav.first.begin(),
                      child_trav.first.end());
  }

  if (_children[1]->isInternal()) {
    auto child_trav =
        _children[1]->traverseAndGenerateBackwardOperations(ws, rm_op);
    auto child_disp_op = child_trav.second;

    child_disp_op->terminate_bot(ws.get_bot2_clv_reverse(_id));
    rsplit_ops.insert(rsplit_ops.end(), child_trav.first.begin(),
                      child_trav.first.end());
  }
  return {rsplit_ops, disp_ops};
}

std::shared_ptr<DispersionOperation> Node::generateDispersionOperations(
    Workspace &ws,
    const std::shared_ptr<MakeRateMatrixOperation> &rm_op) const {
  return std::make_shared<DispersionOperation>(
      ws.get_top_clv(_id), _branch_length, ws.suggest_prob_matrix_index(),
      rm_op);
}

std::shared_ptr<DispersionOperation> Node::generateDispersionOperationsReverse(
    Workspace &ws,
    const std::shared_ptr<MakeRateMatrixOperation> &rm_op) const {
  return std::make_shared<DispersionOperation>(
      ws.get_top_clv_reverse(_id), std::numeric_limits<size_t>::max(),
      _branch_length, ws.suggest_prob_matrix_index(), rm_op,
      /*transpose=*/true);
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

void Node::assignTipData(Workspace &ws,
                         const std::unordered_map<std::string, lagrange_dist_t>
                             &distrib_data) const {
  if (_children.size() == 0) {
    ws.set_tip_clv(ws.get_top_clv(_id), distrib_data.at(_label));
  } else {
    for (auto &c : _children) {
      c->assignTipData(ws, distrib_data);
    }
  }
}

void Node::assignIdRecursive(size_t &id) {
  _id = id++;
  for (auto &c : _children) {
    c->assignIdRecursive(id);
  }
}

void Node::assignId() {
  size_t id = 0;
  assignIdRecursive(id);
}

size_t Node::getId() const { return _id; }
