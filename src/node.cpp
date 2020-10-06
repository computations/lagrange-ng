/*
 * node.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <cstring>
#include <functional>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

#include "node.h"
#include "string_node_object.h"

Node::Node()
    : _branch_length(0.0), _height(0.0), _number(0), _label(""),
      _comment(""), _children{}, _label_map_superdouble{} {}

Node::Node(double bl, int innumber, const string &inname)
    : _branch_length(bl), _height(0.0), _number(innumber), _label(inname),
      _comment(""), _children{}, _label_map_superdouble{} {}

vector<std::shared_ptr<Node>> Node::getChildren() { return _children; }

bool Node::isExternal() const { return _children.size() < 1; }

bool Node::isInternal() const { return _children.size() > 0; }

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

string
Node::getNewick(bool branch_lengths,
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
  } else { // EXTERNAL
    if (_label.size() > 0)
      newick_oss << _label;
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

void Node::setConditionalVector(const string &name,
                                const vector<Superdouble> &v) {
  assocDoubleVector(name, v);
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

void Node::assocDoubleVector(const string &name,
                             const vector<Superdouble> &obj) {
  vector<Superdouble> tvec(obj.size());
  for (unsigned int i = 0; i < obj.size(); i++) {
    tvec[i] = obj[i];
  }
  _label_map_superdouble[name] = tvec;
}

vector<Superdouble> *Node::getDoubleVector(const string &name) {
  return &_label_map_superdouble[name];
}

void Node::deleteDoubleVector(const string &name) {
  if (_label_map_superdouble.count(name) > 0) {
    _label_map_superdouble.erase(name);
  }
}

void Node::initSegVector() { _branch_segments = vector<BranchSegment>(); }

vector<BranchSegment> &Node::getSegVector() { return _branch_segments; }

void Node::initExclDistVector() { _excluded_dists = new vector<vector<int>>(); }

vector<vector<int>> *Node::getExclDistVector() { return _excluded_dists; }

void Node::deleteExclDistVector() { delete _excluded_dists; }

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

std::shared_ptr<Node>
getMRCAWithNode(const std::shared_ptr<Node> &current,
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

std::shared_ptr<Node> getParentWithNode(std::shared_ptr<Node> current,
                                        std::shared_ptr<Node> n) {
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
