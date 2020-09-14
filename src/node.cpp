/*
 * node.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <cstring>
#include <iostream>
#include <sstream>

using namespace std;

#include "node.h"
#include "string_node_object.h"

Node::Node()
    : _branch_length(0.0), _height(0.0), _number(0),
      _label(""), _children{}, _label_map{}, _label_map_superdouble{},
      _comment("") {}

Node::Node(double bl, int innumber, string inname)
    : _branch_length(bl), _height(0.0), _number(innumber),
      _label(inname), _children{}, _label_map{}, _label_map_superdouble{},
      _comment("") {}

vector<std::shared_ptr<Node>> Node::getChildren() { return _children; }

bool Node::isExternal() {
  if (_children.size() < 1)
    return true;
  else
    return false;
}

bool Node::isInternal() {
  if (_children.size() > 0)
    return true;
  else
    return false;
}

int Node::getNumber() { return _number; }

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
  } else {
    return false;
  }
}

bool Node::removeChild(std::shared_ptr<Node> c) {
  if (hasChild(c) == true) {
    for (auto it = _children.begin(); it != _children.end(); it++) {
      if (*it == c) {
        _children.erase(it);
      }
    }
    return true;
  } else {
    return false;
  }
}

std::shared_ptr<Node> Node::getChild(int c) { return _children.at(c); }

string Node::getName() { return _label; }

void Node::setName(string s) { _label = s; }

string Node::getComment() { return _comment; }

void Node::setComment(string s) { _comment = s; }

string Node::getNewick(bool bl) {
  string ret = "";
  for (int i = 0; i < getChildCount(); i++) {
    if (i == 0)
      ret = ret + "(";
    ret = ret + getChild(i)->getNewick(bl);
    if (bl == true) {
      std::ostringstream o;
      o << getChild(i)->getBL();
      ret = ret + ":" + o.str();
    }
    if (i == getChildCount() - 1)
      ret = ret + ")";
    else
      ret = ret + ",";
  }
  if (_label.size() > 0)
    ret = ret + _label;
  return ret;
}

/*
 * should be returning the stringnodeobjects as the names for internal
 * nodes
 */
string Node::getNewick(bool bl, string obj) {
  string ret = "";
  for (int i = 0; i < getChildCount(); i++) {
    if (i == 0)
      ret = ret + "(";
    ret = ret + getChild(i)->getNewick(bl, obj);
    if (bl == true) {
      std::ostringstream o;
      o << getChild(i)->getBL();
      ret = ret + ":" + o.str();
    }
    if (i == getChildCount() - 1)
      ret = ret + ")";
    else
      ret = ret + ",";
  }
  if (isInternal() == true) {
    if (obj == "number") {
      std::ostringstream o;
      o << _number;
      ret = ret + o.str();
    } else {
      if (getObject(obj) != NULL) {
        std::ostringstream o;
        o << (*((StringNodeObject *)(getObject(obj))));
        ret = ret + o.str();
      }
    }
  } else { // EXTERNAL
    if (_label.size() > 0)
      ret = ret + _label;
  }
  return ret;
}

/*
 * should return with branch lengths determined from the string obj
 */
string Node::getNewickOBL(string obj) {
  string ret = "";
  for (int i = 0; i < getChildCount(); i++) {
    if (i == 0)
      ret = ret + "(";
    ret = ret + getChild(i)->getNewickOBL(obj);
    std::ostringstream o;
    double bl = (*(VectorNodeObject<double> *)getChild(i)->getObject(obj))[0];
    o << bl;
    ret = ret + ":" + o.str();

    if (i == getChildCount() - 1)
      ret = ret + ")";
    else
      ret = ret + ",";
  }
  if (_label.size() > 0)
    ret = ret + _label;
  return ret;
}

int Node::getChildCount() { return _children.size(); }

void Node::assocObject(string name, NodeObject &obj) {
  // need to see why this doesn't work
  _label_map[name] = obj.clone();
}

void Node::assocDoubleVector(string name, vector<Superdouble> &obj) {
  if (_label_map_superdouble.count(name) > 0) {
    _label_map_superdouble.erase(name);
  }
  vector<Superdouble> tvec(obj.size());
  for (unsigned int i = 0; i < obj.size(); i++) {
    tvec[i] = obj[i];
  }
  _label_map_superdouble[name] = tvec;
}

vector<Superdouble> *Node::getDoubleVector(string name) {
  return &_label_map_superdouble[name];
}

void Node::deleteDoubleVector(string name) {
  if (_label_map_superdouble.count(name) > 0) {
    _label_map_superdouble.erase(name);
  }
}

void Node::initSegVector() { _branch_segments = new vector<BranchSegment>(); }

vector<BranchSegment> *Node::getSegVector() { return _branch_segments; }

void Node::deleteSegVector() { delete _branch_segments; }

void Node::initExclDistVector() { _excluded_dists = new vector<vector<int>>(); }

vector<vector<int>> *Node::getExclDistVector() { return _excluded_dists; }

void Node::deleteExclDistVector() { delete _excluded_dists; }

/*
 * use the string ones like this
 * StringNodeObject sno("...a node object");
 * tree.getRoot()->assocObject("test",sno);
 * cout << *((StringNodeObject*) (tree.getRoot()->getObject("test"))) << endl;
 *
 * and the vector like
 * VectorNodeObject<int> vno;
 * vno.push_back(1);vno.push_back(2);
 * tree.getRoot()->assocObject("testvno",vno);
 * cout << ((VectorNodeObject<int> *)
 * (tree.getRoot()->getObject("testvno")))->at(0) << endl;
 */

NodeObject *Node::getObject(string name) { return _label_map[name]; }

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
