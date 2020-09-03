/*
 * node.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#include "BranchSegment.h"
#include "node.h"
#include "node_object.h"
#include "string_node_object.h"
#include "superdouble.h"
#include "vector_node_object.h"

Node::Node()
    : _branch_length(0.0), _height(0.0), _number(0), _label(""), _parent(NULL),
      _children(vector<Node *>()), _label_map(map<string, NodeObject *>()),
      _label_map_superdouble(map<string, vector<Superdouble>>()), _comment("") {
}

Node::Node(Node *inparent)
    : _branch_length(0.0), _height(0.0), _number(0), _label(""),
      _parent(inparent), _children(vector<Node *>()),
      _label_map(map<string, NodeObject *>()),
      _label_map_superdouble(map<string, vector<Superdouble>>()), _comment("") {
}

Node::Node(double bl, int innumber, string inname, Node *inparent)
    : _branch_length(bl), _height(0.0), _number(innumber), _label(inname),
      _parent(inparent), _children(vector<Node *>()),
      _label_map(map<string, NodeObject *>()),
      _label_map_superdouble(map<string, vector<Superdouble>>()), _comment("") {
}

vector<Node *> Node::getChildren() { return _children; }

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

bool Node::isRoot() {
  if (_parent == NULL)
    return true;
  else
    return false;
}

bool Node::hasParent() {
  if (_parent == NULL)
    return false;
  else
    return true;
}

void Node::setParent(Node &p) { _parent = &p; }

int Node::getNumber() { return _number; }

void Node::setNumber(int n) { _number = n; }

double Node::getBL() { return _branch_length; }

void Node::setBL(double bl) { _branch_length = bl; }

double Node::getHeight() { return _height; }

void Node::setHeight(double he) { _height = he; }

bool Node::hasChild(Node &test) {
  bool ret = false;
  for (unsigned int i = 0; i < _children.size(); i++) {
    if (_children.at(i) == &test) {
      ret = true;
      break;
    }
  }
  return ret;
}

bool Node::addChild(Node &c) {
  if (hasChild(c) == false) {
    _children.push_back(&c);
    c.setParent(*this);
    return true;
  } else {
    return false;
  }
}

bool Node::removeChild(Node &c) {
  if (hasChild(c) == true) {
    for (unsigned int i = 0; i < _children.size(); i++) {
      if (_children.at(i) == &c) {
        _children.erase(_children.begin() + i);
        break;
      }
    }
    return true;
  } else {
    return false;
  }
}

Node &Node::getChild(int c) { return *_children.at(c); }

string Node::getName() { return _label; }

void Node::setName(string s) { _label = s; }

string Node::getComment() { return _comment; }

void Node::setComment(string s) { _comment = s; }

string Node::getNewick(bool bl) {
  string ret = "";
  for (int i = 0; i < getChildCount(); i++) {
    if (i == 0)
      ret = ret + "(";
    ret = ret + getChild(i).getNewick(bl);
    if (bl == true) {
      std::ostringstream o;
      o << getChild(i).getBL();
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
    ret = ret + getChild(i).getNewick(bl, obj);
    if (bl == true) {
      std::ostringstream o;
      o << getChild(i).getBL();
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
    ret = ret + getChild(i).getNewickOBL(obj);
    std::ostringstream o;
    double bl = (*(VectorNodeObject<double> *)getChild(i).getObject(obj))[0];
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

Node *Node::getParent() { return _parent; }

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

void Node::pruneNode(Node *n) {
  for (auto it = _children.begin(); it != _children.end(); ++it) {
    if (*it == n) {
      delete *it;
      _children.erase(it);
      if (_children.size() == 1) {
        auto child = _children[0];
        double tmp_branch_length = _branch_length + child->_branch_length;
        (*this) = child;
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

Node *Node::getMCRA(const std::vector<Node *> leaves) {
  if (_children.size() == 0) {
    for (auto &n : leaves) {
      if (n == this) {
        return this;
      }
    }
    return nullptr;
  } else {
    Node *mcra = nullptr;
    for (auto &c : _children) {
      auto tmp_mcra = c->getMCRA(leaves);
      if (tmp_mcra != nullptr) {
        // if this is our second match, then we can return with this pointer.
        if (mcra != nullptr) {
          return this;
        }
        mcra = tmp_mcra;
      }
    }
    return mcra;
  }
}
