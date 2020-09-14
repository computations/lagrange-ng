/*
 * tree.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <cstring>
#include <exception>
#include <iostream>

using namespace std;

#include "tree.h"

Tree::Tree()
    : _root(NULL), _nodes(vector<std::shared_ptr<Node>>()),
      _internal_nodes(vector<std::shared_ptr<Node>>()),
      _external_nodes(vector<std::shared_ptr<Node>>()), _internal_node_count(0),
      _external_node_count(0) {
  processRoot();
}

Tree::Tree(std::shared_ptr<Node> inroot)
    : _root(inroot), _nodes(vector<std::shared_ptr<Node>>()),
      _internal_nodes(vector<std::shared_ptr<Node>>()),
      _external_nodes(vector<std::shared_ptr<Node>>()), _internal_node_count(0),
      _external_node_count(0) {
  _root = inroot;
  processRoot();
}

void Tree::addExternalNode(std::shared_ptr<Node> tn) {
  _external_nodes.push_back(tn);
  _external_node_count = _external_nodes.size();
  _nodes.push_back(tn);
}

void Tree::addInternalNode(std::shared_ptr<Node> tn) {
  _internal_nodes.push_back(tn);
  _internal_node_count = _internal_nodes.size();
  _nodes.push_back(tn);
}

std::shared_ptr<Node> Tree::getExternalNode(int num) {
  return _external_nodes.at(num);
}

/*
 * could precompute this, check for run time differences
 */
std::shared_ptr<Node> Tree::getExternalNode(string &name) {
  std::shared_ptr<Node> ret = NULL;
  for (unsigned int i = 0; i < _external_nodes.size(); i++) {
    if (_external_nodes.at(i)->getName() == name)
      ret = _external_nodes.at(i);
  }
  return ret;
}

std::shared_ptr<Node> Tree::getInternalNode(int num) {
  return _internal_nodes.at(num);
}

/*
 * could precompute this, check for run time differences
 */
std::shared_ptr<Node> Tree::getInternalNode(string &name) {
  std::shared_ptr<Node> ret = NULL;
  for (unsigned int i = 0; i < _internal_nodes.size(); i++) {
    if (_internal_nodes.at(i)->getName() == name)
      ret = _internal_nodes.at(i);
  }
  return ret;
}

unsigned int Tree::getExternalNodeCount() { return _external_node_count; }

unsigned int Tree::getInternalNodeCount() { return _internal_node_count; }

std::shared_ptr<Node> Tree::getNode(int num) { return _nodes.at(num); }

unsigned int Tree::getNodeCount() { return _nodes.size(); }

std::shared_ptr<Node> Tree::getRoot() { return _root; }

void Tree::setRoot(std::shared_ptr<Node> inroot) { _root = inroot; }

void Tree::unRoot(std::shared_ptr<Node> inroot) {
  processRoot();
  if (getRoot()->getChildCount() < 3) {
    tritomyRoot(inroot);
  }
  processRoot();
}

/*
 * seems to be working but check for leaks
 */
void Tree::reRoot(std::shared_ptr<Node> outgroup) {
  if (!findNode(outgroup)) {
    throw std::runtime_error{
        "Outgroup is not part of this tree, cannot reroot"};
  }
  processRoot();
  if (getRoot()->getChildCount() < 3) {
    tritomyRoot(outgroup); // not sure if this should actually be the inroot
                           // instead of NULL
  }
  if (_root == outgroup) {
    cout << "you asked to root at the current root" << endl;
  } else {

    std::shared_ptr<Node> tempParent = getParent(outgroup);
    auto newRoot = std::make_shared<Node>(*tempParent);

    newRoot->addChild(std::shared_ptr<Node>(outgroup));
    // outgroup->setParent(newRoot);
    tempParent->removeChild(outgroup);
    tempParent->addChild(std::shared_ptr<Node>(newRoot));
    // newRoot->setParent(tempParent);
    newRoot->setBL(outgroup->getBL() / 2);
    outgroup->setBL(outgroup->getBL() / 2);
    processReRoot(newRoot);
    setRoot(newRoot);
    processRoot();
  }
}

/*
 * seems to be working now
 */
void Tree::tritomyRoot(std::shared_ptr<Node> toberoot) {
  std::shared_ptr<Node> curroot = getRoot();
  if (toberoot == NULL) {
    if (curroot->getChild(0)->isInternal()) {
      std::shared_ptr<Node> currootCH = curroot->getChild(0);
      double nbl = currootCH->getBL();
      curroot->getChild(1)->setBL(curroot->getChild(1)->getBL() + nbl);
      curroot->removeChild(currootCH);
      for (int i = 0; i < currootCH->getChildCount(); i++) {
        curroot->addChild(currootCH->getChild(i));
      }
    } else {
      std::shared_ptr<Node> currootCH = curroot->getChild(1);
      double nbl = currootCH->getBL();
      curroot->getChild(0)->setBL(curroot->getChild(0)->getBL() + nbl);
      curroot->removeChild(currootCH);
      for (int i = 0; i < currootCH->getChildCount(); i++) {
        curroot->addChild(currootCH->getChild(i));
      }
    }
  } else {
    if (curroot->getChild(1) == toberoot) {
      std::shared_ptr<Node> currootCH = curroot->getChild(0);
      double nbl = currootCH->getBL();
      curroot->getChild(1)->setBL(curroot->getChild(1)->getBL() + nbl);
      curroot->removeChild(currootCH);
      for (int i = 0; i < currootCH->getChildCount(); i++) {
        curroot->addChild(currootCH->getChild(i));
      }
    } else {
      std::shared_ptr<Node> currootCH = curroot->getChild(1);
      double nbl = currootCH->getBL();
      curroot->getChild(0)->setBL(curroot->getChild(0)->getBL() + nbl);
      curroot->removeChild(currootCH);
      for (int i = 0; i < currootCH->getChildCount(); i++) {
        curroot->addChild(currootCH->getChild(i));
      }
    }
  }
}

std::shared_ptr<Node> Tree::getMRCA(vector<string> outgroup) {
  if (outgroup.size() == 1) {
    return getExternalNode(outgroup[0]);
  }
  std::shared_ptr<Node> mcra;
  std::vector<std::shared_ptr<Node>> outgroup_nodes;
  outgroup_nodes.reserve(outgroup.size());
  for (auto &o : outgroup) {
    outgroup_nodes.push_back(getExternalNode(o));
  }
  return getMRCA(outgroup_nodes);
}

std::shared_ptr<Node> Tree::getMRCA(vector<std::shared_ptr<Node>> outgroup) {
  return getMRCAWithNode(_root, outgroup);
}

void Tree::setHeightFromRootToNodes() {
  setHeightFromRootToNode(_root, _root->getBL());
}

void Tree::setHeightFromRootToNode(std::shared_ptr<Node> inNode,
                                   double newHeight) {
  if (inNode != _root) {
    newHeight += inNode->getBL();
    inNode->setHeight(newHeight);
  } else {
    inNode->setHeight(newHeight);
  }
  for (int i = 0; i < inNode->getChildCount(); i++) {
    setHeightFromRootToNode(inNode->getChild(i), newHeight);
  }
}

double Tree::getLongestPathRootToTip() const {
  return _root->getMaxHeightRecursive();
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
  if (_root == nullptr)
    return;
  postOrderProcessRoot(_root);
}

void Tree::processReRoot(std::shared_ptr<Node> node) {
  if (node != _root || node->isExternal()) {
    return;
  }
  if (getParent(node) != nullptr) {
    processReRoot(getParent(node));
  }
  // Exchange branch label, length et cetera
  exchangeInfo(getParent(node), node);
  // Rearrange topology
  std::shared_ptr<Node> parent = getParent(node);
  node->addChild(parent);
  parent->removeChild(node);
}

void Tree::exchangeInfo(std::shared_ptr<Node> node1,
                        std::shared_ptr<Node> node2) {
  string swaps;
  double swapd;
  swaps = node1->getName();
  node1->setName(node2->getName());
  node2->setName(swaps);

  swapd = node1->getBL();
  node1->setBL(node2->getBL());
  node2->setBL(swapd);
}

void Tree::postOrderProcessRoot(std::shared_ptr<Node> node) {
  if (node == nullptr)
    return;
  if (node->getChildCount() > 0) {
    for (int i = 0; i < node->getChildCount(); i++) {
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

void Tree::pruneExternalNode(std::shared_ptr<Node> node) {
  if (node->isInternal()) {
    return;
  }
  _root->pruneNode(node);
  processRoot();
}

/*
 * end private
 */

Tree::~Tree() { _root.reset(); }

bool Tree::findNode(std::shared_ptr<Node> n) { return _root->findNode(n); }

std::shared_ptr<Node> Tree::getParent(std::shared_ptr<Node> n) const {
  return getParentWithNode(_root, n);
}
