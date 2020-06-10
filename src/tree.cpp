/*
 * tree.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <exception>

using namespace std;

#include "node.h"
#include "tree.h"

Tree::Tree()
    : _root(NULL), _nodes(vector<Node *>()), _internal_nodes(vector<Node *>()),
      _external_nodes(vector<Node *>()), _internal_node_count(0),
      _external_node_count(0) {
  processRoot();
}

Tree::Tree(Node *inroot)
    : _root(inroot), _nodes(vector<Node *>()), _internal_nodes(vector<Node *>()),
      _external_nodes(vector<Node *>()), _internal_node_count(0),
      _external_node_count(0) {
  _root = inroot;
  processRoot();
}

void Tree::addExternalNode(Node *tn) {
  _external_nodes.push_back(tn);
  _external_node_count = _external_nodes.size();
  _nodes.push_back(tn);
}

void Tree::addInternalNode(Node *tn) {
  _internal_nodes.push_back(tn);
  _internal_node_count = _internal_nodes.size();
  _nodes.push_back(tn);
}

Node *Tree::getExternalNode(int num) { return _external_nodes.at(num); }

/*
 * could precompute this, check for run time differences
 */
Node *Tree::getExternalNode(string &name) {
  Node *ret = NULL;
  for (unsigned int i = 0; i < _external_nodes.size(); i++) {
    if (_external_nodes.at(i)->getName() == name)
      ret = _external_nodes.at(i);
  }
  return ret;
}

Node *Tree::getInternalNode(int num) { return _internal_nodes.at(num); }

/*
 * could precompute this, check for run time differences
 */
Node *Tree::getInternalNode(string &name) {
  Node *ret = NULL;
  for (unsigned int i = 0; i < _internal_nodes.size(); i++) {
    if (_internal_nodes.at(i)->getName() == name)
      ret = _internal_nodes.at(i);
  }
  return ret;
}

int Tree::getExternalNodeCount() { return _external_node_count; }

int Tree::getInternalNodeCount() { return _internal_node_count; }

Node *Tree::getNode(int num) { return _nodes.at(num); }

int Tree::getNodeCount() { return _nodes.size(); }

Node *Tree::getRoot() { return _root; }

void Tree::setRoot(Node *inroot) { _root = inroot; }

void Tree::unRoot(Node &inroot) {
  processRoot();
  if (getRoot()->getChildCount() < 3) {
    tritomyRoot(&inroot);
  }
  processRoot();
}

/*
 * seems to be working but check for leaks
 */
void Tree::reRoot(Node *inroot) {
  processRoot();
  if (getRoot()->getChildCount() < 3) {
    tritomyRoot(inroot); // not sure if this should actually be the inroot
                         // instead of NULL
  }
  // cout << this->root->getNewick(false) << endl;
  if (_root == inroot) {
    cout << "you asked to root at the current root" << endl;
  } else {
    Node *tempParent = inroot->getParent();
    Node *newRoot = new Node(tempParent);
    newRoot->addChild(*inroot);
    inroot->setParent(*newRoot);
    tempParent->removeChild(*inroot);
    tempParent->addChild(*newRoot);
    newRoot->setParent(*tempParent);
    newRoot->setBL(inroot->getBL() / 2);
    inroot->setBL(inroot->getBL() / 2);
    processReRoot(newRoot);
    setRoot(newRoot);
    processRoot();
  }
}

/*
 * seems to be working now
 */
void Tree::tritomyRoot(Node *toberoot) {
  Node *curroot = getRoot();
  if (toberoot == NULL) {
    if (curroot->getChild(0).isInternal()) {
      Node *currootCH = &curroot->getChild(0);
      double nbl = currootCH->getBL();
      curroot->getChild(1).setBL(curroot->getChild(1).getBL() + nbl);
      curroot->removeChild(*currootCH);
      for (int i = 0; i < currootCH->getChildCount(); i++) {
        curroot->addChild(currootCH->getChild(i));
      }
    } else {
      Node *currootCH = &curroot->getChild(1);
      double nbl = currootCH->getBL();
      curroot->getChild(0).setBL(curroot->getChild(0).getBL() + nbl);
      curroot->removeChild(*currootCH);
      for (int i = 0; i < currootCH->getChildCount(); i++) {
        curroot->addChild(currootCH->getChild(i));
      }
    }
  } else {
    if (&curroot->getChild(1) == toberoot) {
      Node *currootCH = &curroot->getChild(0);
      double nbl = currootCH->getBL();
      curroot->getChild(1).setBL(curroot->getChild(1).getBL() + nbl);
      curroot->removeChild(*currootCH);
      for (int i = 0; i < currootCH->getChildCount(); i++) {
        curroot->addChild(currootCH->getChild(i));
      }
    } else {
      Node *currootCH = &curroot->getChild(1);
      double nbl = currootCH->getBL();
      curroot->getChild(0).setBL(curroot->getChild(0).getBL() + nbl);
      curroot->removeChild(*currootCH);
      for (int i = 0; i < currootCH->getChildCount(); i++) {
        curroot->addChild(currootCH->getChild(i));
      }
    }
  }
}

Node *Tree::getMRCA(vector<string> innodes) {
  Node *mrca = NULL;
  if (innodes.size() == 1)
    return getExternalNode(innodes[0]);
  else {
    vector<string> outgroup;
    for (unsigned int i = 0; i < innodes.size(); i++) {
      outgroup.push_back(innodes.at(i));
    }
    Node *cur1 = getExternalNode(outgroup.at(0));
    outgroup.erase(outgroup.begin());
    Node *cur2 = NULL;
    Node *tempmrca = NULL;
    while (outgroup.size() > 0) {
      cur2 = getExternalNode(outgroup.at(0));
      outgroup.erase(outgroup.begin());
      tempmrca = getMRCATraverse(cur1, cur2);
      cur1 = tempmrca;
    }
    mrca = cur1;
  }
  return mrca;
}

Node *Tree::getMRCA(vector<Node *> innodes) {
  Node *mrca = NULL;
  if (innodes.size() == 1)
    return innodes[0];
  else {
    Node *cur1 = innodes.at(0);
    innodes.erase(innodes.begin());
    Node *cur2 = NULL;
    Node *tempmrca = NULL;
    while (innodes.size() > 0) {
      cur2 = innodes.at(0);
      innodes.erase(innodes.begin());
      tempmrca = getMRCATraverse(cur1, cur2);
      cur1 = tempmrca;
    }
    mrca = cur1;
  }
  return mrca;
}

void Tree::setHeightFromRootToNodes() {
  setHeightFromRootToNode(*_root, _root->getBL());
}

void Tree::setHeightFromRootToNode(Node &inNode, double newHeight) {
  if (inNode.isRoot() == false) {
    newHeight += inNode.getBL();
    inNode.setHeight(newHeight);
  } else {
    inNode.setHeight(newHeight);
  }
  for (int i = 0; i < inNode.getChildCount(); i++) {
    setHeightFromRootToNode(inNode.getChild(i), newHeight);
  }
}

/*
 * only makes sense for ultrametric trees
 */
void Tree::setHeightFromTipToNodes() {
  for (int i = 0; i < _external_node_count; i++) {
    double curh = 0.0;
    Node *cur = getExternalNode(i);
    cur->setHeight(curh);
    while (cur->getParent() != NULL) {
      curh += cur->getBL();
      cur = cur->getParent();
      if (cur->getHeight() < curh)
        cur->setHeight(curh);
    }
  }
}

/*
 * private
 */
void Tree::processRoot() {
  _nodes.clear();
  _internal_nodes.clear();
  _external_nodes.clear();
  _internal_node_count = 0;
  _external_node_count = 0;
  if (_root == NULL)
    return;
  postOrderProcessRoot(_root);
}

void Tree::processReRoot(Node *node) {
  if (node->isRoot() || node->isExternal()) {
    return;
  }
  if (node->getParent() != NULL) {
    processReRoot(node->getParent());
  }
  // Exchange branch label, length et cetera
  exchangeInfo(node->getParent(), node);
  // Rearrange topology
  Node *parent = node->getParent();
  node->addChild(*parent);
  parent->removeChild(*node);
  parent->setParent(*node);
}

void Tree::exchangeInfo(Node *node1, Node *node2) {
  string swaps;
  double swapd;
  swaps = node1->getName();
  node1->setName(node2->getName());
  node2->setName(swaps);

  swapd = node1->getBL();
  node1->setBL(node2->getBL());
  node2->setBL(swapd);
}

void Tree::postOrderProcessRoot(Node *node) {
  if (node == NULL)
    return;
  if (node->getChildCount() > 0) {
    for (int i = 0; i < node->getChildCount(); i++) {
      postOrderProcessRoot(&node->getChild(i));
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

void Tree::pruneExternalNode(Node *node) {
  if (node->isInternal()) {
    return;
  }
  /*
   * how this works
   *
   * get the parent = parent
   * get the parent of the parent = mparent
   * remove parent from mparent
   * add !node from parent to mparent
   *
   * doesn't yet take care if node.parent == root
   * or polytomy
   */
  double bl = 0;
  Node *parent = node->getParent();
  Node *other = nullptr;
  for (int i = 0; i < parent->getChildCount(); i++) {
    if (&parent->getChild(i) != node) {
      other = &parent->getChild(i);
    }
  }
  if(other == nullptr){
    throw std::runtime_error{"Failed to find the other node when pruning"};
  }
  bl = other->getBL() + parent->getBL();
  Node *mparent = parent->getParent();
  if (mparent != NULL) {
    mparent->addChild(*other);
    other->setBL(bl);
    for (int i = 0; i < mparent->getChildCount(); i++) {
      if (&mparent->getChild(i) == parent) {
        mparent->removeChild(*parent);
        break;
      }
    }
  }
  delete node;
  processRoot();
}

Node *Tree::getMRCATraverse(Node *curn1, Node *curn2) {
  Node *mrca = NULL;
  // get path to root for first node
  vector<Node *> path1;
  Node *parent = curn1;
  path1.push_back(parent);
  while (parent != NULL) {
    path1.push_back(parent);
    if (parent->getParent() != NULL)
      parent = parent->getParent();
    else
      break;
  }
  // find first match between this node and the first one
  parent = curn2;
  bool x = true;
  while (x == true) {
    for (unsigned int i = 0; i < path1.size(); i++) {
      if (parent == path1.at(i)) {
        mrca = parent;
        x = false;
        break;
      }
    }
    parent = parent->getParent();
  }
  return mrca;
}

/*
 * end private
 */

Tree::~Tree() {
  for (int i = 0; i < _internal_node_count; i++) {
    delete getInternalNode(i);
  }
  for (int i = 0; i < _external_node_count; i++) {
    delete getExternalNode(i);
  }
}
