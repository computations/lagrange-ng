/*
 * tree.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef TREE_H_
#define TREE_H_

#include <string>
#include <vector>

using namespace std;

#include "node.h"

class Tree {
private:
  Node *_root;
  vector<Node *> _nodes;
  vector<Node *> _internal_nodes;
  vector<Node *> _external_nodes;
  int _internal_node_count;
  int _external_node_count;

  void processReRoot(Node *node);
  void exchangeInfo(Node *node1, Node *node2);
  void postOrderProcessRoot(Node *node);
  Node *getMRCATraverse(Node *curn1, Node *curn2);
  void setHeightFromRootToNode(Node &inNode, double newHeight);
  double getGreatestDistance(Node *inNode);

public:
  Tree();
  Tree(Node *root);

  void addExternalNode(Node *tn);
  void addInternalNode(Node *tn);
  void pruneExternalNode(Node *node);
  Node *getExternalNode(int num);
  Node *getExternalNode(string &name);
  Node *getInternalNode(int num);
  Node *getInternalNode(string &name);
  Node *getNode(int num);
  int getNodeCount();
  int getExternalNodeCount();
  int getInternalNodeCount();
  Node *getRoot();
  void setRoot(Node *inroot);
  void unRoot(Node &inroot);
  void reRoot(Node *inroot);
  void tritomyRoot(Node *toberoot);
  Node *getMRCA(vector<string> innodes);
  Node *getMRCA(vector<Node *> innodes);
  void processRoot();
  double getLongestPathRootToTip() const;

  void setHeightFromRootToNodes();
  void setHeightFromTipToNodes();

  ~Tree();
};

#endif /* TREE_H_ */
