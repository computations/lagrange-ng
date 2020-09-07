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
  std::shared_ptr<Node> _root;
  vector<std::shared_ptr<Node>> _nodes;
  vector<std::shared_ptr<Node>> _internal_nodes;
  vector<std::shared_ptr<Node>> _external_nodes;
  unsigned int _internal_node_count;
  unsigned int _external_node_count;

  void processReRoot(std::shared_ptr<Node> node);
  void exchangeInfo(std::shared_ptr<Node> node1, std::shared_ptr<Node> node2);
  void postOrderProcessRoot(std::shared_ptr<Node> node);
  void setHeightFromRootToNode(std::shared_ptr<Node> inNode, double newHeight);
  double getGreatestDistance(std::shared_ptr<Node> inNode);

  bool findNode(std::shared_ptr<Node> n);
  std::shared_ptr<Node> getParent(std::shared_ptr<Node> n);

public:
  Tree();
  Tree(std::shared_ptr<Node> root);

  void addExternalNode(std::shared_ptr<Node> tn);
  void addInternalNode(std::shared_ptr<Node> tn);
  void pruneExternalNode(std::shared_ptr<Node> node);
  std::shared_ptr<Node> getExternalNode(int num);
  std::shared_ptr<Node> getExternalNode(string &name);
  std::shared_ptr<Node> getInternalNode(int num);
  std::shared_ptr<Node> getInternalNode(string &name);
  std::shared_ptr<Node> getNode(int num);
  unsigned int getNodeCount();
  unsigned int getExternalNodeCount();
  unsigned int getInternalNodeCount();
  std::shared_ptr<Node> getRoot();
  void setRoot(std::shared_ptr<Node> inroot);
  void unRoot(std::shared_ptr<Node> inroot);
  void reRoot(std::shared_ptr<Node> inroot);
  void tritomyRoot(std::shared_ptr<Node> toberoot);
  std::shared_ptr<Node> getMRCA(vector<string> innodes);
  std::shared_ptr<Node> getMRCA(vector<std::shared_ptr<Node>> innodes);
  void processRoot();
  double getLongestPathRootToTip() const;

  void setHeightFromRootToNodes();
  void setHeightFromTipToNodes();

  ~Tree();
};

#endif /* TREE_H_ */
