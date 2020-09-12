/*
 * node.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef NODE_H_
#define NODE_H_

#include <map>
#include <string>
#include <vector>
using namespace std;

#include "BranchSegment.h"
#include "node_object.h"
#include "superdouble.h"

class Node {
private:
  double _branch_length; // branch lengths
  double _height;        // could be from tip or from root
  int _number;
  string _label;
  Node *_parent; // This needs to be a bare pointer, for now, so that we can set
                 // when we add a child to a node, we don't end up with 2 shared
                 // pointers pointing to the same memory. There might be a way
                 // to fix this eventually, but for now, we just say that if all
                 // the children are shared pointers, then the only way for the
                 // child to have a parent is for the parent to still be
                 // allocated. This isn't true, strictly speaking, but for now
                 // it is good enough.
  vector<std::shared_ptr<Node>> _children;
  map<string, NodeObject *> _label_map;
  map<string, vector<Superdouble>> _label_map_superdouble;
  string _comment;
  vector<BranchSegment> *_branch_segments;
  vector<vector<int>> *_excluded_dists;

public:
  Node();
  Node(Node *parent);
  Node(double bl, int number, string name, Node *parent);

  vector<std::shared_ptr<Node>> getChildren();
  bool isExternal();
  bool isInternal();
  bool isRoot();
  bool hasParent();
  void setParent(Node *p);
  int getNumber();
  void setNumber(int n);
  double getBL();
  void setBL(double bl);
  double getHeight();
  void setHeight(double he);
  bool hasChild(std::shared_ptr<Node> test);
  bool addChild(std::shared_ptr<Node> c);
  bool removeChild(std::shared_ptr<Node> c);
  std::shared_ptr<Node> getChild(int c);
  string getName();
  string getComment();
  void setName(string s);
  void setComment(string s);
  string getNewick(bool bl);
  string getNewickOBL(string obj);
  string getNewick(bool bl, string obj);
  Node *getParent();
  int getChildCount();
  void assocObject(string name, NodeObject &obj);
  void assocDoubleVector(string name, vector<Superdouble> &obj);
  vector<Superdouble> *getDoubleVector(string name);
  void deleteDoubleVector(string name);
  void initSegVector();
  vector<BranchSegment> *getSegVector();
  void deleteSegVector();
  void initExclDistVector();
  vector<vector<int>> *getExclDistVector();
  void deleteExclDistVector();
  NodeObject *getObject(string name);
  bool findNode(std::shared_ptr<Node> n);

  double getMaxHeightRecursive() const;
  double getMaxHeight() const;
  void setHeightRecursive(double height);
  void setHeightRecursive();
  void pruneNode(std::shared_ptr<Node> n);

  friend std::shared_ptr<Node> getParentWithNode(std::shared_ptr<Node> current,
                                                 std::shared_ptr<Node> n);
  friend std::shared_ptr<Node>
  getMRCAWithNode(const std::shared_ptr<Node> &current,
                  const std::vector<std::shared_ptr<Node>> &nodes);
};
std::shared_ptr<Node>
getMRCAWithNode(const std::shared_ptr<Node> &current,
                const std::vector<std::shared_ptr<Node>> &nodes);

#endif /* NODE_H_ */
