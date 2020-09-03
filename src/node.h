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
#include "vector_node_object.h"

class Node {
private:
  double _branch_length; // branch lengths
  double _height;        // could be from tip or from root
  int _number;
  string _label;
  Node *_parent;
  vector<Node *> _children;
  map<string, NodeObject *> _label_map;
  map<string, vector<Superdouble>> _label_map_superdouble;
  string _comment;
  vector<BranchSegment> *_branch_segments;
  vector<vector<int>> *_excluded_dists;

public:
  Node();
  Node(Node *parent);
  Node(double bl, int number, string name, Node *parent);

  vector<Node *> getChildren();
  bool isExternal();
  bool isInternal();
  bool isRoot();
  bool hasParent();
  void setParent(Node &p);
  int getNumber();
  void setNumber(int n);
  double getBL();
  void setBL(double bl);
  double getHeight();
  void setHeight(double he);
  bool hasChild(Node &test);
  bool addChild(Node &c);
  bool removeChild(Node &c);
  Node &getChild(int c);
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

  double getMaxHeightRecursive() const;
  double getMaxHeight() const;
  void setHeightRecursive(double height);
  void setHeightRecursive();
  void pruneNode(Node *n);

  Node *getMCRA(const std::vector<Node *>);
};

#endif /* NODE_H_ */
