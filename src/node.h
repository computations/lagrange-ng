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
  vector<std::shared_ptr<Node>> _children;
  map<string, NodeObject *> _label_map;
  map<string, vector<Superdouble>> _label_map_superdouble;
  string _comment;
  vector<BranchSegment> *_branch_segments;
  vector<vector<int>> *_excluded_dists;

public:
  Node();
  Node(double bl, int number, string name);

  vector<std::shared_ptr<Node>> getChildren();
  bool isExternal();
  bool isInternal();
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
