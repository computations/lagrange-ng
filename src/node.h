/*
 * node.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef NODE_H_
#define NODE_H_

#include <functional>
#include <map>
#include <string>
#include <vector>
using namespace std;

#include "BranchSegment.h"
#include "superdouble.h"

class Node {
private:
  double _branch_length; // branch lengths
  double _height;        // could be from tip or from root
  int _number;
  string _label;
  string _comment;
  string _split_string;
  string _state_string;
  string _stoch_string;
  vector<std::shared_ptr<Node>> _children;
  map<string, vector<Superdouble>> _label_map_superdouble;
  vector<BranchSegment> _branch_segments;
  vector<vector<int>> *_excluded_dists;

public:
  Node();
  Node(double bl, int number, const string &name);

  vector<std::shared_ptr<Node>> getChildren();
  bool isExternal() const;
  bool isInternal() const;

  int getNumber() const;
  void setNumber(int n);

  double getBL();
  void setBL(double bl);

  double getHeight();
  void setHeight(double he);

  bool hasChild(std::shared_ptr<Node> test);
  bool addChild(std::shared_ptr<Node> c);
  bool removeChild(std::shared_ptr<Node> c);
  std::shared_ptr<Node> getChild(int c) const;

  string getName();
  string getComment();
  void setName(const string &s);
  void setComment(const string &s);

  string getNewick(bool bl) const;
  string getNewick(bool branch_lengths,
                   const std::function<string(const Node &)> &) const;
  string getNewickLambda(const std::function<string(const Node &)> &) const;

  int getChildCount() const;

  void setConditionalVector(const string &name, const vector<Superdouble> &v);
  void setSplitString(const string &splitstring);
  void setStateString(const string &splitstring);
  void setStochString(const string &stochstring);
  string getSplitString() const;
  string getStateString() const;
  string getStochString() const;

  void assocDoubleVector(const string &name, const vector<Superdouble> &obj);
  vector<Superdouble> *getDoubleVector(const string &name);
  void deleteDoubleVector(const string &name);

  void initSegVector();
  vector<BranchSegment> &getSegVector();

  void initExclDistVector();
  vector<vector<int>> *getExclDistVector();
  void deleteExclDistVector();

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
