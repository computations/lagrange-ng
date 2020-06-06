/*
 * PhyloTree.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 */

#ifndef PHYLOTREE_H_
#define PHYLOTREE_H_

#include "AncSplit.h"
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
#include "node.h"
#include "superdouble.h"
#include "tree.h"

#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif

class BioGeoTreeTools {
public:
  Tree *getTreeFromString(string treestring);
  vector<Node *> getAncestors(Tree &tree, Node &node);

  void summarizeSplits(Node *node,
                       unordered_map<vector<int>, vector<AncSplit>> &ans,
                       unordered_map<int, string> &areanamemaprev,
                       RateModel *rm);
  void summarizeAncState(Node *node, vector<Superdouble> &ans,
                         unordered_map<int, string> &areanamemaprev,
                         RateModel *rm);
  string get_string_from_dist_int(int dist,
                                  unordered_map<int, string> &areanamemaprev,
                                  RateModel *rm);
};

#endif /* PHYLOTREE_H_ */
