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
  std::shared_ptr<Tree> getTreeFromString(string treestring);

  void summarizeSplits(std::shared_ptr<Node> node,
                       unordered_map<vector<int>, vector<AncSplit>> &ans,
                       unordered_map<int, string> &areanamemaprev,
                       std::shared_ptr<RateModel> rm);
  void summarizeAncState(std::shared_ptr<Node> node, vector<Superdouble> &ans,
                         unordered_map<int, string> &areanamemaprev,
                         std::shared_ptr<RateModel> rm);
  string get_string_from_dist_int(int dist,
                                  unordered_map<int, string> &areanamemaprev,
                                  std::shared_ptr<RateModel> rm);
};

#endif /* PHYLOTREE_H_ */
