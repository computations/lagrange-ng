/*
 * PhyloTree.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef PHYLOTREE_H_
#define PHYLOTREE_H_

#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "AncSplit.h"
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
                       unordered_map<lagrange_dist_t, vector<AncSplit>> &ans,
                       unordered_map<int, string> &areanamemaprev,
                       const std::unordered_map<int, lagrange_dist_t> &distmap,
                       size_t area_count);
  void summarizeAncState(
      std::shared_ptr<Node> node, vector<Superdouble> &ans,
      unordered_map<int, string> &areanamemaprev,
      const std::unordered_map<int, lagrange_dist_t> &distmap, size_t areasize,
      nlohmann::json &output_json);
  string get_string_from_dist_int(
      int dist, const unordered_map<int, string> &areanamemaprev,
      const unordered_map<int, lagrange_dist_t> &distmap, size_t area_count);
};

#endif /* PHYLOTREE_H_ */
