/*
 * PhyloTree.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unordered_map>

#include "BioGeoTreeTools.h"
#include "RateMatrixUtils.h"
using namespace std;

#include "Utils.h"
#include "string_node_object.h"
#include "tree_reader.h"
#include "vector_node_object.h"

std::shared_ptr<Tree> BioGeoTreeTools::getTreeFromString(string treestring) {
  TreeReader tr;
  return tr.readTree(treestring);
}

void BioGeoTreeTools::summarizeSplits(
    std::shared_ptr<Node> node,
    unordered_map<lagrange_dist_t, vector<AncSplit>> &ans,
    unordered_map<int, string> &areanamemaprev,
    const unordered_map<int, lagrange_dist_t> &distmap, size_t area_count) {
  Superdouble sum(0);
  vector<pair<Superdouble, string>> printstring;

  const auto &first = ans.begin()->second[0];
  Superdouble best(first.getLikelihood());
  lagrange_dist_t bestldist = distmap.at(first.ldescdistint);
  lagrange_dist_t bestrdist = distmap.at(first.rdescdistint);

  for (auto it = ans.begin(); it != ans.end(); it++) {
    vector<AncSplit> tans = (*it).second;
    for (unsigned int i = 0; i < tans.size(); i++) {
      if (tans[i].getLikelihood() > best) {
        best = tans[i].getLikelihood();
        bestldist = distmap.at(tans[i].ldescdistint);
        bestrdist = distmap.at(tans[i].rdescdistint);
      }
      sum += tans[i].getLikelihood();
    }
  }
  Superdouble test2(2);
  for (auto it = ans.begin(); it != ans.end(); it++) {
    vector<AncSplit> tans = (*it).second;
    for (unsigned int i = 0; i < tans.size(); i++) {
      if ((best.getLn() - (tans[i].getLikelihood().getLn())) < test2) {
        string tdisstring = "";
        size_t count1 = 0;
        for (size_t m = 0; m < area_count; m++) {
          if (lagrange_bextr(distmap.at(tans[i].ldescdistint), m) == 1) {
            tdisstring += areanamemaprev[m];
            count1 += 1;
            if (count1 < lagrange_popcount(distmap.at(tans[i].ldescdistint))) {
              tdisstring += "_";
            }
          }
        }
        tdisstring += "|";
        count1 = 0;
        for (size_t m = 0; m < area_count; m++) {
          if (lagrange_bextr(distmap.at(tans[i].rdescdistint), m) == 1) {
            tdisstring += areanamemaprev[m];
            count1 += 1;
            if (count1 < lagrange_popcount(distmap.at(tans[i].rdescdistint))) {
              tdisstring += "_";
            }
          }
        }
        printstring.push_back(
            std::make_pair(tans[i].getLikelihood(), tdisstring));
      }
    }
  }
  sort(printstring.begin(), printstring.end(),
       [](auto a, auto b) { return a < b; });
  Superdouble none(-1);
  vector<pair<Superdouble, string>>::reverse_iterator pit;
  for (pit = printstring.rbegin(); pit != printstring.rend(); pit++) {
    Superdouble lnl(((*pit).first));
    cout << "\t" << (*pit).second << "\t" << double(lnl / sum) << "\t("
         << double(none * lnl.getLn()) << ")" << endl;
  }
  string disstring = "";
  size_t count = 0;
  for (unsigned int m = 0; m < area_count; m++) {
    if (lagrange_bextr(bestldist, m) == 1) {
      disstring += areanamemaprev[m];
      count += 1;
      if (count < lagrange_popcount(bestldist)) {
        disstring += "_";
      }
    }
  }
  disstring += "|";
  count = 0;
  for (unsigned int m = 0; m < area_count; m++) {
    if (lagrange_bextr(bestrdist, m) == 1) {
      disstring += areanamemaprev[m];
      count += 1;
      if (count < lagrange_popcount(bestrdist)) {
        disstring += "_";
      }
    }
  }
  node->setSplitString(disstring);
}

void BioGeoTreeTools::summarizeAncState(
    std::shared_ptr<Node> node, vector<Superdouble> &ans,
    unordered_map<int, string> &areanamemaprev,
    const unordered_map<int, lagrange_dist_t> &distmap, size_t areasize) {
  // Superdouble best(ans[1]);
  Superdouble best(ans[1]);  // use ans[1] because ans[0] is just 0
  Superdouble sum(0);
  map<Superdouble, string> printstring;
  lagrange_dist_t bestancdist = 0;
  Superdouble zero(0);
  // start from 1 because 0 is just the 0 all extinct one
  for (unsigned int i = 1; i < ans.size(); i++) {
    if (ans[i] >= best && ans[i] != zero) {  // added != 0, need to test
      best = ans[i];
      bestancdist = distmap.at(i);
    }
    sum += ans[i];
  }
  Superdouble none(-1);
  Superdouble test2(2);
  for (unsigned int i = 0; i < ans.size(); i++) {
    if (((best.getLn()) - (ans[i].getLn())) < test2) {
      string tdisstring = "";
      size_t count1 = 0;
      for (size_t m = 0; m < areasize; m++) {
        if (lagrange_bextr(distmap.at(i), m) == 1) {
          tdisstring += areanamemaprev[m];
          count1 += 1;
          if (count1 < lagrange_popcount(distmap.at(i))) {
            tdisstring += "_";
          }
        }
      }
      printstring[ans[i]] = tdisstring;
    }
  }

  map<Superdouble, string>::reverse_iterator pit;
  for (pit = printstring.rbegin(); pit != printstring.rend(); pit++) {
    Superdouble lnl(((*pit).first));
    // cout << lnl << endl;
    cout << "\t" << (*pit).second << "\t" << double(lnl / sum) << "\t("
         << double(none * lnl.getLn()) << ")" << endl;
  }
  string disstring = "";
  size_t count = 0;
  for (unsigned int m = 0; m < areasize; m++) {
    if (lagrange_bextr(bestancdist, m) == 1) {
      disstring += areanamemaprev[m];
      count += 1;
      // if(count < calculate_vector_int_sum(&bestldist))
      if (count < lagrange_popcount(bestancdist)) {
        disstring += "_";
      }
    }
  }
  node->setStateString(disstring);
  // cout << -log(best) << " "<< best/sum << endl;
}

string BioGeoTreeTools::get_string_from_dist_int(
    int dist, const unordered_map<int, string> &areanamemaprev,
    const unordered_map<int, lagrange_dist_t> &distmap, size_t area_count) {
  lagrange_dist_t bestancdist = distmap.at(dist);

  string disstring = "";
  size_t count = 0;
  for (unsigned int m = 0; m < area_count; m++) {
    if (lagrange_bextr(bestancdist, m) == 1) {
      disstring += areanamemaprev.at(m);
      count += 1;
      // if(count < calculate_vector_int_sum(&bestldist))
      if (count < lagrange_popcount(bestancdist)) {
        disstring += "_";
      }
    }
  }
  return disstring;
}
