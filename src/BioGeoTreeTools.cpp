/*
 * PhyloTree.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */

#include "BioGeoTreeTools.h"
#include "RateMatrixUtils.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unordered_map>
using namespace std;

#include "string_node_object.h"
#include "tree_reader.h"
#include "vector_node_object.h"

std::shared_ptr<Tree> BioGeoTreeTools::getTreeFromString(string treestring) {
  TreeReader tr;
  return tr.readTree(treestring);
}

void BioGeoTreeTools::summarizeSplits(
    std::shared_ptr<Node> node,
    unordered_map<vector<int>, vector<AncSplit>> &ans,
    unordered_map<int, string> &areanamemaprev,
    const unordered_map<int, vector<int>> &distmap) {
  Superdouble best(0);
  Superdouble sum(0);
  vector<pair<Superdouble, string>> printstring;
  int areasize = (*ans.begin()).first.size();
  vector<int> bestldist;
  vector<int> bestrdist;
  unordered_map<vector<int>, vector<AncSplit>>::iterator it;
  bool first = true;
  for (it = ans.begin(); it != ans.end(); it++) {
    vector<int> dis = (*it).first;
    vector<AncSplit> tans = (*it).second;
    for (unsigned int i = 0; i < tans.size(); i++) {
      if (first == true) {
        first = false;
        best = tans[i].getLikelihood();
        bestldist = distmap.at(tans[i].ldescdistint); // tans[i].getLDescDist();
        bestrdist = distmap.at(tans[i].rdescdistint); // tans[i].getRDescDist();
      } else if (tans[i].getLikelihood() > best) {
        best = tans[i].getLikelihood();
        bestldist = distmap.at(tans[i].ldescdistint); // tans[i].getLDescDist();
        bestrdist = distmap.at(tans[i].rdescdistint); // tans[i].getRDescDist();
      }
      // cout << -log(tans[i].getLikelihood()) << endl;
      sum += tans[i].getLikelihood();
    }
  }
  Superdouble test2(2);
  for (it = ans.begin(); it != ans.end(); it++) {
    vector<int> dis = (*it).first;
    vector<AncSplit> tans = (*it).second;
    for (unsigned int i = 0; i < tans.size(); i++) {
      if ((best.getLn() - (tans[i].getLikelihood().getLn())) < test2) {
        string tdisstring = "";
        int count1 = 0;
        for (int m = 0; m < areasize;
             // tans[i].getLDescDist().size();
             m++) {
          // if(tans[i].getLDescDist()[m] == 1){
          if (distmap.at(tans[i].ldescdistint)[m] == 1) {
            tdisstring += areanamemaprev[m];
            count1 += 1;
            if (count1 <
                calculate_vector_int_sum(distmap.at(tans[i].ldescdistint))) {
              tdisstring += "_";
            }
          }
        }
        tdisstring += "|";
        count1 = 0;
        for (int m = 0; m < areasize;
             // tans[i].getRDescDist().size();
             m++) {
          if (distmap.at(tans[i].rdescdistint)[m] == 1) {
            tdisstring += areanamemaprev[m];
            count1 += 1;
            if (count1 <
                calculate_vector_int_sum(distmap.at(tans[i].rdescdistint))) {
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
  StringNodeObject disstring = "";
  int count = 0;
  for (unsigned int m = 0; m < bestldist.size(); m++) {
    if (bestldist[m] == 1) {
      disstring += areanamemaprev[m];
      count += 1;
      // if(count < calculate_vector_int_sum(&bestldist))
      if (count < accumulate(bestldist.begin(), bestldist.end(), 0))
        disstring += "_";
    }
  }
  disstring += "|";
  count = 0;
  for (unsigned int m = 0; m < bestrdist.size(); m++) {
    if (bestrdist[m] == 1) {
      disstring += areanamemaprev[m];
      count += 1;
      // if(count < calculate_vector_int_sum(&bestrdist))
      if (count < accumulate(bestrdist.begin(), bestrdist.end(), 0))
        disstring += "_";
    }
  }
  string spl = "split";
  node->assocObject(spl, disstring);
  // cout << -log(best) << " "<< best/sum << endl;
}

void BioGeoTreeTools::summarizeAncState(
    std::shared_ptr<Node> node, vector<Superdouble> &ans,
    unordered_map<int, string> &areanamemaprev,
    const unordered_map<int, vector<int>> &distmap, int areasize) {
  // Superdouble best(ans[1]);
  Superdouble best(ans[1]); // use ans[1] because ans[0] is just 0
  Superdouble sum(0);
  map<Superdouble, string> printstring;
  vector<int> bestancdist;
  Superdouble zero(0);
  for (unsigned int i = 1; i < ans.size();
       i++) { // 1 because 0 is just the 0 all extinct one
    if (ans[i] >= best && ans[i] != zero) { // added != 0, need to test
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
      int count1 = 0;
      for (int m = 0; m < areasize; m++) {
        if (distmap.at(i)[m] == 1) {
          tdisstring += areanamemaprev[m];
          count1 += 1;
          if (count1 < calculate_vector_int_sum(distmap.at(i))) {
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
  StringNodeObject disstring = "";
  int count = 0;
  for (unsigned int m = 0; m < bestancdist.size(); m++) {
    if (bestancdist[m] == 1) {
      disstring += areanamemaprev[m];
      count += 1;
      // if(count < calculate_vector_int_sum(&bestldist))
      if (count < accumulate(bestancdist.begin(), bestancdist.end(), 0))
        disstring += "_";
    }
  }
  string spl = "state";
  node->assocObject(spl, disstring);
  // cout << -log(best) << " "<< best/sum << endl;
}

string BioGeoTreeTools::get_string_from_dist_int(
    int dist, unordered_map<int, string> &areanamemaprev,
    const unordered_map<int, vector<int>> &distmap) {
  vector<int> bestancdist = distmap.at(dist);

  StringNodeObject disstring = "";
  int count = 0;
  for (unsigned int m = 0; m < bestancdist.size(); m++) {
    if (bestancdist[m] == 1) {
      disstring += areanamemaprev[m];
      count += 1;
      // if(count < calculate_vector_int_sum(&bestldist))
      if (count < accumulate(bestancdist.begin(), bestancdist.end(), 0))
        disstring += "_";
    }
  }
  return disstring;
}
