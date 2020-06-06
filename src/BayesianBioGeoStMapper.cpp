/*
 *  BayesianBioGeo.cpp
 *  lagrange_cpp
 *
 *  Created by Stephen Smith on 1/13/10.
 *  Copyright 2010 Yale University. All rights reserved.
 *
 */

#include "BayesianBioGeoStMapper.h"

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
#include <stack>
#include <string>
#include <vector>
using namespace std;

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "BioGeoTree.h"
#include "BranchSegment.h"
#include "RateMatrixUtils.h"
#include "RateModel.h"
#include "node.h"
#include "tree.h"

namespace {
inline double MIN(const double &a, const double &b) {
  return b < a ? (b) : double(a);
}
} // namespace

BayesianBioGeoStMapper::BayesianBioGeoStMapper(BioGeoTree *inbgt, Tree *intr,
                                               RateModel *inrm, bool marg,
                                               int gen) {
  gens = gen;
  marginal = marg;
  bgt = inbgt;
  tree = intr;
  rm = inrm;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  sts = "states";
  pts = "points";
}

void BayesianBioGeoStMapper::run_mappings() {
  /*
     """
      Stochastically map a binary character on a set of trees.

      newick_iter
        * a sequence of newick trees (for newick in newick_iter: ...)

      name2state
        * mapping of tip name -> state (0 or 1)

      treelen_prior_func
        * function to get another draw from the treelength prior

      pi0_prior_func
        * ditto for pi0, the stationary frequency of state 0 and
          so-called 'bias' parameter

      outfile, histfile:
        * open files for logging posterior parameter values (pi0,
          treelength, etc) and the actual simulation results (see map())

      sims_per_tree:
        * number of simulations to save for each tree

      ttable:
        * optional mapping of tip names to actual names in name2state

      char_name:
        * optional name of the character, for console logging

      thin:
        * skip to very ith tree

      silent:
        * suppress logging to console
      """

      states = (0, 1)
      prior = (1.0, 1.0)  # priors for states

      print >> outfile, "pi0\ttreelen\tgains\tlosses\tprop0"

      treenum = 1
      for nwck in newick_iter:
          if treenum % thin == 0:
              t = newick.parse(nwck, ttable=ttable)
              phylo.polarize(t)
              nodes = [ n for n in t.descendants(phylo.PREORDER) ]
              for i, n in enumerate(nodes):
                  n.nodeNum = i
              origlens = [ (n.length or 0) for n in nodes ]
              totlen = sum([ n.length for n in nodes if n.back ])

              # this is inefficient
              nodeNum2state = dict([ (n.nodeNum, name2state[n.label]) \
                                     for n in t.leaves() ])

              preorder = range(len(nodes))
              postorder = [ n.nodeNum for n in t.descendants(phylo.POSTORDER) ]

              sims = 0
              while sims < sims_per_tree:
                  treelen = treelen_prior_func()
                  pi0 = pi0_prior_func()

                  try:
                      result = map(nodes, preorder, postorder, nodeNum2state,
                                   states, prior, pi0, treelen,
                                   totlen, origlens)
                  except OverflowError:
                      continue

                  if result:
                      sims += 1
                      changes, times = summarize(result, states)

                      line = "%g\t%g\t%g\t%g\t%g" % \
                             (pi0, treelen, changes[0][1],
                              changes[1][0], times[0])

                      print >> outfile, line
                      print >> histfile, repr(result)

                      if not silent:
                          if outfile is not sys.stdout:
                              if char_name:
                                  print "%s %s: %s" % (char_name, treenum, line)
                              else:
                                  print treenum, sims, line

              # clean up
              for node in nodes: node.unlink()
              gc.collect()

          treenum += 1
*/
  double sum = 0;
  for (int i = 0; i < tree->getNodeCount(); i++) {
    sum += (tree->getNode(i)->getBL());
  }
  double totlen = sum;
  int sims = 0;
  int sims_per_tree = 2;
  while (sims < sims_per_tree) {
    // cout << sims << endl;
    double treelen = 1.0; // treelen_prior_func();
    double pi0 = 1.0;     // pi0_prior_func();
    // cout << "TESTING" << endl;
    initialize_nodes();
    bool result = mapping(treelen, totlen);
    if (result == true)
      cout << "success" << endl;
    sims += 1;
    if (result == true) {
      sims += 1;
      for (int i = 0; i < tree->getNodeCount(); i++) {
        VectorNodeObject<int> *newstates =
            ((VectorNodeObject<int> *)tree->getNode(i)->getObject(sts));
        VectorNodeObject<double> *newpoints =
            ((VectorNodeObject<double> *)tree->getNode(i)->getObject(pts));

        for (unsigned int j = 0; j < newstates->size(); j++) {
          cout << tree->getNode(i)->getName() << " " << j << " "
               << newstates->at(j) << " " << newpoints->at(j) << endl;
        }
      }
      // changes, times = summarize(result, states);

      //			line = "%g\t%g\t%g\t%g\t%g" % \
//					(pi0, treelen, changes[0][1],
      //							changes[1][0],
      //times[0]);

      //			cout << outfile, line;
      //			print >> histfile, repr(result);
    }
  }
}

bool BayesianBioGeoStMapper::mapping(double treelen, double totlen) {
  map<Node *, double> origlens;
  /*
   *     lenscale = treelen/totlen

      for n in nodes:
          n.length = (n.length or 0) * lenscale
  */
  double lenscale = treelen / totlen;
  for (int i = 0; i < tree->getNodeCount(); i++) {
    origlens[tree->getNode(i)] = tree->getNode(i)->getBL();
    tree->getNode(i)->setBL(tree->getNode(i)->getBL() * lenscale);
  }

  /*
Q = rates.Q2(pi0)
#Qmap = rates.Q2Pdict(nodes, Q)
Qmap = rates.binPdict(nodes, pi0)

tfracs = rates.fractionals(nodes, postorder, data,
                         states, state_priors, Qmap)

ancs = rates.sample_ancstates(nodes, preorder, states, tfracs, Qmap)
*/

  sample_ancstates();
  cout << "TESTING2" << endl;
  /*
result = simulate_on_nodes(nodes, ancs, Q, states)

for i, n in enumerate(nodes):
  n.length = origlens[i]

return result
*/
  bool result = simulate_on_nodes();
  if (result == true)
    cout << "success" << endl;
  cout << "TESTING3" << endl;

  for (int i = 0; i < tree->getNodeCount(); i++) {
    tree->getNode(i)->setBL(origlens[tree->getNode(i)]);
  }
  return result;
}

void BayesianBioGeoStMapper::sample_ancstates() {
  /*
   *     """
      Sample ancestral states from their conditional probabilities.
      Return a mapping of nodeNum -> ancstate
      """
      ancstates = {}
      for n in [ nodes[i] for i in preorder ]: #if not nodes[i].istip ]:
          nodefracs = fractionals[n.nodeNum]

          if n.parent:
              P = Qmap[n.length]
              ancst = ancstates[n.parent.nodeNum]
              newstate_Prow = P[ancst]
              nodefracs = nodefracs * newstate_Prow
              nodefracs /= sum(nodefracs)

          rv = uniform()
          v = 0.0
          for state, frac in zip(states, nodefracs):
              v += frac
              if rv < v:
                  break
          ancstates[n.nodeNum] = state

      return ancstates
   */

  /*
   * need to do eval_likelihood and see if i can get the segment
   * conditionals - then sample the ancestral splits and associate
   * them with
   * int si = startstate_node_map[tree->getNode(i)->getParent()];
                  int sj = startstate_node_map[tree->getNode(i)];
   */
  // preorder
  stack<Node *> preorder_nds;
  preorder_nds.push(tree->getRoot());
  while (preorder_nds.empty() == false) {
    Node *tnode = preorder_nds.top();
    preorder_nds.pop();
    for (int i = 0; i < tnode->getChildCount(); i++) {
      preorder_nds.push(&tnode->getChild(i));
    }
    VectorNodeObject<double> *distconds;
    if (tnode != tree->getRoot()) {
      VectorNodeObject<BranchSegment> *tsegs =
          ((VectorNodeObject<BranchSegment> *)tnode->getObject("segments"));
      distconds = tsegs->at(0).distconds;
    } else {
      distconds =
          ((VectorNodeObject<double> *)tnode->getObject("dist_conditionals"));
    }
    VectorNodeObject<double> usedistconds(distconds->size(), 0);
    double sum = 0;
    for (unsigned int j = 0; j < usedistconds.size(); j++) {
      sum += distconds->at(j);
    }
    for (unsigned int j = 0; j < usedistconds.size(); j++) {
      usedistconds[j] = distconds->at(j) / sum;
    }
    if (tnode != tree->getRoot()) {
      vector<vector<double>> P = rm->setup_fortran_P(0, tnode->getBL(), false);
      // int sj = startstate_node_map[tnode->getParent()];
      // because of how states are chosen it is actually the startstatefor the
      // node
      int sj = startstate_node_map[tnode];
      sum = 0;
      for (unsigned int j = 0; j < usedistconds.size(); j++) {
        usedistconds[j] *= P[sj][j];
        sum += usedistconds[j];
      }
      for (unsigned int j = 0; j < usedistconds.size(); j++) {
        usedistconds[j] /= sum;
      }
    }
    double rv = gsl_ran_flat(r, 0.0, 1.0);
    double v = 0.0;
    /*
     * need to deal with the splits
     */
    cout << "---" << endl;
    for (unsigned int j = 0; j < usedistconds.size(); j++) {
      v += usedistconds[j];
      if (rv < v) {
        // startstate_node_map[tnode] = j;
        endstate_node_map[tnode] = j;
        // randomly choose a speciation model
        if (tnode->getChildCount() > 0) {
          double weight;
          vector<int> leftdists;
          vector<int> rightdists;
          iter_ancsplits_just_int(rm, rm->get_int_dists_map()->at(j), leftdists,
                                  rightdists, weight);
          int splsd = 0;
          if (leftdists.size() > 1) {
            double rv2 = gsl_ran_flat(r, 0.0, 1.0);
            double vc2 = 0.0;
            vector<double> tusespls(leftdists.size(), 0);
            /*
             * get the splits numbers
             */
            VectorNodeObject<double> v1 =
                bgt->conditionals(tnode->getChild(0), marginal, false);
            VectorNodeObject<double> v2 =
                bgt->conditionals(tnode->getChild(1), marginal, false);
            double sum2 = 0;
            for (unsigned int m = 0; m < leftdists.size(); m++) {
              int ind1 = leftdists[m];
              int ind2 = rightdists[m];
              double lh_part = v1.at(ind1) * v2.at(ind2);
              tusespls[m] = (lh_part * weight);
              sum2 += (lh_part * weight);
            }
            for (unsigned int m = 0; m < leftdists.size(); m++) {
              tusespls[m] /= sum2;
            }

            for (unsigned int k = 0; k < leftdists.size(); k++) {
              vc2 += tusespls[k];
              if (rv2 < vc2) {
                splsd = k;
                break;
              }
            }
          }
          cout << "sample: " << tnode->getName() << " " << splsd << " " << j
               << " lft: " << leftdists[splsd] << " rt: " << rightdists[splsd]
               << endl;
          startstate_node_map[&tnode->getChild(0)] =
              leftdists[splsd]; // instead of parent
          startstate_node_map[&tnode->getChild(1)] =
              rightdists[splsd]; // instead of parent
        }
        break;
      }
    }

    cout << rv << endl;
    print_vector_double(usedistconds);
    cout << "ancest:" << tnode->getName()
         << " parspl:" << startstate_node_map[tnode]
         << " anc:" << endstate_node_map[tnode] << endl;
    cout << "---" << endl;
  }
}

void BayesianBioGeoStMapper::initialize_nodes() {
  for (int i = 0; i < tree->getNodeCount(); i++) {
    VectorNodeObject<int> *tsts = new VectorNodeObject<int>();
    tree->getNode(i)->assocObject(sts, *tsts);
    VectorNodeObject<double> *tpts = new VectorNodeObject<double>();
    tree->getNode(i)->assocObject(pts, *tpts);
  }
}

void BayesianBioGeoStMapper::clean_nodes() {
  for (int i = 0; i < tree->getNodeCount(); i++) {
    delete tree->getNode(i)->getObject(sts);
    delete tree->getNode(i)->getObject(pts);
  }
}

bool BayesianBioGeoStMapper::simulate_on_nodes() {
  // Like simulate_on_tree, but the vector of nodes is precalculated

  /*   results = {}
     for node in nodes:
         if node.parent:
             si = data[node.parent.nodeNum]
             sj = data[node.nodeNum]
             res = simulate_on_branch(si, sj, Q, states, node.length)
             if res:
                 results[node.nodeNum] = res
             else:
                 if DEBUG: print "*** break on node %d" % node.nodeNum
 ##                 print "*** break on node %d, %g %s %s" %\
 ##                       (node.nodeNum, node.length, si, sj)
                 results = None
                 break
     return results*/
  bool res = false;
  for (int i = 0; i < tree->getNodeCount(); i++) {
    if (tree->getNode(i)->hasParent()) {
      // int si = startstate_node_map[tree->getNode(i)->getParent()];
      // int sj = startstate_node_map[tree->getNode(i)];
      int si = startstate_node_map[tree->getNode(i)];
      int sj = endstate_node_map[tree->getNode(i)];
      cout << "simonnodes: " << tree->getNode(i)->getName() << " start:" << si
           << " end:" << sj << endl;
      VectorNodeObject<int> *newstates =
          ((VectorNodeObject<int> *)tree->getNode(i)->getObject(sts));
      VectorNodeObject<double> *newpoints =
          ((VectorNodeObject<double> *)tree->getNode(i)->getObject(pts));
      res = simulate_on_branch(si, sj, tree->getNode(i)->getBL(), newstates,
                               newpoints);
      // cout << res << endl;
      if (res == false) {
        newstates->clear();
        newpoints->clear();
        break;
      }
    }
  }
  return res;
}

bool BayesianBioGeoStMapper::simulate_on_branch(int starting_state,
                                                int end_state, double brlen,
                                                vector<int> *newstates,
                                                vector<double> *newpoints) {
  /*
Simulate evolution from state si to sj given the rate matrix Q
along the branch of length brlen.  Return a list of (state, len)
tuples that track the history of changes from 0 to brlen, or None
if the simulation was unsuccessful.
   */
  // osi = si; osj = sj
  int osi = starting_state;
  int osj = end_state;

  // point = 0.0
  double point = 0.0;

  // history = [(si, point)]  # a growing list of (state, changepoint) tuples
  //                         # tracking the character along the branch
  newstates->clear();
  newpoints->clear();
  newstates->push_back(osi);
  newpoints->push_back(point);
  /*
if si != sj:  # condition on one change occurring
  lambd = -(Q[si,si])
  U = uniform(0.0, 1.0)
  # see appendix of Nielsen 2001, Genetics
  t = brlen - point
  newpoint = -(1.0/lambd) * log(1.0 - U*(1.0 - exp(-lambd * t)))
  newstate = draw_new_state(Q, si, states)
  history.append((newstate, newpoint))
  si = newstate; point = newpoint
   */
  cout << starting_state << " " << end_state << endl;
  if (starting_state != end_state) {
    double lambd = -(rm->get_Q()[0][osi][osi]);
    double uni = gsl_ran_flat(r, 0.0, 1.0);
    double t = brlen - point;
    double newpoint = -(1.0 / lambd) * log(1.0 - uni * (1.0 - exp(-lambd * t)));
    int newstate = draw_new_state(osi);
    newstates->push_back(newstate);
    newpoints->push_back(newpoint);
    osi = newstate;
    point = newpoint;
  }
  cout << "osi: " << osi << endl;
  /*
while 1:
  lambd = -(Q[si,si])
  rv = expovariate(lambd)
  newpoint = point + rv

  if newpoint <= brlen:  # state change along branch
      newstate = draw_new_state(Q, si, states)
      history.append((newstate, newpoint))
      si = newstate; point = newpoint
  else:
      history.append((si, brlen))
      break
   */
  bool keepgoing = true;
  while (keepgoing) {
    double lambd = -(rm->get_Q()[0][osi][osi]);
    double rv = gsl_ran_exponential(r, lambd);
    cout << rv << endl;
    double newpoint = point + rv;
    if (newpoint <= brlen) {
      int newstate = draw_new_state(osi);
      newstates->push_back(newstate);
      newpoints->push_back(newpoint);
      osi = newstate;
      point = newpoint;
      cout << "newstate: " << newstate << " newpoint: " << newpoint << endl;
    } else {
      newstates->push_back(osi);
      newpoints->push_back(brlen);
      keepgoing = false;
    }
  }
  /*
if si == sj or (not condition_on_success): # success
  #history.append((si, brlen))
  return history
else:
  pass
  */
  // true is success
  if (osi == osj) {
    return true;
  } else {
    return false;
  }
}

/*
 * return the state number
 */
int BayesianBioGeoStMapper::draw_new_state(int starting_state) {
  // cout << " - " << starting_state << endl;
  /*
Given a rate matrix Q, a starting state si, and an ordered
sequence of states, eg (0, 1), draw a new state sj with
probability -(qij/qii)
   */
  // Qrow = Q[si]
  // qii = Qrow[si]
  // qij_probs = [ (x, -(Qrow[x]/qii)) for x in states if x != si ]
  double qii = rm->get_Q()[0][starting_state][starting_state];
  std::map<int, double> qij_probs;
  double sum = 0;
  for (unsigned int i = 1; i < rm->getDists()->size(); i++) {
    // int distint = rm->get_dists_int_map()->at(rm->getDists()->at(i));
    int distint = i;
    cout << "distint: " << distint << endl;
    if (starting_state != distint) {
      qij_probs[distint] = -rm->get_Q()[0][starting_state][distint] / qii;
      sum += -rm->get_Q()[0][starting_state][distint] / qii;
    }
  }
  cout << "sum:" << sum << endl;
  // uni = uniform(0.0, 1.0)
  // val = 0.0
  int state;
  double uni = gsl_ran_flat(r, 0.0, 1.0);
  double val = 0.0;
  std::map<int, double>::iterator it;
  // for sj, prob in qij_probs:
  //    val += prob
  //    if uni < val:
  //        return sj
  for (it = qij_probs.begin(); it != qij_probs.end(); it++) {
    val += (*it).second / sum;
    cout << "its:" << (*it).second / sum << " " << (*it).first << endl;
    if (uni < val) {
      state = (*it).first;
      return state;
    }
  }
  cout << "state: " << state << endl;
  return state;
}
