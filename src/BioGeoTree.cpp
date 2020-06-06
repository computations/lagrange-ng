/*
 * BioGeoTree.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */
#include <algorithm>
#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace std;

#include <armadillo>

using namespace arma;

#include "AncSplit.h"
#include "BioGeoTree.h"
#include "BioGeoTreeTools.h"
#include "BranchSegment.h"
#include "RateMatrixUtils.h"
#include "RateModel.h"

#include "node.h"
#include "tree.h"
#include "vector_node_object.h"

//#include "omp.h"
// octave usage
//#include <octave/oct.h>

Superdouble MAX(const Superdouble &a, const Superdouble &b) {
  return b > a ? b : a;
}

/*
 * sloppy beginning but best for now because of the complicated bits
 */

BioGeoTree::BioGeoTree(Tree *tr, vector<double> ps)
    : tree(tr), periods(ps), age("age"), dc("dist_conditionals"),
      en("excluded_dists"), andc("anc_dist_conditionals"), columns(NULL),
      whichcolumns(NULL), rootratemodel(NULL), distmap(NULL),
      store_p_matrices(false), use_stored_matrices(false), revB("revB"),
      rev(false), rev_exp_number("rev_exp_number"),
      rev_exp_time("rev_exp_time"), stochastic(false),
      stored_EN_matrices(unordered_map<int, unordered_map<double, mat>>()),
      stored_EN_CX_matrices(
          unordered_map<int, unordered_map<double, cx_mat>>()),
      stored_ER_matrices(unordered_map<int, unordered_map<double, mat>>()) {

  /*
   * initialize each node with segments
   */
  cout << "initializing nodes..." << endl;
  for (int i = 0; i < tree->getNodeCount(); i++) {
    if (tree->getNode(i)->getBL() < 0.000001)
      tree->getNode(i)->setBL(0.000001);
    tree->getNode(i)->initSegVector();
    tree->getNode(i)->initExclDistVector();
  }
  /*
   * initialize the actual branch segments for each node
   */
  tree->setHeightFromTipToNodes();
  cout << "initializing branch segments..." << endl;
  for (int i = 0; i < tree->getNodeCount(); i++) {
    if (tree->getNode(i)->hasParent()) {
      vector<double> pers(periods);
      double anc = tree->getNode(i)->getParent()->getHeight();
      double des = tree->getNode(i)->getHeight();
      // assert anc > des:q
      double t = des;
      if (pers.size() > 0) {
        for (unsigned int j = 0; j < pers.size(); j++) {
          double s = 0;
          if (pers.size() == 1)
            s = pers[0];
          for (unsigned int k = 0; k < j + 1; k++) {
            s += pers[k];
          }
          if (t < s) {
            double duration = min(s - t, anc - t);
            if (duration > 0) {
              BranchSegment tseg = BranchSegment(duration, j);
              tree->getNode(i)->getSegVector()->push_back(tseg);
            }
            t += duration; // TODO: make sure that this is all working
          }
          if (t > anc || pers[j] > t) {
            break;
          }
        }
      } else {
        BranchSegment tseg = BranchSegment(tree->getNode(i)->getBL(), 0);
        tree->getNode(i)->getSegVector()->push_back(tseg);
      }
    }
  }
}

void BioGeoTree::set_store_p_matrices(bool i) { store_p_matrices = i; }

void BioGeoTree::set_use_stored_matrices(bool i) { use_stored_matrices = i; }

void BioGeoTree::set_default_model(RateModel *mod) {
  rootratemodel = mod;
  for (int i = 0; i < tree->getNodeCount(); i++) {
    vector<BranchSegment> *tsegs = tree->getNode(i)->getSegVector();
    for (unsigned int j = 0; j < tsegs->size(); j++) {
      tsegs->at(j).setModel(mod);
      vector<Superdouble> *distconds =
          new vector<Superdouble>(rootratemodel->getDists()->size(), 0);
      tsegs->at(j).distconds = distconds;
      vector<Superdouble> *ancdistconds =
          new vector<Superdouble>(rootratemodel->getDists()->size(), 0);
      tsegs->at(j).ancdistconds = ancdistconds;
    }
  }
  vector<Superdouble> *distconds =
      new vector<Superdouble>(rootratemodel->getDists()->size(), 0);
  tree->getRoot()->assocDoubleVector(dc, *distconds);
  delete distconds;
  vector<Superdouble> *ancdistconds =
      new vector<Superdouble>(rootratemodel->getDists()->size(), 0);
  tree->getRoot()->assocDoubleVector(andc, *ancdistconds);
  delete ancdistconds;
}

void BioGeoTree::update_default_model(RateModel *mod) {
  rootratemodel = mod;

  for (int i = 0; i < tree->getNodeCount(); i++) {
    vector<BranchSegment> *tsegs = tree->getNode(i)->getSegVector();
    for (unsigned int j = 0; j < tsegs->size(); j++) {
      tsegs->at(j).setModel(mod);
    }
  }
}

void BioGeoTree::set_tip_conditionals(
    unordered_map<string, vector<int>> distrib_data) {
  int numofleaves = tree->getExternalNodeCount();
  for (int i = 0; i < numofleaves; i++) {
    vector<BranchSegment> *tsegs = tree->getExternalNode(i)->getSegVector();
    RateModel *mod = tsegs->at(0).getModel();
    int ind1 = get_vector_int_index_from_multi_vector_int(
        &distrib_data[tree->getExternalNode(i)->getName()], mod->getDists());
    tsegs->at(0).distconds->at(ind1) = 1.0;
  }
}

void BioGeoTree::set_excluded_dist(vector<int> ind, Node *node) {
  node->getExclDistVector()->push_back(ind);
}

Superdouble BioGeoTree::eval_likelihood(bool marginal) {
  if (rootratemodel->sparse == true) {
    columns = new vector<int>(rootratemodel->getDists()->size());
    whichcolumns = new vector<int>();
  }
  ancdist_conditional_lh(*tree->getRoot(), marginal);
  if (rootratemodel->sparse == true) {
    delete columns;
    delete whichcolumns;
  }
  // return (-(log(calculate_vector_double_sum(*(vector<double>*)
  // tree->getRoot()->getDoubleVector(dc)))));
  return -(calculate_vector_Superdouble_sum(
               *(vector<Superdouble> *)tree->getRoot()->getDoubleVector(dc)))
              .getLn();
}

vector<Superdouble> BioGeoTree::conditionals(Node &node, bool marginal,
                                             bool sparse) {
  vector<Superdouble> distconds;
  vector<BranchSegment> *tsegs = node.getSegVector();

  distconds = *tsegs->at(0).distconds;
  for (unsigned int i = 0; i < tsegs->size(); i++) {
    for (unsigned int j = 0; j < distconds.size(); j++) {
      tsegs->at(i).distconds->at(j) = distconds.at(j);
    }
    RateModel *rm = tsegs->at(i).getModel();
    vector<Superdouble> *v =
        new vector<Superdouble>(rootratemodel->getDists()->size(), 0);
    vector<int> distrange;
    if (tsegs->at(i).get_start_dist_int() != -666) {
      int ind1 = tsegs->at(i).get_start_dist_int();
      distrange.push_back(ind1);
    } else if (tsegs->at(i).getFossilAreas().size() > 0) {
      for (unsigned int j = 0; j < rootratemodel->getDists()->size(); j++) {
        distrange.push_back(j);
      }
      for (unsigned int k = 0; k < distrange.size(); k++) {
        bool flag = true;
        for (unsigned int x = 0; x < tsegs->at(i).getFossilAreas().size();
             x++) {
          if (tsegs->at(i).getFossilAreas()[x] == 1 && distrange.at(x) == 0) {
            flag = false;
          }
        }
        if (flag == true) {
          distrange.erase(distrange.begin() + k);
        }
      }
    } else {
      for (unsigned int j = 0; j < rootratemodel->getDists()->size(); j++) {
        distrange.push_back(j);
      }
    }
    /*
     * marginal
     */
    if (marginal == true) {
      if (sparse == false) {
        vector<vector<double>> p;
        if (use_stored_matrices == false) {
          p = rm->setup_fortran_P(tsegs->at(i).getPeriod(),
                                  tsegs->at(i).getDuration(), store_p_matrices);
        } else {
          p = rm->stored_p_matrices[tsegs->at(i).getPeriod()]
                                   [tsegs->at(i).getDuration()];
        }
        for (unsigned int j = 0; j < distrange.size(); j++) {
          for (unsigned int k = 0; k < distconds.size(); k++) {
            v->at(distrange[j]) += (distconds.at(k) * p[distrange[j]][k]);
          }
        }
      } else { // sparse
        /*
          testing pthread version
        */
        if (rm->get_nthreads() > 0) {
          /*commented so this can compile in MAC
           * vector<vector<double > > p =
          rm->setup_pthread_sparse_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),*whichcolumns);
          for(unsigned int j=0;j<distrange.size();j++){
          for(unsigned int k=0;k<distconds.size();k++){
              v->at(distrange[j]) += (distconds.at(k)*p[distrange[j]][k]);
          }
          }*/
        } else {
          for (unsigned int j = 0; j < distrange.size(); j++) {
            bool inthere = false;
            if (columns->at(distrange[j]) == 1)
              inthere = true;
            vector<double> p;
            if (inthere == true) {
              p = rm->setup_sparse_single_column_P(tsegs->at(i).getPeriod(),
                                                   tsegs->at(i).getDuration(),
                                                   distrange[j]);
            } else {
              p = vector<double>(distconds.size(), 0);
            }
            for (unsigned int k = 0; k < distconds.size(); k++) {
              v->at(distrange[j]) += (distconds.at(k) * p[k]);
            }
          }
        }
      }
    }
    /*
     * joint reconstruction
     * NOT FINISHED YET -- DONT USE
     */
    else {
      if (sparse == false) {
        vector<vector<double>> p =
            rm->setup_fortran_P(tsegs->at(i).getPeriod(),
                                tsegs->at(i).getDuration(), store_p_matrices);
        for (unsigned int j = 0; j < distrange.size(); j++) {
          Superdouble maxnum = 0;
          for (unsigned int k = 0; k < distconds.size(); k++) {
            Superdouble tx = (distconds.at(k) * p[distrange[j]][k]);
            maxnum = MAX(tx, maxnum);
          }
          v->at(distrange[j]) = maxnum;
        }
      } else { // sparse
      }
    }
    for (unsigned int j = 0; j < distconds.size(); j++) {
      distconds[j] = v->at(j);
    }
    if (store_p_matrices == true) {
      tsegs->at(i).seg_sp_alphas = distconds;
    }
    delete v;
  }
  /*
   * if store is true we want to store the conditionals for each node
   * for possible use in ancestral state reconstruction
   */
  if (store_p_matrices == true) {
    tsegs->at(0).alphas = distconds;
  }
  return distconds;
}

void BioGeoTree::ancdist_conditional_lh(Node &node, bool marginal) {
  vector<Superdouble> distconds(rootratemodel->getDists()->size(), 0);
  if (node.isExternal() == false) { // is not a tip
    Node *c1 = &node.getChild(0);
    Node *c2 = &node.getChild(1);
    RateModel *model;
    if (node.hasParent() == true) {
      vector<BranchSegment> *tsegs = node.getSegVector();
      model = tsegs->at(0).getModel();
    } else {
      model = rootratemodel;
    }
    ancdist_conditional_lh(*c1, marginal);
    ancdist_conditional_lh(*c2, marginal);
    bool sparse = rootratemodel->sparse;
    vector<Superdouble> v1;
    vector<Superdouble> v2;
    if (sparse == true) {
      // getcolumns
      vector<BranchSegment> *c1tsegs = c1->getSegVector();
      vector<BranchSegment> *c2tsegs = c2->getSegVector();
      vector<int> lcols =
          get_columns_for_sparse(*c1tsegs->at(0).distconds, rootratemodel);
      vector<int> rcols =
          get_columns_for_sparse(*c2tsegs->at(0).distconds, rootratemodel);
      whichcolumns->clear();
      for (unsigned int i = 0; i < lcols.size(); i++) {
        if (lcols[i] == 1 || rcols[i] == 1) {
          columns->at(i) = 1;
          if (i != 0 &&
              count(whichcolumns->begin(), whichcolumns->end(), i) == 0)
            whichcolumns->push_back(i);
        } else {
          columns->at(i) = 0;
        }
      }
      if (calculate_vector_int_sum(columns) == 0) {
        for (unsigned int i = 0; i < lcols.size(); i++) {
          columns->at(i) = 1;
        }
      }
      columns->at(0) = 0;
    }

    v1 = conditionals(*c1, marginal, sparse);
    v2 = conditionals(*c2, marginal, sparse);

    vector<vector<int>> *dists = rootratemodel->getDists();
    vector<int> leftdists;
    vector<int> rightdists;
    double weight;
    // cl1 = clock();

    for (unsigned int i = 0; i < dists->size(); i++) {

      if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
        Superdouble lh = 0.0;
        vector<vector<int>> *exdist = node.getExclDistVector();
        int cou = count(exdist->begin(), exdist->end(), dists->at(i));
        if (cou == 0) {
          iter_ancsplits_just_int(rootratemodel, dists->at(i), leftdists,
                                  rightdists, weight);
          for (unsigned int j = 0; j < leftdists.size(); j++) {
            int ind1 = leftdists[j];
            int ind2 = rightdists[j];
            Superdouble lh_part = v1.at(ind1) * v2.at(ind2);
            lh += (lh_part * weight);
          }
        }
        distconds.at(i) = lh;
      }
    }
    /// cl2 = clock();
    // ti += cl2-cl1;
  } else {
    vector<BranchSegment> *tsegs = node.getSegVector();
    distconds = *tsegs->at(0).distconds;
  }
  // testing scale
  // if (run_with_scale){
  //    scale_node(&distconds);
  //}
  // testing scale
  if (node.hasParent() == true) {
    vector<BranchSegment> *tsegs = node.getSegVector();
    for (unsigned int i = 0; i < distconds.size(); i++) {
      tsegs->at(0).distconds->at(i) = distconds.at(i);
    }
  } else {
    for (unsigned int i = 0; i < distconds.size(); i++) {
      node.getDoubleVector(dc)->at(i) = distconds.at(i);
      // cout << distconds.at(i) << endl;
    }
  }
}

/*
 * ********************************************
 *
 * adds fossils either at the node or along a branch
 *
 * ********************************************
 */
void BioGeoTree::setFossilatNodeByMRCA(vector<string> nodeNames,
                                       int fossilarea) {
  Node *mrca = tree->getMRCA(nodeNames);
  vector<vector<int>> *dists = rootratemodel->getDists();
  for (unsigned int i = 0; i < dists->size(); i++) {
    if (dists->at(i).at(fossilarea) == 0) {
      vector<vector<int>> *exd = mrca->getExclDistVector();
      exd->push_back(dists->at(i));
    }
  }
}

void BioGeoTree::setFossilatNodeByMRCA_id(Node *id, int fossilarea) {
  vector<vector<int>> *dists = rootratemodel->getDists();
  for (unsigned int i = 0; i < dists->size(); i++) {
    if (dists->at(i).at(fossilarea) == 0) {
      vector<vector<int>> *exd = id->getExclDistVector();
      exd->push_back(dists->at(i));
    }
  }
}

void BioGeoTree::setFossilatBranchByMRCA(vector<string> nodeNames,
                                         int fossilarea, double age) {
  Node *mrca = tree->getMRCA(nodeNames);
  vector<BranchSegment> *tsegs = mrca->getSegVector();
  double startage = mrca->getHeight();
  for (unsigned int i = 0; i < tsegs->size(); i++) {
    if (age > startage && age < (startage + tsegs->at(i).getDuration())) {
      tsegs->at(i).setFossilArea(fossilarea);
    }
    startage += tsegs->at(i).getDuration();
  }
}

void BioGeoTree::setFossilatBranchByMRCA_id(Node *id, int fossilarea,
                                            double age) {
  vector<BranchSegment> *tsegs = id->getSegVector();
  double startage = id->getHeight();

  for (unsigned int i = 0; i < tsegs->size(); i++) {
    if (age > startage && age < (startage + tsegs->at(i).getDuration())) {
      tsegs->at(i).setFossilArea(fossilarea);
    }
    startage += tsegs->at(i).getDuration();
  }
}

/************************************************************
 forward and reverse stuff for ancestral states
 ************************************************************/
// add joint
void BioGeoTree::prepare_ancstate_reverse() { reverse(*tree->getRoot()); }

/*
 * called from prepare_ancstate_reverse and that is all
 */
void BioGeoTree::reverse(Node &node) {
  rev = true;
  vector<Superdouble> *revconds =
      new vector<Superdouble>(rootratemodel->getDists()->size(),
                              0); // need to delete this at some point
  if (&node == tree->getRoot()) {
    for (unsigned int i = 0; i < rootratemodel->getDists()->size(); i++) {
      revconds->at(i) = 1.0; // prior
    }
    node.assocDoubleVector(revB, *revconds);
    delete revconds;
    for (int i = 0; i < node.getChildCount(); i++) {
      reverse(node.getChild(i));
    }
  } else {
    // else if(node.isExternal() == false){
    // calculate A i
    // sum over all alpha k of sister node of the parent times the priors of the
    // speciations (weights) times B of parent j
    vector<Superdouble> *parrev = node.getParent()->getDoubleVector(revB);
    vector<Superdouble> sisdistconds;
    if (&node.getParent()->getChild(0) != &node) {
      vector<BranchSegment> *tsegs =
          node.getParent()->getChild(0).getSegVector();
      sisdistconds = tsegs->at(0).alphas;
    } else {
      vector<BranchSegment> *tsegs =
          node.getParent()->getChild(1).getSegVector();
      sisdistconds = tsegs->at(0).alphas;
    }
    vector<vector<int>> *dists = rootratemodel->getDists();
    vector<int> leftdists;
    vector<int> rightdists;
    double weight;
    // cl1 = clock();
    vector<Superdouble> tempA(rootratemodel->getDists()->size(), 0);
    for (unsigned int i = 0; i < dists->size(); i++) {
      if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
        vector<vector<int>> *exdist = node.getExclDistVector();
        int cou = count(exdist->begin(), exdist->end(), dists->at(i));
        if (cou == 0) {
          iter_ancsplits_just_int(rootratemodel, dists->at(i), leftdists,
                                  rightdists, weight);
          // root has i, curnode has left, sister of cur has right
          for (unsigned int j = 0; j < leftdists.size(); j++) {
            int ind1 = leftdists[j];
            int ind2 = rightdists[j];
            tempA[ind1] += (sisdistconds.at(ind2) * weight * parrev->at(i));
          }
        }
      }
    }

    // now calculate node B
    vector<BranchSegment> *tsegs = node.getSegVector();
    vector<Superdouble> tempmoveA(tempA);
    // for(unsigned int ts=0;ts<tsegs->size();ts++){
    for (int ts = tsegs->size() - 1; ts != -1; ts--) {
      for (unsigned int j = 0; j < dists->size(); j++) {
        revconds->at(j) = 0;
      }
      RateModel *rm = tsegs->at(ts).getModel();
      vector<vector<double>> *p =
          &rm->stored_p_matrices[tsegs->at(ts).getPeriod()]
                                [tsegs->at(ts).getDuration()];
      mat *EN = NULL;
      mat *ER = NULL;
      vector<Superdouble> tempmoveAer(tempA);
      vector<Superdouble> tempmoveAen(tempA);
      if (stochastic == true) {
        // initialize the segment B's
        for (unsigned int j = 0; j < dists->size(); j++) {
          tempmoveAer[j] = 0;
        }
        for (unsigned int j = 0; j < dists->size(); j++) {
          tempmoveAen[j] = 0;
        }
        EN = &stored_EN_matrices[tsegs->at(ts).getPeriod()]
                                [tsegs->at(ts).getDuration()];
        ER = &stored_ER_matrices[tsegs->at(ts).getPeriod()]
                                [tsegs->at(ts).getDuration()];
        // cout << (*EN) << endl;
        cx_mat *EN_CX = NULL;
        EN_CX = &stored_EN_CX_matrices[tsegs->at(ts).getPeriod()]
                                      [tsegs->at(ts).getDuration()];
        // cout << (*EN_CX) << endl;
        // exit(0);
      }
      for (unsigned int j = 0; j < dists->size(); j++) {
        if (accumulate(dists->at(j).begin(), dists->at(j).end(), 0) > 0) {
          for (unsigned int i = 0; i < dists->size(); i++) {
            if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
              // cout << "here " << j << " " << i<< " " << node.getBL() << " "
              // << ts <<" " << tsegs->size() << endl;
              revconds->at(j) +=
                  tempmoveA[i] *
                  ((*p)[i][j]); // tempA needs to change each time
              // cout << "and" << endl;
              if (stochastic == true) {
                tempmoveAer[j] += tempmoveA[i] * (((*ER)(i, j)));
                tempmoveAen[j] += tempmoveA[i] * (((*EN)(i, j)));
              }
            }
          }
        }
      }
      for (unsigned int j = 0; j < dists->size(); j++) {
        tempmoveA[j] = revconds->at(j);
      }
      if (stochastic == true) {
        tsegs->at(ts).seg_sp_stoch_map_revB_time = tempmoveAer;
        tsegs->at(ts).seg_sp_stoch_map_revB_number = tempmoveAen;
      }
    }
    node.assocDoubleVector(revB, *revconds);
    delete revconds;
    for (int i = 0; i < node.getChildCount(); i++) {
      reverse(node.getChild(i));
    }
  }
}

/*
 * calculates the most likely split (not state) -- the traditional result for
 * lagrange
 */

unordered_map<vector<int>, vector<AncSplit>>
BioGeoTree::calculate_ancsplit_reverse(Node &node, bool marg) {
  vector<Superdouble> *Bs = node.getDoubleVector(revB);
  unordered_map<vector<int>, vector<AncSplit>> ret;
  for (unsigned int j = 0; j < rootratemodel->getDists()->size(); j++) {
    vector<int> dist = rootratemodel->getDists()->at(j);
    vector<AncSplit> ans = iter_ancsplits(rootratemodel, dist);
    if (node.isExternal() == false) { // is not a tip
      Node *c1 = &node.getChild(0);
      Node *c2 = &node.getChild(1);
      vector<BranchSegment> *tsegs1 = c1->getSegVector();
      vector<BranchSegment> *tsegs2 = c2->getSegVector();
      for (unsigned int i = 0; i < ans.size(); i++) {
        vector<vector<int>> *exdist = node.getExclDistVector();
        int cou =
            count(exdist->begin(), exdist->end(),
                  (*rootratemodel->get_int_dists_map())[ans[i].ancdistint]);
        if (cou == 0) {
          vector<Superdouble> v1 = tsegs1->at(0).alphas;
          vector<Superdouble> v2 = tsegs2->at(0).alphas;
          Superdouble lh = (v1[ans[i].ldescdistint] * v2[ans[i].rdescdistint] *
                            Bs->at(j) * ans[i].getWeight());
          ans[i].setLikelihood(lh);
          // cout << lh << endl;
        }
      }
    }
    ret[dist] = ans;
  }
  return ret;
}

/*
 * calculates the ancestral area over all the possible splits
 */
vector<Superdouble> BioGeoTree::calculate_ancstate_reverse(Node &node,
                                                           bool marg) {
  if (node.isExternal() == false) { // is not a tip
    vector<Superdouble> *Bs = node.getDoubleVector(revB);
    vector<vector<int>> *dists = rootratemodel->getDists();
    vector<int> leftdists;
    vector<int> rightdists;
    double weight;
    Node *c1 = &node.getChild(0);
    Node *c2 = &node.getChild(1);
    vector<BranchSegment> *tsegs1 = c1->getSegVector();
    vector<BranchSegment> *tsegs2 = c2->getSegVector();
    vector<Superdouble> v1 = tsegs1->at(0).alphas;
    vector<Superdouble> v2 = tsegs2->at(0).alphas;
    vector<Superdouble> LHOODS(dists->size(), 0);
    for (unsigned int i = 0; i < dists->size(); i++) {
      if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
        vector<vector<int>> *exdist = node.getExclDistVector();
        int cou = count(exdist->begin(), exdist->end(), dists->at(i));
        if (cou == 0) {
          iter_ancsplits_just_int(rootratemodel, dists->at(i), leftdists,
                                  rightdists, weight);
          for (unsigned int j = 0; j < leftdists.size(); j++) {
            int ind1 = leftdists[j];
            int ind2 = rightdists[j];
            LHOODS[i] += (v1.at(ind1) * v2.at(ind2) * weight);
          }
          LHOODS[i] *= Bs->at(i);
        }
      }
    }
    return LHOODS;
  }
}

/**********************************************************
 * forward and reverse stuff for stochastic mapping
 **********************************************************/

void BioGeoTree::prepare_stochmap_reverse_all_nodes(int from, int to) {
  stochastic = true;
  int ndists = rootratemodel->getDists()->size();

  // calculate and store local expectation matrix for each branch length
  //#pragma omp parallel for ordered num_threads(8)
  for (int k = 0; k < tree->getNodeCount(); k++) {
    // cout << k << " " << tree->getNodeCount() << endl;
    vector<BranchSegment> *tsegs = tree->getNode(k)->getSegVector();
    for (unsigned int l = 0; l < tsegs->size(); l++) {
      int per = (*tsegs)[l].getPeriod();
      double dur = (*tsegs)[l].getDuration();
      cx_mat eigvec(ndists, ndists);
      eigvec.fill(0);
      cx_mat eigval(ndists, ndists);
      eigval.fill(0);
      bool isImag =
          rootratemodel->get_eigenvec_eigenval_from_Q(&eigval, &eigvec, per);
      mat Ql(ndists, ndists);
      Ql.fill(0);
      Ql(from, to) = rootratemodel->get_Q()[per][from][to];
      mat W(ndists, ndists);
      W.fill(0);
      W(from, from) = 1;
      cx_mat summed(ndists, ndists);
      summed.fill(0);
      cx_mat summedR(ndists, ndists);
      summedR.fill(0);
      for (int i = 0; i < ndists; i++) {
        mat Ei(ndists, ndists);
        Ei.fill(0);
        Ei(i, i) = 1;
        cx_mat Si(ndists, ndists);
        Si = eigvec * Ei * inv(eigvec);
        for (int j = 0; j < ndists; j++) {
          cx_double dij = (eigval(i, i) - eigval(j, j)) * dur;
          mat Ej(ndists, ndists);
          Ej.fill(0);
          Ej(j, j) = 1;
          cx_mat Sj(ndists, ndists);
          Sj = eigvec * Ej * inv(eigvec);
          cx_double Iijt = 0;
          if (abs(dij) > 10) {
            Iijt = (exp(eigval(i, i) * dur) - exp(eigval(j, j) * dur)) /
                   (eigval(i, i) - eigval(j, j));
          } else if (abs(dij) < 10e-20) {
            Iijt = dur * exp(eigval(j, j) * dur) *
                   (1. + dij / 2. + pow(dij, 2.) / 6. + pow(dij, 3.) / 24.);
          } else {
            if (eigval(i, i) == eigval(j, j)) {
              // WAS Iijt = dur*exp(eigval(j,j)*dur)*expm1(dij)/dij;
              if (isImag)
                Iijt = dur * exp(eigval(j, j) * dur) * (exp(dij) - 1.) / dij;
              else
                Iijt = dur * exp(eigval(j, j) * dur) * (expm1(real(dij))) / dij;
            } else {
              // WAS Iijt = -dur*exp(eigval(i,i)*dur)*expm1(-dij)/dij;
              if (isImag)
                Iijt = -dur * exp(eigval(i, i) * dur) * (exp(-dij) - 1.) / dij;
              else
                Iijt =
                    -dur * exp(eigval(i, i) * dur) * (expm1(real(-dij))) / dij;
            }
          }
          summed += (Si * Ql * Sj * Iijt);
          summedR += (Si * W * Sj * Iijt);
        }
      }
      stored_EN_matrices[per][dur] = (real(summed));
      stored_EN_CX_matrices[per][dur] = summed;
      stored_ER_matrices[per][dur] = (real(summedR));
      // for(int i=0;i<ndists;i++){
      //    for (int j=0;j<ndists;j++){
      //        if (real(summed(i,j)) < 0){
      //            cout <<"N:" <<  summed << endl;
      //            cout << endl;
      //            exit(0);
      //        }
      //        if (real(summedR(i,j)) < 0){
      //            cout <<"R:" << summedR << endl;
      //            cout << endl;
      //            exit(0);
      //        }
      //    }
      //}
    }
  }
}

/*
 * called directly after reverse_stochastic
 */

vector<Superdouble> BioGeoTree::calculate_reverse_stochmap(Node &node,
                                                           bool time) {
  if (node.isExternal() == false) { // is not a tip
    vector<BranchSegment> *tsegs = node.getSegVector();
    vector<vector<int>> *dists = rootratemodel->getDists();
    vector<Superdouble> totalExp(dists->size(), 0);
    for (int t = 0; t < tsegs->size(); t++) {
      if (t == 0) {
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_number;
        vector<int> leftdists;
        vector<int> rightdists;
        double weight;
        Node *c1 = &node.getChild(0);
        Node *c2 = &node.getChild(1);
        vector<BranchSegment> *tsegs1 = c1->getSegVector();
        vector<BranchSegment> *tsegs2 = c2->getSegVector();
        vector<Superdouble> v1 = tsegs1->at(0).alphas;
        vector<Superdouble> v2 = tsegs2->at(0).alphas;
        vector<Superdouble> LHOODS(dists->size(), 0);
        for (unsigned int i = 0; i < dists->size(); i++) {
          if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
            vector<vector<int>> *exdist = node.getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists->at(i));
            if (cou == 0) {
              iter_ancsplits_just_int(rootratemodel, dists->at(i), leftdists,
                                      rightdists, weight);
              for (unsigned int j = 0; j < leftdists.size(); j++) {
                int ind1 = leftdists[j];
                int ind2 = rightdists[j];
                LHOODS[i] += (v1.at(ind1) * v2.at(ind2) * weight);
              }
              LHOODS[i] *= Bs.at(i);
            }
          }
        }
        for (int i = 0; i < dists->size(); i++) {
          totalExp[i] = LHOODS[i];
        }
      } else {
        vector<Superdouble> alphs = tsegs->at(t - 1).seg_sp_alphas;
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_number;
        vector<Superdouble> LHOODS(dists->size(), 0);
        for (unsigned int i = 0; i < dists->size(); i++) {
          if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
            vector<vector<int>> *exdist = node.getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists->at(i));
            if (cou == 0) {
              LHOODS[i] =
                  Bs.at(i) * (alphs[i]); // do i do this or do i do from i to j
            }
          }
        }
        for (int i = 0; i < dists->size(); i++) {
          totalExp[i] += LHOODS[i];
        }
      }
    }
    // not sure if this should return a Superdouble or not when doing a bigtree
    return totalExp;
  } else {
    vector<BranchSegment> *tsegs = node.getSegVector();
    vector<vector<int>> *dists = rootratemodel->getDists();
    vector<Superdouble> totalExp(dists->size(), 0);
    for (int t = 0; t < tsegs->size(); t++) {
      if (t == 0) {
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_number;
        vector<Superdouble> LHOODS(dists->size(), 0);
        for (unsigned int i = 0; i < dists->size(); i++) {
          if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
            vector<vector<int>> *exdist = node.getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists->at(i));
            if (cou == 0) {
              LHOODS[i] = Bs.at(i) * (tsegs->at(0).distconds->at(i));
            }
          }
        }
        for (int i = 0; i < dists->size(); i++) {
          totalExp[i] = LHOODS[i];
        }
      } else {
        vector<Superdouble> alphs = tsegs->at(t - 1).seg_sp_alphas;
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_number;
        vector<Superdouble> LHOODS(dists->size(), 0);
        for (unsigned int i = 0; i < dists->size(); i++) {
          if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
            vector<vector<int>> *exdist = node.getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists->at(i));
            if (cou == 0) {
              LHOODS[i] = Bs.at(i) * (alphs[i]);
            }
          }
        }
        for (int i = 0; i < dists->size(); i++) {
          totalExp[i] += LHOODS[i];
        }
      }
    }
    return totalExp;
  }
}

/**********************************************************
 * trash collection
 **********************************************************/
BioGeoTree::~BioGeoTree() {
  for (int i = 0; i < tree->getNodeCount(); i++) {
    vector<BranchSegment> *tsegs = tree->getNode(i)->getSegVector();
    for (unsigned int j = 0; j < tsegs->size(); j++) {
      delete tsegs->at(j).distconds;
      delete tsegs->at(j).ancdistconds;
    }
    tree->getNode(i)->deleteExclDistVector();
    if (rev == true && tree->getNode(i)->isInternal()) {
      tree->getNode(i)->deleteDoubleVector(revB);
    }
    tree->getNode(i)->deleteSegVector();
  }
  tree->getRoot()->deleteDoubleVector(dc);
  tree->getRoot()->deleteDoubleVector(andc);
  tree->getRoot()->deleteDoubleVector(revB);
}

size_t hash<std::vector<int>>::operator()(const std::vector<int> &vec) const {
  std::size_t seed = vec.size();
  for (auto &i : vec) {
    seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}
