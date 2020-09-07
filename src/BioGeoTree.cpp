/*
 * BioGeoTree.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */
#include <algorithm>
#include <cmath>
#include <ctime>
#include <exception>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace std;

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

BioGeoTree::BioGeoTree(std::shared_ptr<Tree> tr, vector<double> ps)
    : _tree(tr), _periods(ps), _columns(nullptr), _which_columns(nullptr),
      _root_ratemodel(nullptr), _store_p_matrices(false),
      _use_stored_matrices(false), _reverse(false), _stocastic(false),
      _stored_EN_matrices{}, _stored_EN_CX_matrices{}, _stored_ER_matrices{} {

  /*
   * initialize each node with segments
   */
  cout << "initializing nodes..." << endl;
  for (unsigned int i = 0; i < _tree->getNodeCount(); i++) {
    if (_tree->getNode(i)->getBL() < 0.000001)
      _tree->getNode(i)->setBL(0.000001);
    _tree->getNode(i)->initSegVector();
    _tree->getNode(i)->initExclDistVector();
  }
  /*
   * initialize the actual branch segments for each node
   */
  _tree->setHeightFromTipToNodes();
  cout << "initializing branch segments..." << endl;

  for (size_t i = 0; i < _tree->getInternalNodeCount(); ++i) {
    auto current_node = _tree->getInternalNode(i);
    for (auto child : current_node->getChildren()) {
      double ancestor_height = current_node->getHeight();
      double descendant_height = child->getHeight();
      double t = descendant_height;

      if (_periods.size() > 0) {
        double s = 0;
        if (_periods.size() == 1) {
          s = _periods[0];
        }
        for (size_t j = 0; j < _periods.size(); j++) {
          s += _periods[j];
          if (t < s) {
            double duration = min(s - t, ancestor_height - t);
            if (duration > 0) {
              BranchSegment tseg = BranchSegment(duration, j);
              child->getSegVector()->push_back(tseg);
            }
            t += duration; // TODO: make sure that this is all working
          }
          if (t > ancestor_height || _periods[j] > t) {
            break;
          }
        }
      } else {
        BranchSegment tseg = BranchSegment(current_node->getBL(), 0);
        current_node->getSegVector()->push_back(tseg);
      }
    }
  }

  /*
  for (int i = 0; i < _tree->getNodeCount(); i++) {
    if (_tree->getNode(i)->hasParent()) {
      double anc = _tree->getNode(i)->getParent()->getHeight();
      double des = _tree->getNode(i)->getHeight();
      // assert anc > des
      double t = des;
      if (_periods.size() > 0) {
        for (unsigned int j = 0; j < _periods.size(); j++) {
          double s = 0;
          if (_periods.size() == 1)
            s = _periods[0];
          for (unsigned int k = 0; k < j + 1; k++) {
            s += _periods[k];
          }
          if (t < s) {
            double duration = min(s - t, anc - t);
            if (duration > 0) {
              BranchSegment tseg = BranchSegment(duration, j);
              _tree->getNode(i)->getSegVector()->push_back(tseg);
            }
            t += duration; // TODO: make sure that this is all working
          }
          if (t > anc || _periods[j] > t) {
            break;
          }
        }
      } else {
        BranchSegment tseg = BranchSegment(_tree->getNode(i)->getBL(), 0);
        _tree->getNode(i)->getSegVector()->push_back(tseg);
      }
    }
  }
  */
}

void BioGeoTree::set_store_p_matrices(bool i) { _store_p_matrices = i; }

void BioGeoTree::set_use_stored_matrices(bool i) { _use_stored_matrices = i; }

void BioGeoTree::set_default_model(std::shared_ptr<RateModel> mod) {
  _root_ratemodel = mod;
  for (unsigned int i = 0; i < _tree->getNodeCount(); i++) {
    vector<BranchSegment> *tsegs = _tree->getNode(i)->getSegVector();
    for (unsigned int j = 0; j < tsegs->size(); j++) {
      tsegs->at(j).setModel(mod);
      vector<Superdouble> *distconds =
          new vector<Superdouble>(_root_ratemodel->getDistsSize(), 0);
      tsegs->at(j).distconds = distconds;
      vector<Superdouble> *ancdistconds =
          new vector<Superdouble>(_root_ratemodel->getDistsSize(), 0);
      tsegs->at(j).ancdistconds = ancdistconds;
    }
  }
  vector<Superdouble> *distconds =
      new vector<Superdouble>(_root_ratemodel->getDistsSize(), 0);
  _tree->getRoot()->assocDoubleVector(_dist_conditionals_key, *distconds);
  delete distconds;
  vector<Superdouble> *ancdistconds =
      new vector<Superdouble>(_root_ratemodel->getDistsSize(), 0);
  _tree->getRoot()->assocDoubleVector(_anc_dist_conditionals_key,
                                      *ancdistconds);
  delete ancdistconds;
}

void BioGeoTree::update_default_model(std::shared_ptr<RateModel> mod) {
  _root_ratemodel = mod;

  for (unsigned int i = 0; i < _tree->getNodeCount(); i++) {
    vector<BranchSegment> *tsegs = _tree->getNode(i)->getSegVector();
    for (unsigned int j = 0; j < tsegs->size(); j++) {
      tsegs->at(j).setModel(mod);
    }
  }
}

void BioGeoTree::set_tip_conditionals(
    unordered_map<string, vector<int>> distrib_data) {
  int numofleaves = _tree->getExternalNodeCount();
  for (int i = 0; i < numofleaves; i++) {
    vector<BranchSegment> *tsegs = _tree->getExternalNode(i)->getSegVector();
    std::shared_ptr<RateModel> mod = tsegs->at(0).getModel();
    int ind1 = get_vector_int_index_from_multi_vector_int(
        distrib_data[_tree->getExternalNode(i)->getName()], mod->getDists());
    tsegs->at(0).distconds->at(ind1) = 1.0;
  }
}

void BioGeoTree::set_excluded_dist(vector<int> ind,
                                   std::shared_ptr<Node> node) {
  node->getExclDistVector()->push_back(ind);
}

Superdouble BioGeoTree::eval_likelihood(bool marginal) {
  if (_root_ratemodel->_sparse == true) {
    _columns = std::make_shared<vector<int>>(_root_ratemodel->getDistsSize());
    _which_columns = std::make_shared<vector<int>>();
  }
  ancdist_conditional_lh(_tree->getRoot(), marginal);
  if (_root_ratemodel->_sparse == true) {
    if (_columns.unique()) {
      _columns.reset();
    }
    if (_which_columns.unique()) {
      _which_columns.reset();
    }
  }
  return -(calculate_vector_Superdouble_sum(
               *(vector<Superdouble> *)_tree->getRoot()->getDoubleVector(
                   _dist_conditionals_key)))
              .getLn();
}

vector<Superdouble> BioGeoTree::conditionals(std::shared_ptr<Node> node,
                                             bool marginal, bool sparse) {
  vector<Superdouble> distconds;
  vector<BranchSegment> *time_segments = node->getSegVector();

  distconds = *time_segments->at(0).distconds;

  for (unsigned int i = 0; i < time_segments->size(); i++) {

    for (unsigned int j = 0; j < distconds.size(); j++) {
      time_segments->at(i).distconds->at(j) = distconds.at(j);
    }

    std::shared_ptr<RateModel> rm = time_segments->at(i).getModel();
    vector<Superdouble> v(_root_ratemodel->getDistsSize(), 0);
    vector<int> distrange;

    if (time_segments->at(i).get_start_dist_int() != -666) {
      int ind1 = time_segments->at(i).get_start_dist_int();
      distrange.push_back(ind1);
    } else if (time_segments->at(i).getFossilAreas().size() > 0) {
      for (unsigned int j = 0; j < _root_ratemodel->getDistsSize(); j++) {
        distrange.push_back(j);
      }
      for (unsigned int k = 0; k < distrange.size(); k++) {
        bool flag = true;
        for (unsigned int x = 0;
             x < time_segments->at(i).getFossilAreas().size(); x++) {
          if (time_segments->at(i).getFossilAreas()[x] == 1 &&
              distrange.at(x) == 0) {
            flag = false;
          }
        }
        if (flag == true) {
          distrange.erase(distrange.begin() + k);
        }
      }
    } else {
      for (unsigned int j = 0; j < _root_ratemodel->getDistsSize(); j++) {
        distrange.push_back(j);
      }
    }
    /*
     * marginal
     */
    if (marginal == true) {
      if (sparse == false) {
        vector<vector<double>> p;
        if (_use_stored_matrices == false) {
          p = rm->setup_fortran_P(time_segments->at(i).getPeriod(),
                                  time_segments->at(i).getDuration(),
                                  _store_p_matrices);
        } else {
          p = rm->stored_p_matrices[time_segments->at(i).getPeriod()]
                                   [time_segments->at(i).getDuration()];
        }
        for (unsigned int j = 0; j < distrange.size(); j++) {
          for (unsigned int k = 0; k < distconds.size(); k++) {
            v[distrange[j]] += (distconds.at(k) * p[distrange[j]][k]);
          }
        }
      } else { // sparse
        /*
          testing pthread version
        */
        if (rm->get_nthreads() > 0) {
        } else {
          for (unsigned int j = 0; j < distrange.size(); j++) {
            bool inthere = false;
            if (_columns->at(distrange[j]) == 1)
              inthere = true;
            vector<double> p;
            if (inthere == true) {
              p = rm->setup_sparse_single_column_P(
                  time_segments->at(i).getPeriod(),
                  time_segments->at(i).getDuration(), distrange[j]);
            } else {
              p = vector<double>(distconds.size(), 0);
            }
            for (unsigned int k = 0; k < distconds.size(); k++) {
              v.at(distrange[j]) += (distconds.at(k) * p[k]);
            }
          }
        }
      }
    }
    for (unsigned int j = 0; j < distconds.size(); j++) {
      distconds[j] = v[j];
    }
    if (_store_p_matrices == true) {
      time_segments->at(i).seg_sp_alphas = distconds;
    }
  }
  /*
   * if store is true we want to store the conditionals for each node
   * for possible use in ancestral state reconstruction
   */
  if (_store_p_matrices == true) {
    time_segments->at(0).alphas = distconds;
  }
  return distconds;
}

void BioGeoTree::ancdist_conditional_lh(std::shared_ptr<Node> node,
                                        bool marginal) {
  vector<Superdouble> distconds(_root_ratemodel->getDistsSize(), 0);
  if (node->isExternal() == false) { // is not a tip
    std::shared_ptr<Node> c1 = node->getChild(0);
    std::shared_ptr<Node> c2 = node->getChild(1);
    /*
    RateModel *model;
    if (node.hasParent() == true) {
      vector<BranchSegment> *tsegs = node.getSegVector();
      model = tsegs->at(0).getModel();
    } else {
      model = _root_ratemodel;
    }
    */
    ancdist_conditional_lh(c1, marginal);
    ancdist_conditional_lh(c2, marginal);
    bool sparse = _root_ratemodel->_sparse;
    vector<Superdouble> v1;
    vector<Superdouble> v2;
    if (sparse == true) {
      // getcolumns
      vector<BranchSegment> *c1tsegs = c1->getSegVector();
      vector<BranchSegment> *c2tsegs = c2->getSegVector();
      vector<int> lcols =
          get_columns_for_sparse(*c1tsegs->at(0).distconds, _root_ratemodel);
      vector<int> rcols =
          get_columns_for_sparse(*c2tsegs->at(0).distconds, _root_ratemodel);
      _which_columns->clear();
      for (unsigned int i = 0; i < lcols.size(); i++) {
        if (lcols[i] == 1 || rcols[i] == 1) {
          _columns->at(i) = 1;
          if (i != 0 &&
              count(_which_columns->begin(), _which_columns->end(), i) == 0)
            _which_columns->push_back(i);
        } else {
          _columns->at(i) = 0;
        }
      }
      if (calculate_vector_int_sum(*_columns) == 0) {
        for (unsigned int i = 0; i < lcols.size(); i++) {
          _columns->at(i) = 1;
        }
      }
      _columns->at(0) = 0;
    }

    v1 = conditionals(c1, marginal, sparse);
    v2 = conditionals(c2, marginal, sparse);

    auto dists = _root_ratemodel->getDists();
    vector<int> leftdists;
    vector<int> rightdists;
    double weight;

    for (unsigned int i = 0; i < dists.size(); i++) {

      if (accumulate(dists.at(i).begin(), dists.at(i).end(), 0) > 0) {
        Superdouble lh = 0.0;
        vector<vector<int>> *exdist = node->getExclDistVector();
        int cou = count(exdist->begin(), exdist->end(), dists.at(i));
        if (cou == 0) {
          iter_ancsplits_just_int(_root_ratemodel, dists.at(i), leftdists,
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
  } else {
    vector<BranchSegment> *tsegs = node->getSegVector();
    distconds = *tsegs->at(0).distconds;
  }
  if (node->hasParent() == true) {
    vector<BranchSegment> *tsegs = node->getSegVector();
    for (unsigned int i = 0; i < distconds.size(); i++) {
      tsegs->at(0).distconds->at(i) = distconds.at(i);
    }
  } else {
    for (unsigned int i = 0; i < distconds.size(); i++) {
      node->getDoubleVector(_dist_conditionals_key)->at(i) = distconds.at(i);
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
  std::shared_ptr<Node> mrca = _tree->getMRCA(nodeNames);
  auto dists = _root_ratemodel->getDists();
  for (unsigned int i = 0; i < dists.size(); i++) {
    if (dists.at(i).at(fossilarea) == 0) {
      vector<vector<int>> *exd = mrca->getExclDistVector();
      exd->push_back(dists.at(i));
    }
  }
}

void BioGeoTree::setFossilatNodeByMRCA_id(std::shared_ptr<Node> id,
                                          int fossilarea) {
  auto dists = _root_ratemodel->getDists();
  for (unsigned int i = 0; i < dists.size(); i++) {
    if (dists.at(i).at(fossilarea) == 0) {
      vector<vector<int>> *exd = id->getExclDistVector();
      exd->push_back(dists.at(i));
    }
  }
}

void BioGeoTree::setFossilatBranchByMRCA(vector<string> nodeNames,
                                         int fossilarea, double age) {
  std::shared_ptr<Node> mrca = _tree->getMRCA(nodeNames);
  vector<BranchSegment> *tsegs = mrca->getSegVector();
  double startage = mrca->getHeight();
  for (unsigned int i = 0; i < tsegs->size(); i++) {
    if (age > startage && age < (startage + tsegs->at(i).getDuration())) {
      tsegs->at(i).setFossilArea(fossilarea);
    }
    startage += tsegs->at(i).getDuration();
  }
}

void BioGeoTree::setFossilatBranchByMRCA_id(std::shared_ptr<Node> id,
                                            int fossilarea, double age) {
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
void BioGeoTree::prepare_ancstate_reverse() { reverse(_tree->getRoot()); }

/*
 * called from prepare_ancstate_reverse and that is all
 */
void BioGeoTree::reverse(std::shared_ptr<Node> node) {
  _reverse = true;
  vector<Superdouble> *revconds =
      new vector<Superdouble>(_root_ratemodel->getDistsSize(),
                              0); // need to delete this at some point
  if (node == _tree->getRoot()) {
    for (unsigned int i = 0; i < _root_ratemodel->getDistsSize(); i++) {
      revconds->at(i) = 1.0; // prior
    }
    node->assocDoubleVector(_reverse_bits_key, *revconds);
    delete revconds;
    for (int i = 0; i < node->getChildCount(); i++) {
      reverse(node->getChild(i));
    }
  } else {
    vector<Superdouble> *parrev =
        node->getParent()->getDoubleVector(_reverse_bits_key);
    vector<Superdouble> sisdistconds;
    if (node->getParent()->getChild(0) != node) {
      vector<BranchSegment> *tsegs =
          node->getParent()->getChild(0)->getSegVector();
      sisdistconds = tsegs->at(0).alphas;
    } else {
      vector<BranchSegment> *tsegs =
          node->getParent()->getChild(1)->getSegVector();
      sisdistconds = tsegs->at(0).alphas;
    }
    auto dists = _root_ratemodel->getDists();
    vector<int> leftdists;
    vector<int> rightdists;
    double weight;
    vector<Superdouble> tempA(_root_ratemodel->getDistsSize(), 0);
    for (unsigned int i = 0; i < dists.size(); i++) {
      if (accumulate(dists.at(i).begin(), dists.at(i).end(), 0) > 0) {
        vector<vector<int>> *exdist = node->getExclDistVector();
        int cou = count(exdist->begin(), exdist->end(), dists.at(i));
        if (cou == 0) {
          iter_ancsplits_just_int(_root_ratemodel, dists.at(i), leftdists,
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
    vector<BranchSegment> *tsegs = node->getSegVector();
    vector<Superdouble> tempmoveA(tempA);
    for (int ts = tsegs->size() - 1; ts != -1; ts--) {
      for (unsigned int j = 0; j < dists.size(); j++) {
        revconds->at(j) = 0;
      }
      std::shared_ptr<RateModel> rm = tsegs->at(ts).getModel();
      vector<vector<double>> *p =
          &rm->stored_p_matrices[tsegs->at(ts).getPeriod()]
                                [tsegs->at(ts).getDuration()];
      blaze::DynamicMatrix<double> *EN = nullptr;
      blaze::DynamicMatrix<double> *ER = nullptr;
      vector<Superdouble> tempmoveAer(tempA);
      vector<Superdouble> tempmoveAen(tempA);
      if (_stocastic == true) {
        // initialize the segment B's
        for (unsigned int j = 0; j < dists.size(); j++) {
          tempmoveAer[j] = 0;
        }
        for (unsigned int j = 0; j < dists.size(); j++) {
          tempmoveAen[j] = 0;
        }
        EN = &_stored_EN_matrices[tsegs->at(ts).getPeriod()]
                                 [tsegs->at(ts).getDuration()];
        ER = &_stored_ER_matrices[tsegs->at(ts).getPeriod()]
                                 [tsegs->at(ts).getDuration()];
      }
      for (unsigned int j = 0; j < dists.size(); j++) {
        if (accumulate(dists.at(j).begin(), dists.at(j).end(), 0) > 0) {
          for (unsigned int i = 0; i < dists.size(); i++) {
            if (accumulate(dists.at(i).begin(), dists.at(i).end(), 0) > 0) {
              revconds->at(j) +=
                  tempmoveA[i] *
                  ((*p)[i][j]); // tempA needs to change each time
              if (_stocastic == true) {
                tempmoveAer[j] += tempmoveA[i] * (((*ER)(i, j)));
                tempmoveAen[j] += tempmoveA[i] * (((*EN)(i, j)));
              }
            }
          }
        }
      }
      for (unsigned int j = 0; j < dists.size(); j++) {
        tempmoveA[j] = revconds->at(j);
      }
      if (_stocastic == true) {
        tsegs->at(ts).seg_sp_stoch_map_revB_time = tempmoveAer;
        tsegs->at(ts).seg_sp_stoch_map_revB_number = tempmoveAen;
      }
    }
    node->assocDoubleVector(_reverse_bits_key, *revconds);
    delete revconds;
    for (int i = 0; i < node->getChildCount(); i++) {
      reverse(node->getChild(i));
    }
  }
}

/*
 * calculates the most likely split (not state) -- the traditional result for
 * lagrange
 */

unordered_map<vector<int>, vector<AncSplit>>
BioGeoTree::calculate_ancsplit_reverse(std::shared_ptr<Node> node) {
  vector<Superdouble> *Bs = node->getDoubleVector(_reverse_bits_key);
  unordered_map<vector<int>, vector<AncSplit>> ret;
  for (unsigned int j = 0; j < _root_ratemodel->getDistsSize(); j++) {
    vector<int> dist = _root_ratemodel->getDists().at(j);
    vector<AncSplit> ans = iter_ancsplits(_root_ratemodel, dist);
    if (node->isExternal() == false) { // is not a tip
      std::shared_ptr<Node> c1 = node->getChild(0);
      std::shared_ptr<Node> c2 = node->getChild(1);
      vector<BranchSegment> *tsegs1 = c1->getSegVector();
      vector<BranchSegment> *tsegs2 = c2->getSegVector();
      for (unsigned int i = 0; i < ans.size(); i++) {
        vector<vector<int>> *exdist = node->getExclDistVector();
        int cou =
            count(exdist->begin(), exdist->end(),
                  _root_ratemodel->get_int_dists_map().at(ans[i].ancdistint));
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
vector<Superdouble>
BioGeoTree::calculate_ancstate_reverse(std::shared_ptr<Node> node) {
  if (node->isExternal() == false) { // is not a tip
    vector<Superdouble> *Bs = node->getDoubleVector(_reverse_bits_key);
    auto dists = _root_ratemodel->getDists();
    vector<int> leftdists;
    vector<int> rightdists;
    double weight;
    std::shared_ptr<Node> c1 = node->getChild(0);
    std::shared_ptr<Node> c2 = node->getChild(1);
    vector<BranchSegment> *tsegs1 = c1->getSegVector();
    vector<BranchSegment> *tsegs2 = c2->getSegVector();
    vector<Superdouble> v1 = tsegs1->at(0).alphas;
    vector<Superdouble> v2 = tsegs2->at(0).alphas;
    vector<Superdouble> LHOODS(dists.size(), 0);
    for (unsigned int i = 0; i < dists.size(); i++) {
      if (accumulate(dists.at(i).begin(), dists.at(i).end(), 0) > 0) {
        vector<vector<int>> *exdist = node->getExclDistVector();
        int cou = count(exdist->begin(), exdist->end(), dists.at(i));
        if (cou == 0) {
          iter_ancsplits_just_int(_root_ratemodel, dists.at(i), leftdists,
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
  throw std::runtime_error{"calculate_ancstate_reverse was called on a tip"};
}

/**********************************************************
 * forward and reverse stuff for stochastic mapping
 **********************************************************/

void BioGeoTree::prepare_stochmap_reverse_all_nodes(int from, int to) {
  _stocastic = true;
  int ndists = _root_ratemodel->getDistsSize();

  for (unsigned int k = 0; k < _tree->getNodeCount(); k++) {
    vector<BranchSegment> *tsegs = _tree->getNode(k)->getSegVector();
    for (unsigned int l = 0; l < tsegs->size(); l++) {
      int per = (*tsegs)[l].getPeriod();
      double dur = (*tsegs)[l].getDuration();

      lagrange_complex_matrix_t eigvec(ndists, ndists, 0.0);
      lagrange_complex_matrix_t eigval(ndists, ndists, 0.0);

      bool isImag =
          _root_ratemodel->get_eigenvec_eigenval_from_Q(eigval, eigvec, per);

      lagrange_complex_matrix_t Ql(ndists, ndists, 0.0);

      Ql(from, to) = _root_ratemodel->get_Q()[per](from, to);

      lagrange_matrix_t W(ndists, ndists, 0.0);
      W(from, from) = 1;

      lagrange_complex_matrix_t summed(ndists, ndists, 0.0);
      lagrange_complex_matrix_t summedR(ndists, ndists, 0.0);
      for (int i = 0; i < ndists; i++) {
        lagrange_matrix_t Ei(ndists, ndists, 0.0);
        Ei(i, i) = 1;

        lagrange_complex_matrix_t Si(ndists, ndists);
        Si = eigvec * Ei * blaze::inv(eigvec);

        for (int j = 0; j < ndists; j++) {
          std::complex<double> dij = (eigval(i, i) - eigval(j, j)) * dur;
          lagrange_matrix_t Ej(ndists, ndists, 0.0);
          Ej(j, j) = 1;

          lagrange_complex_matrix_t Sj(ndists, ndists);
          Sj = eigvec * Ej * blaze::inv(eigvec);
          std::complex<double> Iijt = 0;
          if (abs(dij) > 10) {
            Iijt = (exp(eigval(i, i) * dur) - exp(eigval(j, j) * dur)) /
                   (eigval(i, i) - eigval(j, j));
          } else if (abs(dij) < 10e-20) {
            Iijt = dur * exp(eigval(j, j) * dur) *
                   (1. + dij / 2. + pow(dij, 2.) / 6. + pow(dij, 3.) / 24.);
          } else {
            if (eigval(i, i) == eigval(j, j)) {
              if (isImag)
                Iijt = dur * exp(eigval(j, j) * dur) * (exp(dij) - 1.) / dij;
              else
                Iijt = dur * exp(eigval(j, j) * dur) * (expm1(real(dij))) / dij;
            } else {
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
      _stored_EN_matrices[per][dur] = (real(summed));
      _stored_EN_CX_matrices[per][dur] = summed;
      _stored_ER_matrices[per][dur] = (real(summedR));
    }
  }
}

/*
 * called directly after reverse_stochastic
 */

vector<Superdouble>
BioGeoTree::calculate_reverse_stochmap(std::shared_ptr<Node> node, bool time) {
  if (node->isExternal() == false) { // is not a tip
    vector<BranchSegment> *tsegs = node->getSegVector();
    auto dists = _root_ratemodel->getDists();
    vector<Superdouble> totalExp(dists.size(), 0);
    for (size_t t = 0; t < tsegs->size(); t++) {
      if (t == 0) {
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_number;
        vector<int> leftdists;
        vector<int> rightdists;
        double weight;
        std::shared_ptr<Node> c1 = node->getChild(0);
        std::shared_ptr<Node> c2 = node->getChild(1);
        vector<BranchSegment> *tsegs1 = c1->getSegVector();
        vector<BranchSegment> *tsegs2 = c2->getSegVector();
        vector<Superdouble> v1 = tsegs1->at(0).alphas;
        vector<Superdouble> v2 = tsegs2->at(0).alphas;
        vector<Superdouble> LHOODS(dists.size(), 0);
        for (unsigned int i = 0; i < dists.size(); i++) {
          if (accumulate(dists.at(i).begin(), dists.at(i).end(), 0) > 0) {
            vector<vector<int>> *exdist = node->getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists.at(i));
            if (cou == 0) {
              iter_ancsplits_just_int(_root_ratemodel, dists.at(i), leftdists,
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
        for (size_t i = 0; i < dists.size(); i++) {
          totalExp[i] = LHOODS[i];
        }
      } else {
        vector<Superdouble> alphs = tsegs->at(t - 1).seg_sp_alphas;
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_number;
        vector<Superdouble> LHOODS(dists.size(), 0);
        for (unsigned int i = 0; i < dists.size(); i++) {
          if (accumulate(dists.at(i).begin(), dists.at(i).end(), 0) > 0) {
            vector<vector<int>> *exdist = node->getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists.at(i));
            if (cou == 0) {
              LHOODS[i] =
                  Bs.at(i) * (alphs[i]); // do i do this or do i do from i to j
            }
          }
        }
        for (size_t i = 0; i < dists.size(); i++) {
          totalExp[i] += LHOODS[i];
        }
      }
    }
    // not sure if this should return a Superdouble or not when doing a bigtree
    return totalExp;
  } else {
    vector<BranchSegment> *tsegs = node->getSegVector();
    auto dists = _root_ratemodel->getDists();
    vector<Superdouble> totalExp(dists.size(), 0);
    for (size_t t = 0; t < tsegs->size(); t++) {
      if (t == 0) {
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_number;
        vector<Superdouble> LHOODS(dists.size(), 0);
        for (unsigned int i = 0; i < dists.size(); i++) {
          if (accumulate(dists.at(i).begin(), dists.at(i).end(), 0) > 0) {
            vector<vector<int>> *exdist = node->getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists.at(i));
            if (cou == 0) {
              LHOODS[i] = Bs.at(i) * (tsegs->at(0).distconds->at(i));
            }
          }
        }
        for (size_t i = 0; i < dists.size(); i++) {
          totalExp[i] = LHOODS[i];
        }
      } else {
        vector<Superdouble> alphs = tsegs->at(t - 1).seg_sp_alphas;
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs->at(t).seg_sp_stoch_map_revB_number;
        vector<Superdouble> LHOODS(dists.size(), 0);
        for (unsigned int i = 0; i < dists.size(); i++) {
          if (accumulate(dists.at(i).begin(), dists.at(i).end(), 0) > 0) {
            vector<vector<int>> *exdist = node->getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists.at(i));
            if (cou == 0) {
              LHOODS[i] = Bs.at(i) * (alphs[i]);
            }
          }
        }
        for (size_t i = 0; i < dists.size(); i++) {
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
  for (unsigned int i = 0; i < _tree->getNodeCount(); i++) {
    vector<BranchSegment> *tsegs = _tree->getNode(i)->getSegVector();
    for (unsigned int j = 0; j < tsegs->size(); j++) {
      delete tsegs->at(j).distconds;
      delete tsegs->at(j).ancdistconds;
    }
    _tree->getNode(i)->deleteExclDistVector();
    if (_reverse == true && _tree->getNode(i)->isInternal()) {
      _tree->getNode(i)->deleteDoubleVector(_reverse_bits_key);
    }
    _tree->getNode(i)->deleteSegVector();
  }
  _tree->getRoot()->deleteDoubleVector(_dist_conditionals_key);
  _tree->getRoot()->deleteDoubleVector(_anc_dist_conditionals_key);
  _tree->getRoot()->deleteDoubleVector(_reverse_bits_key);
}

size_t hash<std::vector<int>>::operator()(const std::vector<int> &vec) const {
  std::size_t seed = vec.size();
  for (auto &i : vec) {
    seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}
