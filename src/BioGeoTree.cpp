/*
 * BioGeoTree.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */
#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <numeric>

using namespace std;

#include "BioGeoTree.h"
#include "BioGeoTreeTools.h"
#include "RateMatrixUtils.h"
#include "Utils.h"

Superdouble MAX(const Superdouble &a, const Superdouble &b) {
  return b > a ? b : a;
}

/*
 * sloppy beginning but best for now because of the complicated bits
 */

BioGeoTree::BioGeoTree(
    std::shared_ptr<Tree> tr, const vector<double> &ps, size_t regions,
    const std::unordered_map<std::string, lagrange_dist_t> &distrib_data)
    : _tree(tr),
      _periods(ps),
      _columns(nullptr),
      _which_columns(nullptr),
      _root_ratemodel(nullptr),
      _store_p_matrices(false),
      _use_stored_matrices(false),
      _reverse(false),
      _stocastic(false),
      _stored_EN_matrices{},
      _stored_EN_CX_matrices{},
      _stored_ER_matrices{} {
  if (_periods.size() == 0) {
    throw std::runtime_error{"No periods when creating a biogeotree"};
  }

  for (size_t i = 0; i < _tree->getInternalNodeCount(); ++i) {
    auto current_node = _tree->getInternalNode(i);
    for (auto child : current_node->getChildren()) {
      double ancestor_height = current_node->getHeight();
      double descendant_height = child->getHeight();
      double t = descendant_height;

      double s = 0;
      for (size_t j = 0; j < _periods.size(); j++) {
        s += _periods[j];
        if (t < s) {
          double duration = min(s - t, ancestor_height - t);
          if (duration > 0) {
            BranchSegment tseg = BranchSegment(duration, j);
            child->getSegVector().push_back(tseg);
          }
          t += duration;  // TODO: make sure that this is all working
        }
        if (t > ancestor_height || _periods[j] > t) {
          break;
        }
      }
    }
  }
}

void BioGeoTree::set_store_p_matrices(bool i) { _store_p_matrices = i; }

void BioGeoTree::set_use_stored_matrices(bool i) { _use_stored_matrices = i; }

void BioGeoTree::set_default_model(std::shared_ptr<RateModel> mod) {
  _root_ratemodel = mod;
  for (unsigned int i = 0; i < _tree->getNodeCount(); i++) {
    vector<BranchSegment> &tsegs = _tree->getNode(i)->getSegVector();
    for (unsigned int j = 0; j < tsegs.size(); j++) {
      tsegs[j].distconds =
          make_shared<vector<Superdouble>>(_root_ratemodel->getDistsSize(), 0);
      tsegs[j].ancdistconds =
          make_shared<vector<Superdouble>>(_root_ratemodel->getDistsSize(), 0);
    }
  }
}

void BioGeoTree::update_default_model(std::shared_ptr<RateModel> mod) {
  _root_ratemodel = mod;
}

void BioGeoTree::set_tip_conditionals(
    unordered_map<string, uint64_t> distrib_data) {
  int numofleaves = _tree->getExternalNodeCount();
  for (int i = 0; i < numofleaves; i++) {
    vector<BranchSegment> &tsegs = _tree->getExternalNode(i)->getSegVector();
    auto ratemodel_dists = _root_ratemodel->getDists();
    uint64_t key = distrib_data[_tree->getExternalNode(i)->getName()];
    size_t index = 0;
    for (index = 0; index < ratemodel_dists.size(); index++) {
      if (ratemodel_dists[index] == key) {
        break;
      }
    }
    tsegs[0].distconds->at(index) = 1.0;
  }
}

void BioGeoTree::set_excluded_dist(lagrange_dist_t ind,
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
  // This line is wrong. There should be a prior frequency for the CLVS. This is
  //"fine" because its uniform (this is to say that it is only a linear
  // transformation), but the LH values reported by this method are wrong.
  return -(calculate_vector_Superdouble_sum(
               _tree->getRoot()->getConditionalVector()))
              .getLn();
}

vector<Superdouble> BioGeoTree::conditionals(std::shared_ptr<Node> node,
                                             bool marginal, bool sparse) {
  vector<Superdouble> distconds;
  vector<BranchSegment> &time_segments = node->getSegVector();

  distconds = *time_segments[0].distconds;

  for (unsigned int i = 0; i < time_segments.size(); i++) {
    for (unsigned int j = 0; j < distconds.size(); j++) {
      time_segments[i].distconds->at(j) = distconds.at(j);
    }

    vector<Superdouble> v(_root_ratemodel->getDistsSize(), 0);
    vector<int> distrange;

    if (time_segments[i].get_start_dist_int() != -666) {
      int ind1 = time_segments[i].get_start_dist_int();
      distrange.push_back(ind1);
    } else if (time_segments[i].getFossilAreas().size() > 0) {
      for (unsigned int j = 0; j < _root_ratemodel->getDistsSize(); j++) {
        distrange.push_back(j);
      }
      for (unsigned int k = 0; k < distrange.size(); k++) {
        bool flag = true;
        for (unsigned int x = 0; x < time_segments[i].getFossilAreas().size();
             x++) {
          if (time_segments[i].getFossilAreas()[x] == 1 &&
              distrange.at(x) == 0) {
            flag = false;
            break;
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
          p = _root_ratemodel->setup_fortran_P(time_segments[i].getPeriod(),
                                               time_segments[i].getDuration(),
                                               _store_p_matrices);
        } else {
          p = _root_ratemodel
                  ->stored_p_matrices[time_segments[i].getPeriod()]
                                     [time_segments[i].getDuration()];
        }
        for (unsigned int j = 0; j < distrange.size(); j++) {
          for (unsigned int k = 0; k < distconds.size(); k++) {
            v[distrange[j]] += (distconds.at(k) * p[distrange[j]][k]);
          }
        }
      } else {  // sparse
        /*
          testing pthread version
        */
        if (_root_ratemodel->get_nthreads() > 0) {
        } else {
          for (unsigned int j = 0; j < distrange.size(); j++) {
            bool inthere = false;
            if (_columns->at(distrange[j]) == 1) inthere = true;
            vector<double> p;
            if (inthere == true) {
              p = _root_ratemodel->setup_sparse_single_column_P(
                  time_segments[i].getPeriod(), time_segments[i].getDuration(),
                  distrange[j]);
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
      time_segments[i].seg_sp_alphas = distconds;
    }
  }
  /*
   * if store is true we want to store the conditionals for each node
   * for possible use in ancestral state reconstruction
   */
  if (_store_p_matrices == true) {
    time_segments[0].alphas = distconds;
  }
  return distconds;
}

void BioGeoTree::ancdist_conditional_lh(std::shared_ptr<Node> node,
                                        bool marginal) {
  vector<Superdouble> distconds(_root_ratemodel->getDistsSize(), 0);
  if (node->isExternal() == false) {  // is not a tip
    std::shared_ptr<Node> c1 = node->getChild(0);
    std::shared_ptr<Node> c2 = node->getChild(1);
    ancdist_conditional_lh(c1, marginal);
    ancdist_conditional_lh(c2, marginal);
    bool sparse = _root_ratemodel->_sparse;
    vector<Superdouble> v1;
    vector<Superdouble> v2;
    if (sparse == true) {
      // getcolumns
      vector<BranchSegment> &c1tsegs = c1->getSegVector();
      vector<BranchSegment> &c2tsegs = c2->getSegVector();
      vector<int> lcols =
          _root_ratemodel->get_columns_for_sparse(*c1tsegs[0].distconds);
      vector<int> rcols =
          _root_ratemodel->get_columns_for_sparse(*c2tsegs[0].distconds);
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
    vector<int> leftdists_index;
    vector<int> rightdists_index;
    double weight;

    for (unsigned int i = 0; i < dists.size(); i++) {
      if (lagrange_popcount(dists[i]) > 0) {
        Superdouble lh = 0.0;
        auto exdist = node->getExclDistVector();
        auto iter = std::find(exdist->begin(), exdist->end(), dists[i]);
        if (iter == exdist->end()) {
          _root_ratemodel->iter_ancsplits_just_int(dists.at(i), leftdists_index,
                                                   rightdists_index, weight);
          for (unsigned int j = 0; j < leftdists_index.size(); j++) {
            int ind1 = leftdists_index.at(j);
            int ind2 = rightdists_index.at(j);
            Superdouble lh_part = v1.at(ind1) * v2.at(ind2);
            lh += (lh_part * weight);
          }
        }
        distconds.at(i) = lh;
      }
    }
  } else {
    vector<BranchSegment> &tsegs = node->getSegVector();
    distconds = *tsegs[0].distconds;
  }
  if (_tree->getParent(node) != nullptr) {
    vector<BranchSegment> &tsegs = node->getSegVector();
    for (unsigned int i = 0; i < distconds.size(); i++) {
      tsegs[0].distconds->at(i) = distconds.at(i);
    }
  } else {
    node->setConditionalVector(distconds);
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
    if (lagrange_bextr(dists.at(i), fossilarea) == 0) {
      auto exd = mrca->getExclDistVector();
      exd->push_back(dists.at(i));
    }
  }
}

void BioGeoTree::setFossilatNodeByMRCA_id(std::shared_ptr<Node> id,
                                          int fossilarea) {
  auto dists = _root_ratemodel->getDists();
  auto exd = id->getExclDistVector();
  for (unsigned int i = 0; i < dists.size(); i++) {
    if (lagrange_bextr(dists.at(i), fossilarea) == 0) {
      exd->push_back(dists.at(i));
    }
  }
}

void BioGeoTree::setFossilatBranchByMRCA(vector<string> nodeNames,
                                         int fossilarea, double age) {
  std::shared_ptr<Node> mrca = _tree->getMRCA(nodeNames);
  vector<BranchSegment> &tsegs = mrca->getSegVector();
  double startage = mrca->getHeight();
  for (unsigned int i = 0; i < tsegs.size(); i++) {
    if (age > startage && age < (startage + tsegs[i].getDuration())) {
      tsegs[i].setFossilArea(fossilarea);
    }
    startage += tsegs[i].getDuration();
  }
}

void BioGeoTree::setFossilatBranchByMRCA_id(std::shared_ptr<Node> id,
                                            int fossilarea, double age) {
  vector<BranchSegment> &tsegs = id->getSegVector();
  double startage = id->getHeight();

  for (unsigned int i = 0; i < tsegs.size(); i++) {
    if (age > startage && age < (startage + tsegs[i].getDuration())) {
      tsegs[i].setFossilArea(fossilarea);
    }
    startage += tsegs[i].getDuration();
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
  vector<Superdouble> revconds =
      vector<Superdouble>(_root_ratemodel->getDistsSize(), 0);
  if (node == _tree->getRoot()) {
    for (unsigned int i = 0; i < _root_ratemodel->getDistsSize(); i++) {
      revconds[i] = 1.0;  // prior
    }
    node->setReverseBits(revconds);
    for (int i = 0; i < node->getChildCount(); i++) {
      reverse(node->getChild(i));
    }
  } else {
    auto parrev = _tree->getParent(node)->getReverseBits();
    vector<Superdouble> sisdistconds;
    if (_tree->getParent(node)->getChild(0) != node) {
      vector<BranchSegment> &tsegs =
          _tree->getParent(node)->getChild(0)->getSegVector();
      sisdistconds = tsegs[0].alphas;
    } else {
      vector<BranchSegment> &tsegs =
          _tree->getParent(node)->getChild(1)->getSegVector();
      sisdistconds = tsegs[0].alphas;
    }
    auto dists = _root_ratemodel->getDists();
    vector<int> leftdists;
    vector<int> rightdists;
    double weight;
    vector<Superdouble> tempA(_root_ratemodel->getDistsSize(), 0);
    for (unsigned int i = 0; i < dists.size(); i++) {
      if (lagrange_popcount(dists[i]) > 0) {
        auto exdist = node->getExclDistVector();
        int cou = count(exdist->begin(), exdist->end(), dists.at(i));
        if (cou == 0) {
          _root_ratemodel->iter_ancsplits_just_int(dists.at(i), leftdists,
                                                   rightdists, weight);
          // root has i, curnode has left, sister of cur has right
          for (unsigned int j = 0; j < leftdists.size(); j++) {
            int ind1 = leftdists[j];
            int ind2 = rightdists[j];
            tempA[ind1] += (sisdistconds.at(ind2) * weight * parrev[i]);
          }
        }
      }
    }

    // now calculate node B
    vector<BranchSegment> &tsegs = node->getSegVector();
    vector<Superdouble> tempmoveA(tempA);
    for (int ts = tsegs.size() - 1; ts != -1; ts--) {
      for (unsigned int j = 0; j < dists.size(); j++) {
        revconds[j] = 0;
      }
      vector<vector<double>> *p =
          &_root_ratemodel->stored_p_matrices[tsegs[ts].getPeriod()]
                                             [tsegs[ts].getDuration()];
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
        EN = &_stored_EN_matrices[tsegs[ts].getPeriod()]
                                 [tsegs[ts].getDuration()];
        ER = &_stored_ER_matrices[tsegs[ts].getPeriod()]
                                 [tsegs[ts].getDuration()];
      }
      for (unsigned int j = 0; j < dists.size(); j++) {
        if (lagrange_popcount(dists[j]) == 0) {
          continue;
        }

        for (unsigned int i = 0; i < dists.size(); i++) {
          if (lagrange_popcount(dists[i]) == 0) {
            continue;
          }
          revconds[j] +=
              tempmoveA[i] * ((*p)[i][j]);  // tempA needs to change each time
          if (_stocastic == true) {
            tempmoveAer[j] += tempmoveA[i] * (((*ER)(i, j)));
            tempmoveAen[j] += tempmoveA[i] * (((*EN)(i, j)));
          }
        }
      }
      for (unsigned int j = 0; j < dists.size(); j++) {
        tempmoveA[j] = revconds[j];
      }
      if (_stocastic == true) {
        tsegs[ts].seg_sp_stoch_map_revB_time = tempmoveAer;
        tsegs[ts].seg_sp_stoch_map_revB_number = tempmoveAen;
      }
    }
    node->setReverseBits(revconds);
    for (int i = 0; i < node->getChildCount(); i++) {
      reverse(node->getChild(i));
    }
  }
}

/*
 * calculates the most likely split (not state) -- the traditional result for
 * lagrange
 */

unordered_map<lagrange_dist_t, vector<AncSplit>>
BioGeoTree::calculate_ancsplit_reverse(std::shared_ptr<Node> node) {
  auto Bs = node->getReverseBits();
  unordered_map<lagrange_dist_t, vector<AncSplit>> ret;
  for (unsigned int j = 0; j < _root_ratemodel->getDistsSize(); j++) {
    lagrange_dist_t dist = _root_ratemodel->getDists().at(j);
    vector<AncSplit> ans = _root_ratemodel->iter_ancsplits(dist);
    if (node->isExternal() == false) {  // is not a tip
      std::shared_ptr<Node> c1 = node->getChild(0);
      std::shared_ptr<Node> c2 = node->getChild(1);
      vector<BranchSegment> &tsegs1 = c1->getSegVector();
      vector<BranchSegment> &tsegs2 = c2->getSegVector();
      for (unsigned int i = 0; i < ans.size(); i++) {
        auto exdist = node->getExclDistVector();
        int cou =
            count(exdist->begin(), exdist->end(),
                  _root_ratemodel->get_int_dists_map().at(ans[i].ancdistint));
        if (cou == 0) {
          vector<Superdouble> v1 = tsegs1[0].alphas;
          vector<Superdouble> v2 = tsegs2[0].alphas;
          Superdouble lh = (v1[ans[i].ldescdistint] * v2[ans[i].rdescdistint] *
                            Bs[j] * ans[i].getWeight());
          ans[i].setLikelihood(lh);
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
vector<Superdouble> BioGeoTree::calculate_ancstate_reverse(
    std::shared_ptr<Node> node) {
  if (node->isExternal() == false) {  // is not a tip
    auto Bs = node->getReverseBits();
    auto dists = _root_ratemodel->getDists();
    vector<int> leftdists;
    vector<int> rightdists;
    double weight;
    std::shared_ptr<Node> c1 = node->getChild(0);
    std::shared_ptr<Node> c2 = node->getChild(1);
    vector<BranchSegment> &tsegs1 = c1->getSegVector();
    vector<BranchSegment> &tsegs2 = c2->getSegVector();
    vector<Superdouble> v1 = tsegs1[0].alphas;
    vector<Superdouble> v2 = tsegs2[0].alphas;
    vector<Superdouble> LHOODS(dists.size(), 0);
    for (unsigned int i = 0; i < dists.size(); i++) {
      if (lagrange_popcount(dists[i]) > 0) {
        auto exdist = node->getExclDistVector();
        int cou = count(exdist->begin(), exdist->end(), dists.at(i));
        if (cou == 0) {
          _root_ratemodel->iter_ancsplits_just_int(dists.at(i), leftdists,
                                                   rightdists, weight);
          for (unsigned int j = 0; j < leftdists.size(); j++) {
            int ind1 = leftdists[j];
            int ind2 = rightdists[j];
            LHOODS[i] += (v1.at(ind1) * v2.at(ind2) * weight);
          }
          LHOODS[i] *= Bs[i];
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
    vector<BranchSegment> &tsegs = _tree->getNode(k)->getSegVector();
    for (unsigned int l = 0; l < tsegs.size(); l++) {
      int per = tsegs[l].getPeriod();
      double dur = tsegs[l].getDuration();

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

vector<Superdouble> BioGeoTree::calculate_reverse_stochmap(
    std::shared_ptr<Node> node, bool time) {
  if (node->isExternal() == false) {  // is not a tip
    vector<BranchSegment> &tsegs = node->getSegVector();
    auto dists = _root_ratemodel->getDists();
    vector<Superdouble> totalExp(dists.size(), 0);
    for (size_t t = 0; t < tsegs.size(); t++) {
      if (t == 0) {
        vector<Superdouble> Bs;
        if (time) {
          Bs = tsegs[t].seg_sp_stoch_map_revB_time;
        } else {
          Bs = tsegs[t].seg_sp_stoch_map_revB_number;
        }
        vector<int> leftdists;
        vector<int> rightdists;
        double weight;
        std::shared_ptr<Node> c1 = node->getChild(0);
        std::shared_ptr<Node> c2 = node->getChild(1);
        vector<BranchSegment> &tsegs1 = c1->getSegVector();
        vector<BranchSegment> &tsegs2 = c2->getSegVector();
        vector<Superdouble> v1 = tsegs1[0].alphas;
        vector<Superdouble> v2 = tsegs2[0].alphas;
        vector<Superdouble> LHOODS(dists.size(), 0);
        for (unsigned int i = 0; i < dists.size(); i++) {
          if (lagrange_popcount(dists[i]) > 0) {
            auto exdist = node->getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists.at(i));
            if (cou == 0) {
              continue;
            }
            _root_ratemodel->iter_ancsplits_just_int(dists.at(i), leftdists,
                                                     rightdists, weight);
            for (unsigned int j = 0; j < leftdists.size(); j++) {
              int ind1 = leftdists[j];
              int ind2 = rightdists[j];
              LHOODS[i] += (v1.at(ind1) * v2.at(ind2) * weight);
            }
            LHOODS[i] *= Bs.at(i);
          }
        }
        for (size_t i = 0; i < dists.size(); i++) {
          totalExp[i] = LHOODS[i];
        }
      } else {
        vector<Superdouble> alphs = tsegs[t - 1].seg_sp_alphas;
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs[t].seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs[t].seg_sp_stoch_map_revB_number;
        vector<Superdouble> LHOODS(dists.size(), 0);
        for (unsigned int i = 0; i < dists.size(); i++) {
          if (lagrange_popcount(dists[i]) > 0) {
            auto exdist = node->getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists.at(i));
            if (cou == 0) {
              LHOODS[i] =
                  Bs.at(i) * (alphs[i]);  // do i do this or do i do from i to j
            }
          }
        }
        for (size_t i = 0; i < dists.size(); i++) {
          totalExp[i] += LHOODS[i];
        }
      }
    }
    // not sure if this should return a Superdouble or not when doing a
    // bigtree
    return totalExp;
  } else {
    vector<BranchSegment> &tsegs = node->getSegVector();
    auto dists = _root_ratemodel->getDists();
    vector<Superdouble> totalExp(dists.size(), 0);
    for (size_t t = 0; t < tsegs.size(); t++) {
      if (t == 0) {
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs[t].seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs[t].seg_sp_stoch_map_revB_number;
        vector<Superdouble> LHOODS(dists.size(), 0);
        for (unsigned int i = 0; i < dists.size(); i++) {
          if (lagrange_popcount(dists[i]) > 0) {
            auto exdist = node->getExclDistVector();
            int cou = count(exdist->begin(), exdist->end(), dists.at(i));
            if (cou == 0) {
              LHOODS[i] = Bs.at(i) * (tsegs[0].distconds->at(i));
            }
          }
        }
        for (size_t i = 0; i < dists.size(); i++) {
          totalExp[i] = LHOODS[i];
        }
      } else {
        vector<Superdouble> alphs = tsegs[t - 1].seg_sp_alphas;
        vector<Superdouble> Bs;
        if (time)
          Bs = tsegs[t].seg_sp_stoch_map_revB_time;
        else
          Bs = tsegs[t].seg_sp_stoch_map_revB_number;
        vector<Superdouble> LHOODS(dists.size(), 0);
        for (unsigned int i = 0; i < dists.size(); i++) {
          if (lagrange_popcount(dists[i]) > 0) {
            auto exdist = node->getExclDistVector();
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

size_t hash<std::vector<int>>::operator()(const std::vector<int> &vec) const {
  std::size_t seed = vec.size();
  for (auto &i : vec) {
    seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}
