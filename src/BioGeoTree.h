/*
 * BioGeoTree.h
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */

#ifndef BIOGEOTREE_H_
#define BIOGEOTREE_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

#include "AncSplit.h"
#include "BranchSegment.h"
#include "RateModel.h"
#include "node.h"
#include "tree.h"
#include "vector_node_object.h"

class BioGeoTree {
private:
  std::shared_ptr<Tree> _tree;
  vector<double> _periods;
  const string _dist_conditionals_key = "dist_conditionals";
  const string _anc_dist_conditionals_key = "anc_dist_conditionals";
  std::shared_ptr<vector<int>> _columns;
  std::shared_ptr<vector<int>> _which_columns;
  std::shared_ptr<RateModel> _root_ratemodel;
  bool _store_p_matrices;
  bool _use_stored_matrices;

  // reverse bits
  const string _reverse_bits_key = "revB";
  bool _reverse;
  // end reverse bits

  // stochastic mapping bits
  bool _stocastic;
  // map of period int and then branch length Superdouble
  unordered_map<int, map<double, blaze::DynamicMatrix<double>>>
      _stored_EN_matrices;
  unordered_map<int, map<double, blaze::DynamicMatrix<std::complex<double>>>>
      _stored_EN_CX_matrices;
  unordered_map<int, map<double, blaze::DynamicMatrix<double>>>
      _stored_ER_matrices;
  // end mapping bits

public:
  BioGeoTree(std::shared_ptr<Tree> tr, vector<double> ps);
  void set_store_p_matrices(bool);
  void set_use_stored_matrices(bool);
  void set_default_model(std::shared_ptr<RateModel> mod);
  void update_default_model(std::shared_ptr<RateModel> mod);
  Superdouble eval_likelihood(bool marg);
  void set_excluded_dist(vector<int> ind, std::shared_ptr<Node> node);
  void set_tip_conditionals(unordered_map<string, vector<int>> distrib_data);
  vector<Superdouble> conditionals(std::shared_ptr<Node> node, bool marg,
                                   bool sparse);
  // void ancdist_conditional_lh(bpp::Node & node, bool marg);
  void ancdist_conditional_lh(std::shared_ptr<Node> node, bool marg);

  /*
          fossil data
   */
  void setFossilatNodeByMRCA(vector<string> nodeNames, int fossilarea);
  void setFossilatNodeByMRCA_id(std::shared_ptr<Node> id, int fossilarea);
  void setFossilatBranchByMRCA(vector<string> nodeNames, int fossilarea,
                               double age);
  void setFossilatBranchByMRCA_id(std::shared_ptr<Node> id, int fossilarea,
                                  double age);

  /*
          for calculating forward and reverse
   */
  void prepare_ancstate_reverse();
  void reverse(std::shared_ptr<Node>);
  unordered_map<vector<int>, vector<AncSplit>>
  calculate_ancsplit_reverse(std::shared_ptr<Node> node);
  vector<Superdouble> calculate_ancstate_reverse(std::shared_ptr<Node> node);
  ~BioGeoTree();
  // need to override these at some point
  BioGeoTree(const BioGeoTree &L); // copy constructor
  BioGeoTree &operator=(const BioGeoTree &L);

  /*
   * for calculating forward and reverse for expected values (stochastic
   * mapping)
   */
  void prepare_stochmap_reverse_all_nodes(int, int);
  vector<Superdouble> calculate_reverse_stochmap(std::shared_ptr<Node>, bool);
  vector<Superdouble>
  calculate_reverse_stochmap_TEST(std::shared_ptr<Node> node, bool time);

  /*
          for timing things
   */
  Superdouble ti;
  Superdouble ti2;
  Superdouble ti3;
};

#endif /* BIOGEOTREE_H_ */
