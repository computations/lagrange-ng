/*
 *  BayesianBioGeo.h
 *  lagrange_cpp
 *
 *  Created by Stephen Smith on 1/13/10.
 *  Copyright 2010 Yale University. All rights reserved.
 *
 */

#ifndef BAYESIANBIOGEOSTMAPPER_H_
#define BAYESIANBIOGEOSTMAPPER_H_

#include <map>
#include <string>
#include <vector>
using namespace std;

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "BioGeoTree.h"
#include "RateModel.h"
#include "node.h"
#include "tree.h"

class BayesianBioGeoStMapper {
private:
  BioGeoTree *bgt;
  Tree *tree;
  RateModel *rm;
  int gens;
  bool marginal;
  vector<double> params;
  const gsl_rng_type *T;
  gsl_rng *r;
  string sts;
  string pts;

  map<Node *, int> startstate_node_map;
  map<Node *, int> endstate_node_map;

  void initialize_nodes();
  void clean_nodes();
  bool mapping(double treelen, double totlen);
  void sample_ancstates();
  bool simulate_on_nodes();
  bool simulate_on_branch(int starting_state, int finished_state, double brlen,
                          vector<int> *newstates, vector<double> *newpoints);
  int draw_new_state(int starting_state);

public:
  BayesianBioGeoStMapper(BioGeoTree *inbgt, Tree *intr, RateModel *inrm,
                         bool marg, int gen);
  void run_mappings();
};

#endif /* BAYESIANBIOGEO_H_ */
