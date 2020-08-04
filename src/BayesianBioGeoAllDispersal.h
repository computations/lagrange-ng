/*
 *  BayesianBioGeo.h
 *  lagrange_cpp
 *
 *  Created by Stephen Smith on 1/13/10.
 *  Copyright 2010 Yale University. All rights reserved.
 *
 */

#ifndef BAYESIANBIOGEOALLDISPERSAL_H_
#define BAYESIANBIOGEOALLDISPERSAL_H_

#include <vector>
using namespace std;

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "BioGeoTree.h"
#include "RateModel.h"

class BayesianBioGeoAllDispersal {
private:
  BioGeoTree *_tree;
  std::shared_ptr<RateModel> _rate_model;
  int _generations;
  bool _marginal;
  vector<double> _params;
  vector<double> _prev_params;
  vector<vector<vector<double>>> _dispersal_mask;
  const gsl_rng_type *_rng_type;
  gsl_rng *_rng;
  double calculate_pdf(double value);
  double calculate_sliding(double value, double sliding);
  double calculate_sliding_log(double value, double sliding, double *hastings);

public:
  BayesianBioGeoAllDispersal(BioGeoTree *intree,
                             std::shared_ptr<RateModel> inrm, bool marg,
                             int gen);
  void run_global_dispersal_extinction();
};

#endif /* BAYESIANBIOGEO_H_ */
