/*
 *  BayesianBioGeo.h
 *  lagrange_cpp
 *
 *  Created by Stephen Smith on 1/13/10.
 *  Copyright 2010 Yale University. All rights reserved.
 *
 */

#ifndef BAYESIANBIOGEO_H_
#define BAYESIANBIOGEO_H_

#include <vector>
using namespace std;

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "BioGeoTree.h"
#include "RateModel.h"

class BayesianBioGeo {
private:
  BioGeoTree *_tree;
  RateModel *_rate_model;
  int _generations;
  bool _marginal;
  vector<double> _params;
  vector<double> _prev_params;
  const gsl_rng_type *_rng_type;
  gsl_rng *_rng;
  double calculate_pdf(double value);
  double calculate_sliding(double value, double sliding);
  double calculate_sliding_log(double value, double sliding, double *hastings);

public:
  BayesianBioGeo(BioGeoTree *intree, RateModel *inrm, bool marg, int gen);
  void run_global_dispersal_extinction();
};

#endif /* BAYESIANBIOGEO_H_ */
