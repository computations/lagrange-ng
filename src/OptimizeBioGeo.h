/*
 * OptimizeBioGeo.h
 *
 *  Created on: Aug 18, 2009
 *      Author: smitty
 */

#ifndef OPTIMIZEBIOGEO_H_
#define OPTIMIZEBIOGEO_H_

#include <vector>
using namespace std;

#include "BioGeoTree.h"
#include "RateModel.h"

#include <gsl/gsl_vector.h>

class OptimizeBioGeo {
private:
  BioGeoTree *_tree;
  RateModel *_rate_model;
  const size_t _max_iterations = 100;
  const double _abs_tol = 0.0001;
  bool _marginal;
  double
  GetLikelihoodWithOptimizedDispersalExtinction(const gsl_vector *variables);
  static double
  GetLikelihoodWithOptimizedDispersalExtinction_gsl(const gsl_vector *variables,
                                                    void *obj);

public:
  OptimizeBioGeo(BioGeoTree *intree, RateModel *inrm, bool marg);
  vector<double> optimize_global_dispersal_extinction();
};

#endif /* OPTIMIZEBIOGEO_H_ */
