/*
 * OptimizeBioGeo.h
 *
 *  Created on: Aug 18, 2009
 *      Author: smitty
 */

#ifndef OPTIMIZEBIOGEOALLDISPERSAL_H_
#define OPTIMIZEBIOGEOALLDISPERSAL_H_

#include <vector>
using namespace std;

#include "BioGeoTree.h"
#include "RateModel.h"

#include <gsl/gsl_vector.h>

class OptimizeBioGeoAllDispersal {
private:
  BioGeoTree *tree;
  RateModel *rm;
  vector<vector<vector<double>>> D_mask;
  int maxiterations;
  double stoppingprecision;
  bool marginal;
  int nareas;
  double
  GetLikelihoodWithOptimizedDispersalExtinction(const gsl_vector *variables);
  static double
  GetLikelihoodWithOptimizedDispersalExtinction_gsl(const gsl_vector *variables,
                                                    void *obj);

public:
  OptimizeBioGeoAllDispersal(BioGeoTree *intree, RateModel *inrm, bool marg);
  vector<double> optimize_global_dispersal_extinction();
};

#endif /* OPTIMIZEBIOGEOALLDISPERSAL_H_ */
