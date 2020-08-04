/*
 * OptimizeBioGeo.h
 *
 *  Created on: Aug 18, 2009
 *      Author: smitty
 */

#ifndef OPTIMIZEBIOGEOALLDISPERSAL_NLOPT_H_
#define OPTIMIZEBIOGEOALLDISPERSAL_NLOPT_H_

#include <vector>
using namespace std;

#include "BioGeoTree.h"
#include "RateModel.h"

double get_likelihood_with_optimized_dispersal_extinction(unsigned n,
                                                          const double *x,
                                                          double *g,
                                                          void *state);
vector<double> optimize_dispersal_extinction_all_nlopt(BioGeoTree *init_tree,
                                                       std::shared_ptr<RateModel> init_rm);

#endif /* OPTIMIZEBIOGEOALLDISPERSAL_H_ */
