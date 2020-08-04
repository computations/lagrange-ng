/*
 * OptimizeBioGeo.cpp
 *
 *  Created on: Aug 18, 2009
 *      Author: smitty
 */

#include <limits>
#include <math.h>
#include <vector>
using namespace std;

#include "BioGeoTree.h"
#include "OptimizeBioGeoAllDispersal_nlopt.h"
#include "RateModel.h"

#include <nlopt.h>

BioGeoTree *nlopt_intree;
std::shared_ptr<RateModel> nlopt_inrm;
vector<vector<vector<double>>> nlopt_D_mask;

double get_likelihood_with_optimized_dispersal_extinction(unsigned n,
                                                          const double *x,
                                                          double *g,
                                                          void *state) {
  (void)n;
  (void)g;
  (void)state;
  double dispersal = x[0];
  double extinction = x[1];
  int count = 2;
  for (unsigned int i = 0; i < nlopt_D_mask.size(); i++) {
    for (unsigned int j = 0; j < nlopt_D_mask[i].size(); j++) {
      nlopt_D_mask[i][j][j] = 1.0;
      for (unsigned int k = 0; k < nlopt_D_mask[i][j].size(); k++) {
        if (k != j) {
          nlopt_D_mask[i][j][k] = x[count];
          if (nlopt_D_mask[i][j][k] <= 0) {
            return std::numeric_limits<double>::max();
          }
          count += 1;
        }
      }
    }
  }
  double like;
  if (dispersal <= 0 || extinction <= 0)
    return std::numeric_limits<double>::max();
  if (dispersal > 100 || extinction > 100)
    return std::numeric_limits<double>::max();
  nlopt_inrm->setup_D_provided(dispersal, nlopt_D_mask);
  nlopt_inrm->setup_E(extinction);
  nlopt_inrm->setup_Q();
  nlopt_intree->update_default_model(nlopt_inrm);
  like = nlopt_intree->eval_likelihood(true);
  if (like < 0 || like == std::numeric_limits<double>::infinity())
    like = std::numeric_limits<double>::max();
  return like;
}

vector<double>
optimize_dispersal_extinction_all_nlopt(BioGeoTree *init_tree,
                                        std::shared_ptr<RateModel> init_rm) {
  nlopt_intree = init_tree;
  nlopt_inrm = init_rm;
  int nareas = init_rm->get_num_areas();
  int nperiods = init_rm->get_num_periods();
  vector<double> cols(nareas, 1);
  vector<vector<double>> rows(nareas, cols);
  nlopt_D_mask = vector<vector<vector<double>>>(nperiods, rows);
  int numparams = 2 + (nareas * nareas * nperiods) - nareas;
  double init_x[numparams];

  double f;
  double low[numparams];
  // double up[numparams];

  for (int i = 0; i < numparams; i++) {
    init_x[i] = 0.01;
    low[i] = 1e-7;
    // up[i] = HUGE_VAL;
  }

  nlopt_opt opt = nlopt_create(NLOPT_LN_SBPLX, numparams);
  nlopt_set_lower_bounds(opt, low);
  nlopt_set_min_objective(
      opt, get_likelihood_with_optimized_dispersal_extinction, NULL);
  nlopt_set_xtol_rel(opt, 1e-7);
  nlopt_result res = nlopt_optimize(opt, init_x, &f);
  cout << "result ouput status: " << res << endl;
  vector<double> results;
  for (int i = 0; i < numparams; i++) {
    results.push_back(init_x[i]);
  }
  return results;
}
