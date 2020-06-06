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
#include "OptimizeBioGeoAllDispersal.h"
#include "RateModel.h"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

OptimizeBioGeoAllDispersal::OptimizeBioGeoAllDispersal(BioGeoTree *intree,
                                                       RateModel *inrm,
                                                       bool marg)
    : tree(intree), rm(inrm), maxiterations(10000), stoppingprecision(0.001),
      marginal(marg) {
  nareas = rm->get_num_areas();
  vector<double> cols(nareas, 1);
  vector<vector<double>> rows(nareas, cols);
  D_mask = vector<vector<vector<double>>>(rm->get_num_periods(), rows);
}

double
OptimizeBioGeoAllDispersal::GetLikelihoodWithOptimizedDispersalExtinction(
    const gsl_vector *variables) {
  double dispersal = gsl_vector_get(variables, 0);
  double extinction = gsl_vector_get(variables, 1);
  int count = 2;
  for (unsigned int i = 0; i < D_mask.size(); i++) {
    for (unsigned int j = 0; j < D_mask[i].size(); j++) {
      D_mask[i][j][j] = 1.0;
      for (unsigned int k = 0; k < D_mask[i][j].size(); k++) {
        if (k != j) {
          D_mask[i][j][k] = gsl_vector_get(variables, count);
          if (D_mask[i][j][k] <= 0) {
            return 100000000;
          }
          count += 1;
        }
      }
    }
  }
  double like;
  if (dispersal <= 0 || extinction <= 0)
    return 100000000;
  if (dispersal > 100 || extinction > 100)
    return 100000000;
  rm->setup_D_provided(dispersal, D_mask);
  rm->setup_E(extinction);
  rm->setup_Q();
  tree->update_default_model(rm);
  like = tree->eval_likelihood(marginal);
  if (like < 0 || like == std::numeric_limits<double>::infinity())
    like = 100000000;
  // cout << "dis: "<< dispersal << " ext: " << extinction << " like: "<< like
  // << endl;
  return like;
}

double
OptimizeBioGeoAllDispersal::GetLikelihoodWithOptimizedDispersalExtinction_gsl(
    const gsl_vector *variables, void *obj) {
  double temp;
  temp = ((OptimizeBioGeoAllDispersal *)obj)
             ->GetLikelihoodWithOptimizedDispersalExtinction(variables);
  return temp;
}

/*
 * USES THE SIMPLEX ALGORITHM
 *
 */
vector<double>
OptimizeBioGeoAllDispersal::optimize_global_dispersal_extinction() {
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  size_t np = 2 + (nareas * nareas * rm->get_num_periods()) - nareas;
  size_t iter = 0, i;
  int status;
  double size;
  /* Initial vertex size vector */
  ss = gsl_vector_alloc(np);
  /* Set all step sizes to .01 */ // Note that it was originally 1
  gsl_vector_set_all(ss, .1);
  /* Starting point */
  // cout<<"Now in OPtimizaRateWithGivenTipVariance in OptimizationFn"<<endl;
  x = gsl_vector_alloc(np);
  for (int i = 0; i < np; i++) {
    gsl_vector_set(x, i, 0.01);
  }
  OptimizeBioGeoAllDispersal *pt;
  pt = (this);
  double (*F)(const gsl_vector *, void *);
  F = &OptimizeBioGeoAllDispersal::
          GetLikelihoodWithOptimizedDispersalExtinction_gsl;
  /* Initialize method and iterate */
  gsl_multimin_function minex_func;
  minex_func.f = *F;
  minex_func.params = pt;
  minex_func.n = np;
  s = gsl_multimin_fminimizer_alloc(T, np);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
  do {
    // cout<<"Now on iteration "<<iter<<endl;
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status != 0) { // 0 Means it's a success
      printf("error: %s\n", gsl_strerror(status));
      break;
    }
    size = gsl_multimin_fminimizer_size(s);
    // status = gsl_multimin_test_size (size, 1e-2);
    status = gsl_multimin_test_size(
        size, stoppingprecision); // since we want more precision
    if (status == GSL_SUCCESS) {
      // printf ("converged to minimum at\n");
    }
    // printf ("%5d ", iter);
    for (i = 0; i < np; i++) {
      // printf ("%10.3e ", gsl_vector_get (s->x, i));
    }
    // printf ("f() = %7.3f size = %.3f\n", s->fval, size);
  } while (status == GSL_CONTINUE && iter < maxiterations);
  vector<double> results;
  // results.push_back(gsl_vector_get(s->x,0));
  // results.push_back(gsl_vector_get(s->x,1));
  for (int i = 0; i < np; i++) {
    results.push_back(gsl_vector_get(s->x, i));
  }
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  // cout << "dis: " << results[0] << " ext: " << results[1] << endl;
  return results;
}
