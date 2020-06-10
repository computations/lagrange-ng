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
#include "OptimizeBioGeo.h"
#include "RateModel.h"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

OptimizeBioGeo::OptimizeBioGeo(BioGeoTree *intree, RateModel *inrm, bool marg)
    : _tree(intree), _rate_model(inrm), _marginal(marg) {}

double OptimizeBioGeo::GetLikelihoodWithOptimizedDispersalExtinction(
    const gsl_vector *variables) {
  double dispersal = gsl_vector_get(variables, 0);
  double extinction = gsl_vector_get(variables, 1);
  double like;
  if (dispersal <= 0 || extinction <= 0)
    return 100000000;
  if (dispersal > 100 || extinction > 100)
    return 100000000;
  _rate_model->setup_D(dispersal);
  _rate_model->setup_E(extinction);
  _rate_model->setup_Q();
  _tree->update_default_model(_rate_model);
  like = _tree->eval_likelihood(_marginal);
  if (like < 0 || like == std::numeric_limits<double>::infinity())
    like = 100000000;
  return like;
}

double OptimizeBioGeo::GetLikelihoodWithOptimizedDispersalExtinction_gsl(
    const gsl_vector *variables, void *obj) {
  double temp;
  temp = ((OptimizeBioGeo *)obj)
             ->GetLikelihoodWithOptimizedDispersalExtinction(variables);
  return temp;
}

/*
 * USES THE SIMPLEX ALGORITHM
 *
 */
vector<double> OptimizeBioGeo::optimize_global_dispersal_extinction() {
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  size_t np = 2;
  size_t iter = 0, i;
  int status;
  double size;
  /* Initial vertex size vector */
  ss = gsl_vector_alloc(np);
  /* Set all step sizes to .01 */ // Note that it was originally 1
  gsl_vector_set_all(ss, .01);
  /* Starting point */
  x = gsl_vector_alloc(np);
  gsl_vector_set(x, 0, 0.01);
  gsl_vector_set(x, 1, 0.01);
  OptimizeBioGeo *pt;
  pt = (this);
  double (*F)(const gsl_vector *, void *);
  F = &OptimizeBioGeo::GetLikelihoodWithOptimizedDispersalExtinction_gsl;
  /* Initialize method and iterate */
  gsl_multimin_function minex_func;
  minex_func.f = *F;
  minex_func.params = pt;
  minex_func.n = np;
  s = gsl_multimin_fminimizer_alloc(T, np);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status != 0) { // 0 Means it's a success
      printf("error: %s\n", gsl_strerror(status));
      break;
    }
    size = gsl_multimin_fminimizer_size(s);
    status =
        gsl_multimin_test_size(size, _abs_tol); // since we want more precision
    if (status == GSL_SUCCESS) {
    }
  } while (status == GSL_CONTINUE && iter < _max_iterations);
  vector<double> results;
  results.push_back(gsl_vector_get(s->x, 0));
  results.push_back(gsl_vector_get(s->x, 1));
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  return results;
}
