/*
 *  BayesianBioGeo.cpp
 *  lagrange_cpp
 *
 *  Created by Stephen Smith on 1/13/10.
 *  Copyright 2010 Yale University. All rights reserved.
 *
 */

#include "BayesianBioGeo.h"

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
#include <vector>
using namespace std;

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "BioGeoTree.h"
#include "RateModel.h"

namespace {
inline double MIN(const double &a, const double &b) {
  return b < a ? (b) : double(a);
}
} // namespace

BayesianBioGeo::BayesianBioGeo(BioGeoTree *intree, RateModel *inrm, bool marg,
                               int gen)
    : _tree(intree), _rate_model(inrm), _generations(gen), _marginal(marg) {
  gsl_rng_env_setup();
  _rng_type = gsl_rng_default;
  _rng = gsl_rng_alloc(_rng_type);
}

void BayesianBioGeo::run_global_dispersal_extinction() {
  double prevlike = 0;
  double prevprior = 0;
  double prevpost = 0;
  double curlike = 0;
  double curpost = 0;
  double curprior = 0;

  vector<double> sliding(2);
  vector<double> trials(2);
  vector<double> success(2);
  sliding[0] = 0.5;
  sliding[1] = 0.5;
  for (unsigned int i = 0; i < trials.size(); i++) {
    trials[i] = 0;
  }
  for (unsigned int i = 0; i < success.size(); i++) {
    success[i] = 0;
  }

  int rot = 0;
  double hastings = 1;

  ofstream outfile("test.txt");

  _params = vector<double>(2);
  _prev_params = vector<double>(2);
  for (unsigned int i = 0; i < _params.size(); i++) {
    _params[i] = 0.1;
  }
  _rate_model->setup_D(0.1);
  _rate_model->setup_E(0.1);
  _rate_model->setup_Q();
  _tree->update_default_model(_rate_model);
  prevlike = -_tree->eval_likelihood(_marginal);
  prevprior = 1;
  for (unsigned int i = 0; i < _params.size(); i++) {
    prevprior *= calculate_pdf(_params[i]);
  }
  cout << "intial like: " << prevlike << endl;
  int iter = 0;
  while (iter < _generations) {
    double dispersal = _params[0];
    double extinction = _params[1];
    _rate_model->setup_D(dispersal);
    _rate_model->setup_E(extinction);
    _rate_model->setup_Q();
    _tree->update_default_model(_rate_model);
    curlike = -_tree->eval_likelihood(_marginal);
    /*
     * calcprior
     */
    curprior = 1;
    for (unsigned int i = 0; i < _params.size(); i++) {
      curprior *= calculate_pdf(_params[i]);
    }

    /*
     * check to keep
     */
    double testr = gsl_ran_flat(_rng, 0, 1);
    double test =
        MIN(1, hastings * (curprior / prevprior) * (exp(curlike - prevlike)));

    if (iter > 100)
      trials[rot] += 1;
    if (testr < test) {
      prevprior = curprior;
      prevlike = curlike;
      prevpost = curpost;
      _prev_params[rot] = _params[rot];
      if (iter > 100)
        success[rot] += 1;
    }
    /*
     * pick next params
     */
    _params[rot] =
        calculate_sliding_log(_prev_params[rot], sliding[rot], &hastings);
    if (iter % 5 == 0) {
      rot += 1;
    }
    if (rot == 2) {
      rot = 0;
    }
    if (iter % 100 == 0 && iter > 100) {
      cout << iter << " " << prevlike;
      for (unsigned int i = 0; i < _params.size(); i++) {
        cout << " " << _prev_params[i];
      }
      for (unsigned int i = 0; i < _params.size(); i++) {
        cout << " " << success[i] / trials[i];
      }
      cout << endl;
      outfile << iter << "\t" << prevlike;
      for (unsigned int i = 0; i < _params.size(); i++) {
        outfile << "\t" << _prev_params[i];
      }
      outfile << endl;
    }
    iter++;
  }
  outfile.close();
}

double BayesianBioGeo::calculate_pdf(double value) {
  return gsl_ran_flat_pdf(value, 0, 100.);
}

double BayesianBioGeo::calculate_sliding(double value, double sliding) {
  // return abs(gsl_ran_flat(r,value - (sliding/2),value + (sliding/2)));
  return abs(value + gsl_ran_gaussian(_rng, sliding));
  // cauchy
}

double BayesianBioGeo::calculate_sliding_log(double value, double sliding,
                                             double *hastings) {
  double newv = log(value) + gsl_ran_gaussian(_rng, sliding);
  newv = abs(exp(newv));
  *hastings = 1 * newv / value;
  return newv;
}
