/*
 *  BayesianBioGeo.cpp
 *  lagrange_cpp
 *
 *  Created by Stephen Smith on 1/13/10.
 *  Copyright 2010 Yale University. All rights reserved.
 *
 */

#include "BayesianBioGeoAllDispersal.h"

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

BayesianBioGeoAllDispersal::BayesianBioGeoAllDispersal(BioGeoTree *intree,
                                                       RateModel *inrm,
                                                       bool marg, int gen)
    : _tree(intree), _rate_model(inrm), _generations(gen), _marginal(marg) {
  int nareas = _rate_model->get_num_areas();
  nareas = _rate_model->get_num_areas();
  vector<double> cols(nareas, 1);
  vector<vector<double>> rows(nareas, cols);
  _dispersal_mask = vector<vector<vector<double>>>(_rate_model->get_num_periods(), rows);
  gsl_rng_env_setup();
  _rng_type = gsl_rng_default;
  _rng = gsl_rng_alloc(_rng_type);
}

void BayesianBioGeoAllDispersal::run_global_dispersal_extinction() {
  double prevlike = 0;
  double prevprior = 1;
  double curlike = 0;
  double curprior = 0;

  int nareas = _rate_model->get_num_areas();
  int nparams = 2 + (nareas * nareas * _rate_model->get_num_periods()) -
                (_rate_model->get_num_periods() * nareas);

  vector<double> sliding(nparams);
  vector<double> trials(nparams);
  vector<double> success(nparams);
  for (unsigned int i = 0; i < sliding.size(); i++) {
    sliding[i] = 0.001;
  }
  for (unsigned int i = 0; i < trials.size(); i++) {
    trials[i] = 0;
  }
  for (unsigned int i = 0; i < success.size(); i++) {
    success[i] = 0;
  }

  size_t rot = 1;
  double hastings = 1;

  ofstream outfile("test.txt");

  _params = vector<double>(nparams);
  _prev_params = vector<double>(nparams);
  for (unsigned int i = 0; i < _params.size(); i++) {
    _params[i] = 0.01;
  }
  _params[0] = 1.;
  _params[1] = 5.28047e-07;
  _rate_model->setup_D(0.1);
  _rate_model->setup_E(0.1);
  _rate_model->setup_Q();
  _tree->update_default_model(_rate_model);
  prevlike = -_tree->eval_likelihood(_marginal);
  int iter = 0;
  while (iter < _generations) {
    double dispersal = _params[0];
    double extinction = _params[1];
    int count = 2;
    for (unsigned int i = 0; i < _dispersal_mask.size(); i++) {
      for (unsigned int j = 0; j < _dispersal_mask[i].size(); j++) {
        _dispersal_mask[i][j][j] = 0.0;
        for (unsigned int k = 0; k < _dispersal_mask[i][j].size(); k++) {
          if (k != j) {
            _dispersal_mask[i][j][k] = _params[count];
            count += 1;
          }
        }
      }
    }
    _rate_model->setup_D_provided(dispersal, _dispersal_mask);
    _rate_model->setup_E(extinction);
    _rate_model->setup_Q();
    _tree->update_default_model(_rate_model);
    curlike = -_tree->eval_likelihood(_marginal);

    /*
     * calcprior
     */
    curprior = 1;

    /*
     * check to keep
     */
    double testr = gsl_ran_flat(_rng, 0, 1);
    double test =
        MIN(1, hastings * (curprior / prevprior) * (exp(curlike - prevlike)));

    if (iter > 1000)
      trials[rot] += 1;
    if (testr < test) {
      prevprior = curprior;
      prevlike = curlike;
      _prev_params = _params;
      if (iter > 1000)
        success[rot] += 1;
    }
    /*
     * pick next params
     */
    _params[rot] =
        calculate_sliding_log(_prev_params[rot], sliding[rot], &hastings);
    if (iter % 10 == 0) {
      rot += 1;
    }
    if (rot == _params.size()) {
      rot = 1;
    }
    if (iter % 100 == 0 && iter > 1) {
      cout << iter << " " << prevlike;
      for (unsigned int i = 0; i < 2; i++) {
        cout << " " << _prev_params[i];
      }
      cout << endl;
      for (unsigned int i = 0; i < _dispersal_mask.size(); i++) {
        for (unsigned int j = 0; j < _dispersal_mask[i].size(); j++) {
          for (unsigned int k = 0; k < _dispersal_mask[i][j].size(); k++) {
          }
        }
      }
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

double BayesianBioGeoAllDispersal::calculate_pdf(double value) {
  return gsl_ran_flat_pdf(value, 0, 100.);
}

double BayesianBioGeoAllDispersal::calculate_sliding(double value,
                                                     double sliding) {
  return abs(value + gsl_ran_gaussian(_rng, sliding));
}

double BayesianBioGeoAllDispersal::calculate_sliding_log(double value,
                                                         double sliding,
                                                         double *hastings) {
  double newv = log(value) + gsl_ran_gaussian(_rng, sliding);
  newv = abs(exp(newv));
  *hastings = 1 * newv / value;
  return newv;
}
