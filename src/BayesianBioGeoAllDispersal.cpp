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
    : tree(intree), rm(inrm), gens(gen), marginal(marg) {
  int nareas = rm->get_num_areas();
  nareas = rm->get_num_areas();
  vector<double> cols(nareas, 1);
  vector<vector<double>> rows(nareas, cols);
  D_mask = vector<vector<vector<double>>>(rm->get_num_periods(), rows);
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
}

void BayesianBioGeoAllDispersal::run_global_dispersal_extinction() {
  double prevlike = 0;
  double prevprior = 1;
  double prevpost = 0;
  double curlike = 0;
  double curpost = 0;
  double curprior = 0;

  int nareas = rm->get_num_areas();
  int nparams = 2 + (nareas * nareas * rm->get_num_periods()) -
                (rm->get_num_periods() * nareas);

  vector<double> sliding(nparams);
  vector<double> trials(nparams);
  vector<double> success(nparams);
  for (unsigned int i = 0; i < sliding.size(); i++) {
    sliding[i] = 0.001;
  }
  // sliding[0] = 0.005;sliding[1] = 0.005;
  for (unsigned int i = 0; i < trials.size(); i++) {
    trials[i] = 0;
  }
  for (unsigned int i = 0; i < success.size(); i++) {
    success[i] = 0;
  }

  int rot = 1;
  double hastings = 1;

  ofstream outfile("test.txt");

  params = vector<double>(nparams);
  prevparams = vector<double>(nparams);
  for (unsigned int i = 0; i < params.size(); i++) {
    params[i] = 0.01;
  }
  params[0] = 1.;
  params[1] = 5.28047e-07;
  rm->setup_D(0.1);
  rm->setup_E(0.1);
  rm->setup_Q();
  tree->update_default_model(rm);
  prevlike = -tree->eval_likelihood(marginal);
  int iter = 0;
  while (iter < gens) {
    double dispersal = params[0];
    double extinction = params[1];
    int count = 2;
    for (unsigned int i = 0; i < D_mask.size(); i++) {
      for (unsigned int j = 0; j < D_mask[i].size(); j++) {
        D_mask[i][j][j] = 0.0;
        for (unsigned int k = 0; k < D_mask[i][j].size(); k++) {
          if (k != j) {
            D_mask[i][j][k] = params[count];
            /*if(D_mask[i][j][k] <= 0){
                    return 100000000;
            }*/
            count += 1;
          }
        }
      }
    }
    rm->setup_D_provided(dispersal, D_mask);
    rm->setup_E(extinction);
    rm->setup_Q();
    tree->update_default_model(rm);
    curlike = -tree->eval_likelihood(marginal);

    /*
     * calcprior
     */
    curprior = 1;
    // for(unsigned int i=0;i<params.size();i++){curprior *=
    // calculate_pdf(params[i]);}

    /*
     * check to keep
     */
    double testr = gsl_ran_flat(r, 0, 1);
    double test =
        MIN(1, hastings * (curprior / prevprior) * (exp(curlike - prevlike)));

    if (iter > 1000)
      trials[rot] += 1;
    //		cout << "- "<< prevlike << " " <<  curlike << " " << (curprior) << "
    //"<< test << " " << testr << endl;
    if (testr < test) {
      prevprior = curprior;
      prevlike = curlike;
      prevpost = curpost;
      prevparams = params;
      if (iter > 1000)
        success[rot] += 1;
    }
    /*
     * pick next params
     */
    // params[rot] = calculate_sliding(prevparams[rot],sliding[rot]);
    params[rot] =
        calculate_sliding_log(prevparams[rot], sliding[rot], &hastings);
    if (iter % 10 == 0) {
      rot += 1;
    }
    if (rot == params.size()) {
      rot = 1;
    }
    if (iter % 100 == 0 && iter > 1) {
      cout << iter << " " << prevlike;
      for (unsigned int i = 0; i < 2; i++) {
        cout << " " << prevparams[i];
      }
      cout << endl;
      for (unsigned int i = 0; i < D_mask.size(); i++) {
        for (unsigned int j = 0; j < D_mask[i].size(); j++) {
          for (unsigned int k = 0; k < D_mask[i][j].size(); k++) {
            //							cout << D_mask[i][j][k] << "
            //";
          }
          //					cout << endl;
        }
        //				cout << endl;
      }
      //			for(unsigned int i=0;i<params.size();i++){cout
      //<< " " << success[i]/trials[i];} 			cout << endl;
      outfile << iter << "\t" << prevlike;
      for (unsigned int i = 0; i < params.size(); i++) {
        outfile << "\t" << prevparams[i];
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
  // return abs(gsl_ran_flat(r,value - (w/2),value + (w/2)));
  return abs(value + gsl_ran_gaussian(r, sliding));
  // gsl_gaus
  // cauchy
}

double BayesianBioGeoAllDispersal::calculate_sliding_log(double value,
                                                         double sliding,
                                                         double *hastings) {
  double newv = log(value) + gsl_ran_gaussian(r, sliding);
  newv = abs(exp(newv));
  *hastings = 1 * newv / value;
  return newv;
}
