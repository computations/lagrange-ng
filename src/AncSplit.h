/*
 * AncSplit.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 */

/*
  The AncSplit (ancestral split) class is used to store the likelihoods and
  distributions associated with an ancstral split.

  This should only be used for ancestral state calculation as there is no
  need to store the likelihoods for each state (other than distconds) when
  calculating the likeklihood.
 */

#ifndef ANCSPLIT_H_
#define ANCSPLIT_H_

#include <vector>
#include <memory>

#include "RateModel.h"
#include "superdouble.h"

class AncSplit {
private:
  std::shared_ptr<RateModel> _model;
  double _weight;
  Superdouble _likelihood;

public:
  AncSplit(std::shared_ptr<RateModel> mod, int, int, int, Superdouble);
  std::shared_ptr<RateModel> getModel();
  double getWeight();
  Superdouble getLikelihood();
  void setLikelihood(Superdouble li);
  int ancdistint;
  int ldescdistint;
  int rdescdistint;
};

#endif /* ANCSPLIT_H_ */
