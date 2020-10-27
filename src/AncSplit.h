/*
 * AncSplit.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
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

#include <memory>
#include <vector>

#include "superdouble.h"

class AncSplit {
 private:
  double _weight;
  Superdouble _likelihood;

 public:
  AncSplit(int, int, int, double);
  double getWeight() const;
  Superdouble getLikelihood() const;
  void setLikelihood(Superdouble li);
  int ancdistint;
  int ldescdistint;
  int rdescdistint;
};

#endif /* ANCSPLIT_H_ */
