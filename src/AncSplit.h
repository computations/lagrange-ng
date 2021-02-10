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

#include "Common.h"

class AncSplit {
 private:
  double _weight;
  double _likelihood;

 public:
  AncSplit(lagrange_dist_t, lagrange_dist_t, lagrange_dist_t, double);
  double getWeight() const;
  double getLikelihood() const;
  void setLikelihood(double li);
  lagrange_dist_t ancdistint;
  lagrange_dist_t ldescdistint;
  lagrange_dist_t rdescdistint;
};

#endif /* ANCSPLIT_H_ */
