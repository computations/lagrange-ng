/*
 * AncSplit.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */

/*
  The AncSplit (ancestral split) class is used to store the likelihoods and
  distributions associated with an ancstral split.

  This should only be used for ancestral state calculation as there is no
  need to store the likelihoods for each state (other than distconds) when
  calculating the likeklihood.
 */

#include "AncSplit.h"
#include "RateModel.h"
#include "superdouble.h"
#include <vector>
using namespace std;

AncSplit::AncSplit(RateModel *mod, int dist, int ldesc, int rdesc,
                   Superdouble we)
    : _model(mod), _weight(we), _likelihood(0), ancdistint(dist),
      ldescdistint(ldesc), rdescdistint(rdesc) {}

RateModel *AncSplit::getModel() { return _model; }

double AncSplit::getWeight() { return _weight; }

Superdouble AncSplit::getLikelihood() { return _likelihood; }

void AncSplit::setLikelihood(Superdouble li) { _likelihood = li; }
