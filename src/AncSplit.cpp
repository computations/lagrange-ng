/*
 * AncSplit.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
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

#include "AncSplit.h"
#include "Common.h"

AncSplit::AncSplit(lagrange_dist_t dist, lagrange_dist_t ldesc,
                   lagrange_dist_t rdesc, double we)
    : _weight(we),
      _likelihood(0.0),
      anc_dist(dist),
      l_dist(ldesc),
      r_dist(rdesc) {}

auto AncSplit::getWeight() const -> double { return _weight; }

void AncSplit::setLikelihood(double li) { _likelihood = li; }

auto AncSplit::getLikelihood() const -> double { return _likelihood; }

void AncSplit::setLWR(double lwr) { _lwr = lwr; }

auto AncSplit::getLWR() const -> double { return _lwr; }
