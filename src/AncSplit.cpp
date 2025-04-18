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

#include "AncSplit.hpp"

#include "Common.hpp"

namespace lagrange {
AncSplit::AncSplit(Range dist, Range ldesc, Range rdesc, double weight) :
    _weight(weight),
    _likelihood(0.0),
    anc_dist(dist),
    l_dist(ldesc),
    r_dist(rdesc) {}

auto AncSplit::getWeight() const -> double { return _weight; }

void AncSplit::setLikelihood(double likelihood) { _likelihood = likelihood; }

auto AncSplit::getLikelihood() const -> double { return _likelihood; }

void AncSplit::setLWR(double lwr) { _lwr = lwr; }

auto AncSplit::getLWR() const -> double { return _lwr; }
}  // namespace lagrange
