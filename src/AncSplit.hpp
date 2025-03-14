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

#ifndef ANCSPLIT_H
#define ANCSPLIT_H

#include <unordered_map>
#include <vector>

#include "Common.hpp"

namespace lagrange {

class AncSplit {
 private:
  double _weight;
  double _likelihood;
  double _lwr;

 public:
  AncSplit(Range, Range, Range, double);

  [[nodiscard]] auto getWeight() const -> double;

  [[nodiscard]] auto getLikelihood() const -> double;

  void setLikelihood(double likelihood);

  void setLWR(double lwr);

  [[nodiscard]] auto getLWR() const -> double;

  Range anc_dist;
  Range l_dist;
  Range r_dist;
};

using SplitReturn = std::unordered_map<Range, std::vector<AncSplit>>;

using SplitReturnList = std::unordered_map<size_t, SplitReturn>;

}  // namespace lagrange
#endif /* ANCSPLIT_H_ */
