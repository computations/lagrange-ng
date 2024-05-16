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
  AncSplit(Dist, Dist, Dist, double);
  auto getWeight() const -> double;
  auto getLikelihood() const -> double;
  void setLikelihood(double li);

  void setLWR(double lwr);
  auto getLWR() const -> double;

  Dist anc_dist;
  Dist l_dist;
  Dist r_dist;
};

using lagrange_split_return_t =
    std::unordered_map<Dist, std::vector<AncSplit>>;

using lagrange_split_list_t = std::vector<lagrange_split_return_t>;

}  // namespace lagrange
#endif /* ANCSPLIT_H_ */
