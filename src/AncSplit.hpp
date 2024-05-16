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

class AncSplit {
 private:
  double _weight;
  double _likelihood;
  double _lwr;

 public:
  AncSplit(lagrange_dist_t, lagrange_dist_t, lagrange_dist_t, double);
  auto getWeight() const -> double;
  auto getLikelihood() const -> double;
  void setLikelihood(double li);

  void setLWR(double lwr);
  auto getLWR() const -> double;

  lagrange_dist_t anc_dist;
  lagrange_dist_t l_dist;
  lagrange_dist_t r_dist;
};

using lagrange_split_return_t =
    std::unordered_map<lagrange_dist_t, std::vector<AncSplit>>;

using lagrange_split_list_t = std::vector<lagrange_split_return_t>;

#endif /* ANCSPLIT_H_ */
