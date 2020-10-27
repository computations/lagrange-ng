/*
 * BranchSegment.h
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef BRANCHSEGMENT_H_
#define BRANCHSEGMENT_H_

#include <vector>
#include <memory>
using namespace std;

#include "superdouble.h"
#include "vector_node_object.h"
#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif

class BranchSegment {
private:
  double _duration;
  int _period;
  vector<int> _fossil_area_indices;
  int _start_dist;

public:
  BranchSegment(double dur, int per);
  void clearStartDist();
  double getDuration();
  int getPeriod();
  void set_start_dist_int(int d);
  int get_start_dist_int();
  vector<int> getFossilAreas();
  void setFossilArea(int area);
#ifdef BIGTREE
  std::shared_ptr<VectorNodeObject<mpfr_class>> distconds;
  VectorNodeObject<mpfr_class> alphas;
  std::shared_ptr<VectorNodeObject<mpfr_class>>
      ancdistconds; // for ancestral state reconstructions
  VectorNodeObject<mpfr_class> seg_sp_alphas;
  VectorNodeObject<mpfr_class> seg_sp_stoch_map_revB_time;
  VectorNodeObject<mpfr_class> seg_sp_stoch_map_revB_number;
#else
  std::shared_ptr<vector<Superdouble>> distconds;
  vector<Superdouble> alphas; // alpha for the entire branch -- stored in the
                              // 0th segment for anc calc
  vector<Superdouble> seg_sp_alphas; // alpha for this specific segment, stored
                                     // for the stoch map
  vector<Superdouble>
      seg_sp_stoch_map_revB_time; // segment specific rev B, combining the tempA
                                  // and the ENLT
  vector<Superdouble>
      seg_sp_stoch_map_revB_number; // segment specific rev B, combining the
                                    // tempA and the ENLT
  std::shared_ptr<vector<Superdouble>>
      ancdistconds; // for ancestral state reconstructions
#endif
};

#endif /* BRANCHSEGMENT_H_ */
