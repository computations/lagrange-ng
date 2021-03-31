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

#include <memory>
#include <vector>
using namespace std;

class BranchSegment {
 private:
  double _duration;
  int _period;
  vector<int> _fossil_area_indices;
  int _start_dist;

 public:
  BranchSegment(double dur, int per);
  void clearStartDist();
  double getDuration() const;
  int getPeriod() const;
  void set_start_dist_int(int d);
  int get_start_dist_int();
  vector<int> getFossilAreas();
  void setFossilArea(int area);
};

#endif /* BRANCHSEGMENT_H_ */
