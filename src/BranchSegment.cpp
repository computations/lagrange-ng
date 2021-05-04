/*
 * BranchSegment.cpp
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include "BranchSegment.h"

BranchSegment::BranchSegment(double dur, int per)
    : _duration(dur), _period(per), _fossil_area_indices{}, _start_dist(-666) {}

void BranchSegment::clearStartDist() {
  _start_dist = -666;  // null is -666
}

double BranchSegment::getDuration() const { return _duration; }

int BranchSegment::getPeriod() const { return _period; }

void BranchSegment::set_start_dist_int(int d) { _start_dist = d; }

int BranchSegment::get_start_dist_int() { return _start_dist; }

std::vector<int> BranchSegment::getFossilAreas() {
  return _fossil_area_indices;
}

void BranchSegment::setFossilArea(int area) {
  _fossil_area_indices.push_back(area);
}
