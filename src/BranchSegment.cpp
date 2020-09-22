/*
 * BranchSegment.cpp
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 */

#include "BranchSegment.h"

using namespace std;

BranchSegment::BranchSegment(double dur, int per)
    : _duration(dur), _period(per),
      _fossil_area_indices(vector<int>()), _start_dist(-666),
      distconds(nullptr), ancdistconds(nullptr) {}

void BranchSegment::clearStartDist() {
  _start_dist = -666; // null is -666
}

double BranchSegment::getDuration() { return _duration; }

int BranchSegment::getPeriod() { return _period; }

void BranchSegment::set_start_dist_int(int d) { _start_dist = d; }

int BranchSegment::get_start_dist_int() { return _start_dist; }

vector<int> BranchSegment::getFossilAreas() { return _fossil_area_indices; }

void BranchSegment::setFossilArea(int area) {
  _fossil_area_indices.push_back(area);
}
