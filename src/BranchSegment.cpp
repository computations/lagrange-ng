/*
 * BranchSegment.cpp
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 */

#include "BranchSegment.h"
#include "RateModel.h"

#include <vector>
using namespace std;

BranchSegment::BranchSegment(double dur, int per)
    : _duration(dur), _period(per), _model(NULL),
      _fossil_area_indices(vector<int>()), _start_dist(-666), distconds(NULL),
      ancdistconds(NULL) {}

void BranchSegment::setModel(std::shared_ptr<RateModel> mod) { _model = mod; }

void BranchSegment::clearStartDist() {
  _start_dist = -666; // null is -666
}

double BranchSegment::getDuration() { return _duration; }

int BranchSegment::getPeriod() { return _period; }

void BranchSegment::set_start_dist_int(int d) { _start_dist = d; }

int BranchSegment::get_start_dist_int() { return _start_dist; }

std::shared_ptr<RateModel> BranchSegment::getModel() { return _model; }

vector<int> BranchSegment::getFossilAreas() { return _fossil_area_indices; }

void BranchSegment::setFossilArea(int area) {
  _fossil_area_indices.push_back(area);
}
