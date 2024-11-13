#include "Periods.hpp"

namespace lagrange {
void Periods::setMaxAreas(RangeSize m) { _max_areas = m; }

RangeSize Periods::maxAreas() const { return _max_areas; }

void Periods::setRegionCount(RangeSize r) { _regions = r; }

RangeSize Periods::regions() const { return _regions; }
}  // namespace lagrange
