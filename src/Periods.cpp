#include "Periods.hpp"

namespace lagrange {
void Periods::setMaxAreas(size_t m) { _max_areas = m; }

size_t Periods::maxAreas() const { return _max_areas; }

void Periods::setRegionCount(size_t r) { _regions = r; }

size_t Periods::regions() const { return _regions; }
}  // namespace lagrange
