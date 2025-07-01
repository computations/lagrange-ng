#include "Periods.hpp"

namespace lagrange {
void PeriodTimes::setMaxAreas(size_t m) { _max_areas = m; }

auto PeriodTimes::maxAreas() const -> size_t { return _max_areas; }

void PeriodTimes::setRegionCount(size_t r) { _regions = r; }

auto PeriodTimes::regions() const -> size_t { return _regions; }
}  // namespace lagrange
