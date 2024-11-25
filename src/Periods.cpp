#include "Periods.hpp"

namespace lagrange {
void Periods::setMaxAreas(size_t m) { _max_areas = m; }

auto Periods::maxAreas() const -> size_t { return _max_areas; }

void Periods::setRegionCount(size_t r) { _regions = r; }

auto Periods::regions() const -> size_t { return _regions; }
}  // namespace lagrange
