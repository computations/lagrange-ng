#include "Common.hpp"

#include "Utils.hpp"

namespace lagrange {
size_t NodeReservation::rangeCountForComputation() const {
  return lagrange_popcount(_child_region_mask);
}

}  // namespace lagrange
