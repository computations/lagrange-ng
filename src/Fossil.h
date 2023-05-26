#ifndef FOSSIL_H
#define FOSSIL_H

#include <string>
#include <vector>

#include "Common.h"
#include "Utils.h"

enum class fossil_type { n, b };

struct Fossil {
  std::string mrca;
  lagrange_dist_t area;
  fossil_type type;
  double age;
};

#endif  // !FOSSIL_H
