#ifndef FOSSIL_H
#define FOSSIL_H

#include <memory>
#include <string>

#include "Common.hpp"
#include "MRCA.hpp"
#include "Utils.hpp"

enum class fossil_type { NODE, BRANCH, FIXED };

struct Fossil {
  std::string mrca_name;
  LagrangeOption<std::shared_ptr<MRCAEntry>> clade;
  double age;
  lagrange_dist_t area;
  fossil_type type;
};

#endif  // !FOSSIL_H
