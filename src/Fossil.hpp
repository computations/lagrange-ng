#ifndef FOSSIL_H
#define FOSSIL_H

#include <memory>
#include <string>

#include "Common.hpp"
#include "MRCA.hpp"
#include "Utils.hpp"

namespace lagrange {

enum class fossil_type { NODE, BRANCH, FIXED };

struct Fossil {
  std::string mrca_name;
  Option<std::shared_ptr<MRCAEntry>> clade;
  double age;
  Dist area;
  fossil_type type;
};

}  // namespace lagrange
#endif  // !FOSSIL_H
