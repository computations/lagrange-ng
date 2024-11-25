#ifndef FOSSIL_H
#define FOSSIL_H

#include <memory>

#include "Common.hpp"
#include "MRCA.hpp"
#include "Utils.hpp"

namespace lagrange {

enum class FossilType : uint8_t {
  NODE,
  BRANCH,
  FIXED,
  EXCLUDE,
  INCLUDE,
  UNKOWN
};

struct Fossil {
  MRCALabel mrca_name;
  Option<std::shared_ptr<MRCAEntry>> clade;
  double age;
  Range area;
  FossilType type;
};

}  // namespace lagrange
#endif  // !FOSSIL_H
