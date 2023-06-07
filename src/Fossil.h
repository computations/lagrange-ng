#ifndef FOSSIL_H
#define FOSSIL_H

#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "Common.h"
#include "MRCA.h"
#include "Utils.h"

enum class fossil_type { node, branch };

struct Fossil {
  std::string mrca_name;
  lagrange_option_t<std::shared_ptr<MRCAEntry>> clade;
  double age;
  lagrange_dist_t area;
  fossil_type type;
};

#endif  // !FOSSIL_H
