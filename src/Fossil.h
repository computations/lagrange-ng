#ifndef FOSSIL_H
#define FOSSIL_H

#include <memory>
#include <string>
#include <vector>

#include "Common.h"
#include "MRCA.h"
#include "Utils.h"

enum class fossil_type { n, b };

struct Fossil {
  std::string mrca_name;
  lagrange_option_t<std::shared_ptr<MRCAEntry>> clade;
  lagrange_dist_t area;
  fossil_type type;
  double age;
};

#endif  // !FOSSIL_H
