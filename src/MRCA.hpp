#ifndef MRCA_H
#define MRCA_H

#include <string>
#include <vector>

namespace lagrange {
struct MRCAEntry {
  std::vector<std::string> clade;
};
}  // namespace lagrange

#endif  // !MRCA_H
