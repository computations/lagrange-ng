#ifndef MRCA_H
#define MRCA_H

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace lagrange {

using Clade = std::vector<std::string>;

struct MRCAEntry {
  Clade clade;

  bool in(const std::string &key) const {
    return std::any_of(clade.begin(),
                       clade.end(),
                       [&key](const std::string &s) { return s == key; });
  }
};

using MRCALabel = std::string;
using MRCAMap = std::unordered_map<MRCALabel, std::shared_ptr<MRCAEntry>>;
}  // namespace lagrange

#endif  // !MRCA_H
