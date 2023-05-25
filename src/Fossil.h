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

class FossilList {
 public:
  void add_fossil(const std::string &mrca, const std::string &type,
                  const std::string &area, const std::string age) {
    _fossils.emplace_back(Fossil{.mrca = mrca, .area = area});
  }

 private:
  std::vector<Fossil> _fossils;
};
#endif  // !FOSSIL_H
