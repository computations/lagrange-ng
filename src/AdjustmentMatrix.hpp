#ifndef ADJUSTMENT_MATRIX_H
#define ADJUSTMENT_MATRIX_H

#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include "CSV.hpp"

namespace lagrange {

enum class AdjustmentMatrixType {
  symmetric,
  nonsymmetric,
  invalid,
};

struct AdjustmentArc {
  size_t from;
  size_t to;
  double dist;
};

class AdjustmentMatrix {
 public:
  AdjustmentMatrix(const std::filesystem::path& mat_filename,
                   const std::vector<std::string>& area_names);

  AdjustmentMatrix(std::istringstream&& instream,
                   const std::vector<std::string>& area_names);

  std::shared_ptr<double[]> to_matrix() const;

  size_t compute_size() const;
  size_t compute_index(size_t from, size_t to) const;

  AdjustmentMatrixType type() { return _type; }

 private:
  static AdjustmentMatrixType determine_matrix_symmetry(
      const std::vector<AdjustmentArc>& arcs, size_t region_count);

  void read_arcs(CSVReader&, const std::vector<std::string>&);

  size_t _region_count;
  AdjustmentMatrixType _type;
  std::vector<AdjustmentArc> _arcs;
};

}  // namespace lagrange
#endif
