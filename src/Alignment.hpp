#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <filesystem>
#include <istream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "Common.hpp"
#include "Utils.hpp"

namespace lagrange {

using TaxaName = std::string;

/**
 * Simple data struct contianing the "alignment", which is to say a map from
 * taxa name to extant range.
 */
struct Alignment {
  void apply_max_areas(size_t max_areas) {
    auto dist_map = invert_dist_map(region_count, max_areas);
    for (auto& kv : data) { kv.second = dist_map[kv.second]; }
  }

  std::unordered_map<TaxaName, Range> data;
  size_t region_count;
  size_t taxa_count;
};

enum class AlignmentFileType : uint8_t { FASTA, PHYLIP };

auto read_fasta(std::istream& instream) -> Alignment;
auto read_phylip(std::istream& instream) -> Alignment;

/**
 * Parsing functions alignments. There are two versions: the version with a
 * known filetype, and the version with an unknown filetype, which is
 * represented by an Option.
 */
auto read_alignment(std::istream& instream,
                    AlignmentFileType type) -> Alignment;
auto read_alignment(const std::filesystem::path& filename,
                    AlignmentFileType type) -> Alignment;

/* Error types */
class AlignmentReadError : std::runtime_error {
 public:
  AlignmentReadError(const std::string& msg) : std::runtime_error{msg} {}
};

class AlignmentFiletypeError : std::runtime_error {
 public:
  AlignmentFiletypeError(const std::string& msg) : std::runtime_error{msg} {}
};
}  // namespace lagrange

#endif
