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
  std::unordered_map<TaxaName, Range> data;
  size_t region_count;
  size_t taxa_count;

  Range rangeUnion() const;
  size_t usedRanges() const;
  bool allRegionsValid() const;
};

enum class AlignmentFileType { FASTA, PHYLIP };

Alignment read_fasta(std::istream& instream);
Alignment read_phylip(std::istream& instream);

/**
 * Parsing functions alignments. There are two versions: the version with a
 * known filetype, and the version with an unknown filetype, which is
 * represented by an Option.
 */
Alignment read_alignment(std::istream& instream, AlignmentFileType type);
Alignment read_alignment(const std::filesystem::path& infile,
                         AlignmentFileType type);

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
