#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <filesystem>
#include <istream>
#include <string>
#include <unordered_map>

#include "Common.hpp"
#include "Utils.hpp"

namespace lagrange {
struct Alignment {
  std::unordered_map<std::string, Dist> data;
  size_t region_count;
  size_t taxa_count;
};

enum class AlignmentFileType { FASTA, PHYLIP };

Alignment read_fasta(std::istream& instream);
Alignment read_phylip(std::istream& instream);

Alignment read_alignment(std::istream& instream, AlignmentFileType type);
Alignment read_alignment(const std::filesystem::path& infile,
                         Option<AlignmentFileType> type);
}  // namespace lagrange

#endif
