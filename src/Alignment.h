#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <istream>
#include <string>
#include <unordered_map>

#include "Common.h"
#include "Utils.h"

struct Alignment {
  std::unordered_map<std::string, lagrange_dist_t> data;
  size_t region_count;
  size_t taxa_count;
};

enum class AlignmentFileType { fasta, phylip };

enum class AlignmentFileType { fasta, phylip };

Alignment read_fasta(std::istream& instream);
Alignment read_phylip(std::istream& instream);

Alignment read_alignment(std::istream& instream, AlignmentFileType type);
Alignment read_alignment(const std::string& infile,
                         lagrange_option_t<AlignmentFileType> type);

#endif
