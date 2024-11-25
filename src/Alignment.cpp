#include "Alignment.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <logger.hpp>
#include <sstream>
#include <string>

#include "Utils.hpp"

namespace lagrange {

/**
 * Clean the taxa name by removing whitespace from the front and back.
 */
template <typename T>
std::string clean_taxa_name(T start_itr, T end_itr) {
  while (std::isspace(*start_itr) != 0) { start_itr += 1; }
  while (std::isspace(*end_itr) != 0) { end_itr -= 1; }

  return {start_itr, end_itr};
}

std::tuple<TaxaName, Range, size_t> process_fasta_block(
    const std::string& block) {
  auto taxa_name_start = std::find(block.begin(), block.end(), '>') + 1;
  auto taxa_name_end = std::find(taxa_name_start, block.end(), ' ');
  TaxaName taxa_name = clean_taxa_name(taxa_name_start, taxa_name_end);

  std::string range_string;
  for (auto cur_itr = taxa_name_end + 1; cur_itr != block.end(); ++cur_itr) {
    if (std::isspace(*cur_itr) != 0) { continue; }
    range_string.push_back(*cur_itr);
  }

  Range range = convert_dist_binary_string_to_dist(range_string);

  return {taxa_name, range, range_string.size()};
}

std::string get_fasta_block(std::istream& instream) {
  std::string line;
  std::stringstream block_string;
  while (instream) {
    if (instream.peek() == '>' && !block_string.view().empty()) { break; }
    std::string line;
    std::getline(instream, line);
    block_string << line << " ";
  }
  return block_string.str();
}

Alignment read_fasta(std::istream& instream) {
  Alignment alignment;

  std::optional<size_t> region_count;

  std::string line;
  TaxaName taxa_name;

  bool good = true;
  while (instream) {
    auto block = get_fasta_block(instream);
    auto res = process_fasta_block(block);
    if (!region_count) { region_count = std::get<2>(res); }
    if (region_count.value() != std::get<2>(res)) {
      good = false;
      LOG(ERROR,
          "taxa {} has a different range size then the rest of the taxa",
          std::get<0>(res))  // NOLINT;
    }
    alignment.data[std::get<0>(res)] = std::get<1>(res);
  }
  if (!good) { throw AlignmentReadError{"Ranges in datafile vary in size"}; }

  alignment.region_count = region_count.value();
  alignment.taxa_count = alignment.data.size();
  return alignment;
}

Alignment read_phylip(std::istream& instream) {
  Alignment alignment;

  TaxaName taxa_name;
  std::string data_string;

  instream >> alignment.taxa_count;
  instream >> alignment.region_count;

  bool good = true;
  while (instream >> taxa_name) {
    instream >> data_string;
    if (data_string.size() != alignment.region_count) {
      LOG(ERROR,
          "Range for taxa '{}' has a different size than specified in the "
          "phylip header",
          taxa_name)  // NOLINT;
      good = false;
    }
    auto taxa_name_cleaned =
        clean_taxa_name(taxa_name.begin(), taxa_name.end());
    alignment.data[taxa_name_cleaned] =
        convert_dist_binary_string_to_dist(data_string);
  }
  if (!good) { throw AlignmentReadError{"Range size differs in phylip file"}; }

  return alignment;
}

Alignment read_alignment(std::istream& instream, AlignmentFileType type) {
  if (type == AlignmentFileType::FASTA) { return read_fasta(instream); }
  return read_phylip(instream);
}

/**
 * Reads the alignment from the specified file. Guesses the type based on the
 * file extension.
 */
Alignment read_alignment(const std::filesystem::path& filename,
                         AlignmentFileType type) {
  if (!std::filesystem::exists(filename)) {
    throw AlignmentReadError{"Failed to find the alignment file "
                             + filename.string()};
  }
  std::ifstream alignment_file(filename);
  return read_alignment(alignment_file, type);
}
}  // namespace lagrange
