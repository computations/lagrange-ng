#include "Alignment.hpp"

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
std::string clean_taxa_name(const std::string& str) {
  auto start_itr = str.begin();
  auto end_itr = str.end() - 1;
  while (std::isspace(*start_itr)) { start_itr += 1; }
  while (std::isspace(*end_itr)) { end_itr -= 1; }

  return {start_itr, end_itr + 1};
}

Alignment read_fasta(std::istream& instream) {
  Alignment alignment;

  size_t line_number = 1;

  size_t region_count = 0;

  std::string line;
  TaxaName taxa_name;
  std::stringstream data_string;

  bool good = true;
  while (std::getline(instream, line)) {
    if (line[0] == '>') {
      if (!taxa_name.empty()) {
        auto tmp = data_string.str();
        if (region_count == 0) {
          region_count = tmp.size();
        } else if (region_count != tmp.size()) {
          LOG(ERROR,
              "The range size for taxa '%s' differs in size",
              taxa_name.c_str());
          good = false;
        }
        alignment.data[taxa_name] = convert_dist_binary_string_to_dist(tmp);

        taxa_name.clear();
        data_string = std::stringstream();
      }

      taxa_name = clean_taxa_name(std::string(line.begin() + 1, line.end()));
    } else {
      data_string << line;
    }
    line_number++;
  }
  if (!good) { throw AlignmentReadError{"Ranges in datafile vary in size"}; }

  alignment.data[taxa_name] =
      convert_dist_binary_string_to_dist(data_string.str());

  alignment.region_count = region_count;
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
          "Range for taxa '%s' has a different size than specified in the "
          "phylip header",
          taxa_name.c_str());
      good = false;
    }
    alignment.data[clean_taxa_name(taxa_name)] =
        convert_dist_binary_string_to_dist(data_string);
  }
  if (!good) { throw AlignmentReadError{"Range size differs in phylip file"}; }

  return alignment;
}

Alignment read_alignment(std::istream& infile, AlignmentFileType type) {
  if (type == AlignmentFileType::FASTA) {
    return read_fasta(infile);
  } else {
    return read_phylip(infile);
  }
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
