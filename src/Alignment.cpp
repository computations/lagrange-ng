#include "Alignment.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "Utils.hpp"

namespace lagrange {

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
  std::string taxa_name;
  std::stringstream data_string;

  while (std::getline(instream, line)) {
    if (line[0] == '>') {
      if (!taxa_name.empty()) {
        auto tmp = data_string.str();
        region_count = tmp.size();
        alignment.data[taxa_name] =
            lagrange_convert_dist_binary_string_to_dist(tmp);

        taxa_name.clear();
        data_string = std::stringstream();
      }

      taxa_name = clean_taxa_name(std::string(line.begin() + 1, line.end()));
    } else {
      data_string << line;
    }
    line_number++;
  }

  alignment.data[taxa_name] =
      lagrange_convert_dist_binary_string_to_dist(data_string.str());

  alignment.region_count = region_count;
  alignment.taxa_count = alignment.data.size();
  return alignment;
}

Alignment read_phylip(std::istream& instream) {
  Alignment alignment;

  std::string taxa_name;
  std::string data_string;

  instream >> alignment.taxa_count;
  instream >> alignment.region_count;

  while (instream >> taxa_name) {
    instream >> data_string;
    alignment.data[clean_taxa_name(taxa_name)] =
        lagrange_convert_dist_binary_string_to_dist(data_string);
  }

  return alignment;
}

Alignment read_alignment(std::istream& infile, AlignmentFileType type) {
  if (type == AlignmentFileType::FASTA) {
    return read_fasta(infile);
  } else {
    return read_phylip(infile);
  }
}

Alignment read_alignment(const std::filesystem::path& filename,
                         Option<AlignmentFileType> type) {
  if (!std::filesystem::exists(filename)) {
    throw std::runtime_error{"Failed to find the alignment file "
                             + filename.string()};
  }
  std::ifstream alignment_file(filename);
  if (type.hasValue()) { return read_alignment(alignment_file, type.get()); }

  auto extension = get_file_extension(filename);
  if (extension == "fasta" || extension == "fas") {
    return read_alignment(alignment_file, AlignmentFileType::FASTA);
  } else if (extension == "phylip" || extension == "phy") {
    return read_alignment(alignment_file, AlignmentFileType::PHYLIP);
  }
  throw std::runtime_error{"Failed to recognize alignment file type"};
}
}  // namespace lagrange
