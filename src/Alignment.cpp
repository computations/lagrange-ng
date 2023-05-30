#include "Alignment.h"

#include <fstream>
#include <sstream>
#include <string>

#include "Utils.h"

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
  std::string line;
  std::string taxa_name;
  std::stringstream data_string;

  while (std::getline(instream, line)) {
    if (line[0] == '>') {
      if (!taxa_name.empty()) {
        alignment[taxa_name] =
            lagrange_convert_dist_binary_string_to_dist(data_string.str());

        taxa_name.clear();
        data_string.clear();
      }

      taxa_name = clean_taxa_name(std::string(line.begin() + 1, line.end()));
    } else {
      data_string << line;
    }
    line_number++;
  }

  alignment[taxa_name] =
      lagrange_convert_dist_binary_string_to_dist(data_string.str());

  return alignment;
}

Alignment read_phylip(std::istream& instream) {
  Alignment alignment;

  std::string taxa_name;
  std::string data_string;

  std::string header_line;

  // The header is worthless, throw it away
  std::getline(instream, header_line);

  while (instream >> taxa_name) {
    instream >> data_string;
    alignment[clean_taxa_name(taxa_name)] =
        lagrange_convert_dist_binary_string_to_dist(data_string);
  }

  return alignment;
}

Alignment read_alignment(std::istream& infile, AlignmentFileType type) {
  if (type == AlignmentFileType::fasta) {
    return read_fasta(infile);
  } else {
    return read_phylip(infile);
  }
}

Alignment read_alignment(const std::string& infile,
                         lagrange_option_t<AlignmentFileType> type) {
  std::ifstream alignment_file(infile);
  if (type.has_value()) { return read_alignment(alignment_file, type.get()); }

  auto extension = get_file_extension(infile);
  if (extension == "fasta") {
    return read_alignment(alignment_file, AlignmentFileType::fasta);
  } else if (extension == "phylip" || extension == "phy") {
    return read_alignment(alignment_file, AlignmentFileType::phylip);
  }
  throw std::runtime_error{"Failed to recognize alignment file type"};
}
