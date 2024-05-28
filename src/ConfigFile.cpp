#include "ConfigFile.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

#include "Fossil.hpp"
#include "logger.hpp"

namespace lagrange {
enum class ConfigLexemeType { VALUE, EQUALS_SIGN, END };

void set_mrcas_for_fossils(ConfigFile &config) {
  for (auto &f : config.fossils) { f.clade = config.mrcas.at(f.mrca_name); }
}

void check_prefix(ConfigFile &config) {
  if (config.prefix.empty()) { config.prefix = config.tree_filename; }
  if (config.prefix.filename().string()[0] == '.') {
    MESSAGE(WARNING,
            "The current prefix starts with a dot. The results will be hidden "
            "on unix-like systems");
  }
}

void finalize(ConfigFile &config) {
  set_mrcas_for_fossils(config);
  check_prefix(config);
}

class ConfigLexer {
 public:
  explicit ConfigLexer(std::string input) :
      _input{std::move(input)},
      _current_index{0} {};

  auto consume() -> ConfigLexemeType {
    auto token = peak();
    if (token == ConfigLexemeType::VALUE) {
      std::stringstream builder;
      char quote_char = 0;
      while (char tmp = _input[_current_index]) {
        /* open the quote */
        if (!quote_char && (tmp == '"' || tmp == '\'')) {
          quote_char = tmp;
          _current_index++;
          continue;
        }

        /* close the quote */
        if (quote_char && tmp == quote_char) {
          quote_char = 0;
          _current_index++;
          continue;
        }

        if (!quote_char && (isPunct(tmp) || std::isspace(tmp))) { break; }
        builder << tmp;
        _current_index++;
      }

      _value = builder.str();
      skipWhitespace();
      return token;
    }
    _current_index++;
    skipWhitespace();
    return token;
  }

  auto peak() -> ConfigLexemeType {
    size_t tmp_index = _current_index;
    char current_char = _input[tmp_index++];

    if (_current_index == _input.size()) { return ConfigLexemeType::END; }
    if (current_char == '=') { return ConfigLexemeType::EQUALS_SIGN; }
    return ConfigLexemeType::VALUE;
  }

  auto consumeValueAsString() -> std::string {
    std::string tmp;
    std::swap(tmp, _value);
    return tmp;
  }

  auto consumeValueAsFloat() -> double {
    auto f_str = consumeValueAsString();
    size_t pos = 0;
    double val = std::stod(f_str, &pos);
    if (pos != f_str.size()) {
      throw ConfigFileLexingError{std::string("Float conversion failed around")
                                  + describePosition()};
    }
    return val;
  }

  auto consumeValueAsDist() -> Dist {
    return convert_dist_binary_string_to_dist(consumeValueAsString());
  }

  auto consumeValueAsSizeT() -> size_t {
    auto f_str = consumeValueAsString();
    size_t pos = 0;
    size_t val = std::stoull(f_str, &pos);
    if (pos != f_str.size()) {
      throw ConfigFileLexingError{std::string("Float conversion failed around")
                                  + describePosition()};
    }
    return val;
  }

  auto describePosition() const -> std::string {
    std::stringstream builder;
    builder << "position " << _current_index;
    return builder.str();
  }

  void expect(ConfigLexemeType token_type) {
    auto ret = consumeTokenPos();
    if (ret.first != token_type) {
      throw ConfigFileLexingError{
          std::string("Got the wrong token type at position ")
          + std::to_string(ret.second + 1) + " was expecting "
          + describeToken(token_type)};
    }
  }

  void consumeUntil(ConfigLexemeType token_type) {
    while (token_type != consume()) {}
  }

  auto atEnd() -> bool { return _input.size() == _current_index; }

 private:
  bool isPunct(char c) { return c == '='; }

  auto consumeTokenPos() -> std::pair<ConfigLexemeType, size_t> {
    auto start_index = _current_index;
    auto token = consume();
    return {token, start_index};
  }

  void skipWhitespace() {
    while (_current_index < _input.size()) {
      char c = _input[_current_index];
      if (std::isspace(c) == 0) { break; }
      _current_index++;
    }
  }

  static auto describeToken(ConfigLexemeType token_type) -> std::string {
    switch (token_type) {
      case ConfigLexemeType::VALUE:
        return {"either a number or a string"};
      case ConfigLexemeType::EQUALS_SIGN:
        return {"equals sign"};
      case ConfigLexemeType::END:
        return {"end of line"};
      default:
        return {"unknown token"};
    }
  }

  std::string _input;
  std::string _value;
  size_t _current_index;
};

std::string parse_filename(ConfigLexer &lexer) {
  lexer.expect(ConfigLexemeType::VALUE);
  auto filename_value = lexer.consumeValueAsString();
  return filename_value;
}

std::vector<std::string> parse_list(ConfigLexer &lexer) {
  std::vector<std::string> values;

  do {
    lexer.expect(ConfigLexemeType::VALUE);
    values.push_back(lexer.consumeValueAsString());
  } while (lexer.peak() != ConfigLexemeType::END);

  return values;
}

double parse_double(ConfigLexer &lexer) {
  lexer.expect(ConfigLexemeType::VALUE);

  return lexer.consumeValueAsFloat();
}

size_t parse_size_t(ConfigLexer &lexer) {
  lexer.expect(ConfigLexemeType::VALUE);

  return lexer.consumeValueAsSizeT();
}

FossilType determine_fossil_type(std::string fossil_type_string) {
  std::transform(fossil_type_string.cbegin(),
                 fossil_type_string.cend(),
                 fossil_type_string.begin(),
                 [](char c) -> char { return std::tolower(c); });

  if (fossil_type_string == "n" || fossil_type_string == "node") {
    return FossilType::NODE;
  } else if (fossil_type_string == "b" || fossil_type_string == "branch") {
    return FossilType::BRANCH;
  } else if (fossil_type_string == "f" || fossil_type_string == "fixed") {
    return FossilType::FIXED;
  } else if (fossil_type_string == "i" || fossil_type_string == "include") {
    return FossilType::INCLUDE;
  } else if (fossil_type_string == "e" || fossil_type_string == "exclude") {
    return FossilType::EXCLUDE;
  }
  return FossilType::UNKOWN;
}

Fossil parse_fossil(ConfigLexer &lexer) {
  lexer.expect(ConfigLexemeType::VALUE);
  auto fossil_type_string = lexer.consumeValueAsString();

  FossilType ft = determine_fossil_type(fossil_type_string);

  lexer.expect(ConfigLexemeType::VALUE);
  auto mrca_name = lexer.consumeValueAsString();

  lexer.expect(ConfigLexemeType::EQUALS_SIGN);

  lexer.expect(ConfigLexemeType::VALUE);

  auto fossil_area = lexer.consumeValueAsDist();

  return {.mrca_name = mrca_name,
          .clade = {},
          .age = ft == FossilType::BRANCH ? parse_double(lexer) : 0.0,
          .area = fossil_area,
          .type = ft};
}

ConfigFile parse_config_file(std::istream &instream) {
  ConfigFile config;
  std::string line;
  size_t line_number = 1;
  while (getline(instream, line)) {
    ConfigLexer lexer(line);
    try {
      if (lexer.peak() == ConfigLexemeType::END) { continue; }

      lexer.expect(ConfigLexemeType::VALUE);
      auto config_value = lexer.consumeValueAsString();

      if (config_value == "treefile") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config.tree_filename = parse_filename(lexer);
      } else if (config_value == "datafile") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config.data_filename = parse_filename(lexer);
      } else if (config_value == "ratematrix") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config.rate_matrix_filename = parse_filename(lexer);
      } else if (config_value == "areanames") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config.area_names = parse_list(lexer);
      } else if (config_value == "prefix") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config.prefix = parse_filename(lexer);
      } else if (config_value == "periods") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        auto tmp_values = parse_list(lexer);

        std::vector<double> time_points;
        time_points.resize(tmp_values.size());

        std::transform(
            tmp_values.begin(),
            tmp_values.end(),
            time_points.begin(),
            [](const std::string &s) -> double { return std::stof(s); });

        config.periods = Periods(time_points);
      } else if (config_value == "mrca") {
        lexer.expect(ConfigLexemeType::VALUE);
        auto mrca_name = lexer.consumeValueAsString();

        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        auto mrca_entries = parse_list(lexer);
        config.mrcas[mrca_name] = std::shared_ptr<MRCAEntry>(
            new MRCAEntry{.clade = std::move(mrca_entries)});
      } else if (config_value == "ancstate") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config.anc_states = parse_list(lexer);
      } else if (config_value == "fossil") {
        config.fossils.push_back(parse_fossil(lexer));
      } else if (config_value == "calctype") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto calc_type_value = lexer.consumeValueAsString();
        if (std::tolower(calc_type_value[0]) == 'm') {
          config.marginal = true;
        } else {
          config.marginal = false;
        }
      } else if (config_value == "report") {
        lexer.expect(ConfigLexemeType::VALUE);
        auto report_value = lexer.consumeValueAsString();
        if (report_value != "split") { config.all_splits = false; }
      } else if (config_value == "splits") {
        if (lexer.peak() == ConfigLexemeType::VALUE) {
          config.split_nodes = parse_list(lexer);
        } else {
          config.all_splits = true;
        }
      } else if (config_value == "states") {
        if (lexer.peak() == ConfigLexemeType::VALUE) {
          config.state_nodes = parse_list(lexer);
        } else {
          config.all_states = true;
        }
      } else if (config_value == "dispersion") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config.extinction = parse_double(lexer);
      } else if (config_value == "extinction") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config.extinction = parse_double(lexer);
      } else if (config_value == "lh-epsilon") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config.extinction = parse_double(lexer);
      } else if (config_value == "workers") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config.workers = parse_size_t(lexer);
      } else if (config_value == "threads-per-worker") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config.threads_per_worker = parse_size_t(lexer);
      } else if (config_value == "maxareas") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config.maxareas = parse_size_t(lexer);
      } else if (config_value == "expm-mode") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto expm_type = lexer.consumeValueAsString();

        if (expm_type == "krylov") {
          config.expm_mode = LagrangeEXPMComputationMode::KRYLOV;
        } else if (expm_type == "pade") {
          config.expm_mode = LagrangeEXPMComputationMode::PADE;
        } else if (expm_type == "adaptive") {
          config.expm_mode = LagrangeEXPMComputationMode::ADAPTIVE;
        }
      } else if (config_value == "mode") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto mode_type = lexer.consumeValueAsString();
        if (mode_type == "optimize") {
          config.run_mode = LagrangeOperationMode::OPTIMIZE;
        } else if (mode_type == "evaluate") {
          config.run_mode = LagrangeOperationMode::EVALUATE;
        }
      } else if (config_value == "alignment-type") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto align_type_str = lexer.consumeValueAsString();
        if (align_type_str == "fasta") {
          config.alignment_file_type = AlignmentFileType::FASTA;
        }
        if (align_type_str == "phylip") {
          config.alignment_file_type = AlignmentFileType::PHYLIP;
        }
      } else if (config_value == "logfile") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        config.log_filename = lexer.consumeValueAsString();
      } else {
        std::stringstream oss;
        oss << "Option '" << config_value << "' on line " << line_number
            << " was not recognized";
        throw ConfigFileParsingError{oss.str()};
      }

      lexer.expect(ConfigLexemeType::END);
    } catch (const std::exception &e) {
      std::ostringstream oss;
      oss << "There was a problem parsing line " << line_number
          << " of the config file:\n"
          << e.what();
      throw ConfigFileParsingError{oss.str()};
    }

    line_number++;
  }

  finalize(config);

  return config;
}

std::filesystem::path ConfigFile::resultsFilename() const {
  auto results_filename = prefix;
  results_filename += ".results.json";
  return results_filename;
}

std::filesystem::path ConfigFile::NodeTreeFilename() const {
  auto node_tree_filename = prefix;
  node_tree_filename += ".nodes.tre";
  return node_tree_filename;
}

std::filesystem::path ConfigFile::scaledTreeFilename() const {
  auto scaled_tree_filename = prefix;
  scaled_tree_filename += ".scaled.tre";
  return scaled_tree_filename;
}

bool ConfigFile::computeStates() const {
  return all_states || !state_nodes.empty();
}

bool ConfigFile::computeSplits() const {
  return all_splits || !split_nodes.empty();
}
}  // namespace lagrange
