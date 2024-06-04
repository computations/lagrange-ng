#include "ConfigFile.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

#include "Fossil.hpp"
#include "IO.hpp"
#include "logger.hpp"

namespace lagrange {
enum class ConfigLexemeType { VALUE, EQUALS_SIGN, END };

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

  auto consumeValueAsLowerString() -> std::string {
    auto str = consumeValueAsString();
    std::transform(str.cbegin(), str.cend(), str.begin(), [](char c) -> char {
      return std::tolower(c);
    });
    return str;
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

  auto consumeValueAsDist() -> Range {
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

std::string parse_filename(ConfigLexer& lexer) {
  lexer.expect(ConfigLexemeType::VALUE);
  auto filename_value = lexer.consumeValueAsString();
  return filename_value;
}

std::vector<std::string> parse_list(ConfigLexer& lexer) {
  std::vector<std::string> values;

  do {
    lexer.expect(ConfigLexemeType::VALUE);
    values.push_back(lexer.consumeValueAsString());
  } while (lexer.peak() != ConfigLexemeType::END);

  return values;
}

double parse_double(ConfigLexer& lexer) {
  lexer.expect(ConfigLexemeType::VALUE);

  return lexer.consumeValueAsFloat();
}

size_t parse_size_t(ConfigLexer& lexer) {
  lexer.expect(ConfigLexemeType::VALUE);

  return lexer.consumeValueAsSizeT();
}

FossilType determine_fossil_type(const std::string& fossil_type_string) {
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

Fossil parse_fossil(ConfigLexer& lexer) {
  lexer.expect(ConfigLexemeType::VALUE);
  auto fossil_type_string = lexer.consumeValueAsLowerString();

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

OutputType determine_output_file_type(const std::string& type_string) {
  if (type_string == "csv") { return OutputType::CSV; }
  if (type_string == "json") { return OutputType::JSON; }
  return OutputType::UNKNOWN;
}

ConfigFile ConfigFile::parse_config_file(std::istream& instream) {
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
        config._tree_filename = parse_filename(lexer);
      } else if (config_value == "datafile") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config._data_filename = parse_filename(lexer);
      } else if (config_value == "ratematrix") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config._rate_matrix_filename = parse_filename(lexer);
      } else if (config_value == "areanames") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config._area_names = parse_list(lexer);
      } else if (config_value == "prefix") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        config._prefix = parse_filename(lexer);
      } else if (config_value == "periods") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        auto tmp_values = parse_list(lexer);

        std::vector<double> time_points;
        time_points.resize(tmp_values.size());

        std::transform(
            tmp_values.begin(),
            tmp_values.end(),
            time_points.begin(),
            [](const std::string& s) -> double { return std::stof(s); });

        config._periods = Periods(time_points);
      } else if (config_value == "mrca") {
        lexer.expect(ConfigLexemeType::VALUE);
        auto mrca_name = lexer.consumeValueAsString();

        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        auto mrca_entries = parse_list(lexer);
        config._mrcas[mrca_name] = std::shared_ptr<MRCAEntry>(
            new MRCAEntry{.clade = std::move(mrca_entries)});
      } else if (config_value == "fossil") {
        config._fossils.push_back(parse_fossil(lexer));
      } else if (config_value == "calctype") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto calc_type_value = lexer.consumeValueAsString();
        if (std::tolower(calc_type_value[0]) == 'm') {
          config._marginal = true;
        } else {
          config._marginal = false;
        }
      } else if (config_value == "report") {
        lexer.expect(ConfigLexemeType::VALUE);
        auto report_value = lexer.consumeValueAsString();
        if (report_value != "split") { config._all_splits = false; }
      } else if (config_value == "splits") {
        if (lexer.peak() == ConfigLexemeType::VALUE) {
          config._split_nodes = parse_list(lexer);
        } else {
          config._all_splits = true;
        }
      } else if (config_value == "states") {
        if (lexer.peak() == ConfigLexemeType::VALUE) {
          config._state_nodes = parse_list(lexer);
        } else {
          config._all_states = true;
        }
      } else if (config_value == "dispersion") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config._extinction = parse_double(lexer);
      } else if (config_value == "extinction") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config._extinction = parse_double(lexer);
      } else if (config_value == "lh-epsilon") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config._extinction = parse_double(lexer);
      } else if (config_value == "workers") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config._workers = parse_size_t(lexer);
      } else if (config_value == "threads-per-worker") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config._threads_per_worker = parse_size_t(lexer);
      } else if (config_value == "maxareas") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);

        config.max_areas(parse_size_t(lexer));
      } else if (config_value == "expm-mode") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto expm_type = lexer.consumeValueAsString();

        if (expm_type == "krylov") {
          config._expm_mode = LagrangeEXPMComputationMode::KRYLOV;
        } else if (expm_type == "pade") {
          config._expm_mode = LagrangeEXPMComputationMode::PADE;
        } else if (expm_type == "adaptive") {
          config._expm_mode = LagrangeEXPMComputationMode::ADAPTIVE;
        }
      } else if (config_value == "mode") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto mode_type = lexer.consumeValueAsString();
        if (mode_type == "optimize") {
          config._run_mode = LagrangeOperationMode::OPTIMIZE;
        } else if (mode_type == "evaluate") {
          config._run_mode = LagrangeOperationMode::EVALUATE;
        }
      } else if (config_value == "alignment-type") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto align_type_str = lexer.consumeValueAsString();
        if (align_type_str == "fasta") {
          config._alignment_file_type = AlignmentFileType::FASTA;
        }
        if (align_type_str == "phylip") {
          config._alignment_file_type = AlignmentFileType::PHYLIP;
        }
      } else if (config_value == "logfile") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        config._log_filename = lexer.consumeValueAsString();
      } else if (config_value == "output-type") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto value = lexer.consumeValueAsLowerString();
        config.output_file_type(determine_output_file_type(value));
      } else {
        std::stringstream oss;
        oss << "Option '" << config_value << "' on line " << line_number
            << " was not recognized";
        throw ConfigFileParsingError{oss.str()};
      }

      lexer.expect(ConfigLexemeType::END);
    } catch (const std::exception& e) {
      std::ostringstream oss;
      oss << "There was a problem parsing line " << line_number
          << " of the config file:\n"
          << e.what();
      throw ConfigFileParsingError{oss.str()};
    }

    line_number++;
  }

  return config;
}

const std::filesystem::path& ConfigFile::tree_filename() const {
  return _tree_filename;
}

void ConfigFile::tree_filename(const std::filesystem::path& path) {
  _tree_filename = path;
}

const std::filesystem::path& ConfigFile::data_filename() const {
  return _data_filename;
}

void ConfigFile::data_filename(const std::filesystem::path& path) {
  _data_filename = path;
}

const std::filesystem::path& ConfigFile::log_filename() const {
  return _log_filename;
}

void ConfigFile::log_filename(const std::filesystem::path& path) {
  _log_filename = path;
}

const std::filesystem::path& ConfigFile::prefix() const { return _prefix; }

void ConfigFile::prefix(const std::filesystem::path& path) { _prefix = path; }

bool ConfigFile::has_rate_matrix_filename() const {
  return _rate_matrix_filename.hasValue();
}

OutputType ConfigFile::output_file_type() const {
  if (_output_file_type.hasValue()) { return _output_file_type.get(); }
  return OutputType::JSON;
}

void ConfigFile::output_file_type(OutputType type) { _output_file_type = type; }

const std::filesystem::path& ConfigFile::rate_matrix_filename() const {
  return _rate_matrix_filename.get();
}

void ConfigFile::rate_matrix_filename(const std::filesystem::path& path) {
  _rate_matrix_filename = path;
}

bool ConfigFile::has_max_areas() const { return _max_areas.hasValue(); }

size_t ConfigFile::max_areas() const { return _max_areas.get(); }

void ConfigFile::max_areas(size_t m) { _max_areas = m; }

const std::vector<std::string>& ConfigFile::area_names() const {
  return _area_names;
}

void ConfigFile::area_names(const std::vector<std::string>& area_names) {
  _area_names = area_names;
}

const Periods& ConfigFile::periods() const { return _periods; }

void ConfigFile::periods(const Periods& periods) { _periods = periods; }

const std::shared_ptr<MRCAEntry>& ConfigFile::mrca(
    const MRCALabel& label) const {
  return _mrcas.at(label);
}

const MRCAMap& ConfigFile::mrcas() const { return _mrcas; }

void ConfigFile::add_mrca(const MRCALabel& label,
                          const std::shared_ptr<MRCAEntry>& entry) {
  _mrcas[label] = entry;
}

const std::vector<Fossil>& ConfigFile::fossils() const { return _fossils; }

void ConfigFile::add_fossil(const Fossil& f) { _fossils.push_back(f); }

bool ConfigFile::marginal() const { return _marginal; }

void ConfigFile::marginal(bool m) { _marginal = m; }

bool ConfigFile::compute_all_splits() const { return _all_splits; }

void ConfigFile::compute_all_splits(bool b) { _all_splits = b; }

bool ConfigFile::compute_all_states() const { return _all_states; }

void ConfigFile::compute_all_states(bool b) { _all_states = b; }

const std::vector<MRCALabel>& ConfigFile::state_nodes() const {
  return _state_nodes;
}

void ConfigFile::state_nodes(const std::vector<MRCALabel>& labels) {
  _state_nodes = labels;
}

const std::vector<MRCALabel>& ConfigFile::split_nodes() const {
  return _split_nodes;
}

void ConfigFile::split_nodes(const std::vector<MRCALabel>& labels) {
  _split_nodes = labels;
}

PeriodParams ConfigFile::period_params() const {
  return {.dispersion_rate = _dispersion, .extinction_rate = _extinction};
}

void ConfigFile::period_params(double d, double e) {
  _dispersion = d;
  _extinction = e;
}

double ConfigFile::lh_epsilon() const { return _lh_epsilon; }

void ConfigFile::lh_epsilon(double e) { _lh_epsilon = e; }

bool ConfigFile::has_expm_mode() const { return _expm_mode.hasValue(); }

const LagrangeEXPMComputationMode& ConfigFile::expm_mode() const {
  return _expm_mode.get();
}

void ConfigFile::expm_mode(const LagrangeEXPMComputationMode& mode) {
  _expm_mode = mode;
}

AlignmentFileType ConfigFile::alignment_file_type() const {
  if (_alignment_file_type.hasValue()) { return _alignment_file_type.get(); }

  auto extension = data_filename().extension();
  if (extension == ".fasta" || extension == ".fas") {
    return AlignmentFileType::FASTA;
  } else if (extension == ".phylip" || extension == ".phy") {
    return AlignmentFileType::PHYLIP;
  }
  throw AlignmentFiletypeError{"Failed to recognize alignment file type"};
}

void ConfigFile::alignment_file_type(const AlignmentFileType& type) {
  _alignment_file_type = type;
}

size_t ConfigFile::region_count() const { return _region_count.get(); }

void ConfigFile::region_count(size_t r) { _region_count = r; }

size_t ConfigFile::workers() const { return _workers; }

void ConfigFile::workers(size_t w) { _workers = w; }

size_t ConfigFile::threads_per_worker() const { return _threads_per_worker; }

void ConfigFile::threads_per_worker(size_t t) { _threads_per_worker = t; }

const LagrangeOperationMode& ConfigFile::run_mode() const {
  return _run_mode.get();
}

void ConfigFile::run_mode(const LagrangeOperationMode& mode) {
  _run_mode = mode;
}

std::filesystem::path ConfigFile::jsonResultsFilename() const {
  auto results_filename = _prefix;
  results_filename += ".results.json";
  return results_filename;
}

std::filesystem::path ConfigFile::nodeTreeFilename() const {
  auto node_tree_filename = _prefix;
  node_tree_filename += ".nodes.tre";
  return node_tree_filename;
}

std::filesystem::path ConfigFile::scaledTreeFilename() const {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".scaled.tre";
  return scaled_tree_filename;
}

std::filesystem::path ConfigFile::splitsCSVResultsFilename() const {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".splits.csv";
  return scaled_tree_filename;
}

std::filesystem::path ConfigFile::statesCSVResultsFilename() const {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".states.csv";
  return scaled_tree_filename;
}

std::filesystem::path ConfigFile::periodsCSVResultsFilename() const {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".periods.csv";
  return scaled_tree_filename;
}

std::filesystem::path ConfigFile::distributionsCSVResultsFilename() const {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".distributions.csv";
  return scaled_tree_filename;
}

bool ConfigFile::computeStates() const {
  return _all_states || !_state_nodes.empty();
}

bool ConfigFile::computeSplits() const {
  return _all_splits || !_split_nodes.empty();
}
}  // namespace lagrange
