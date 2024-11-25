#include "ConfigFile.hpp"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>

#include "Fossil.hpp"

namespace lagrange {
enum class ConfigLexemeType : uint8_t { VALUE, EQUALS_SIGN, END };

class ConfigLexer {
 public:
  explicit ConfigLexer(std::string input) : _input{std::move(input)} {};

  auto consume() -> ConfigLexemeType {
    auto token = peak();
    if (token == ConfigLexemeType::VALUE) {
      std::stringstream builder;
      char quote_char = 0;
      while (char tmp = _input[_current_index]) {
        /* open the quote */
        if ((quote_char == 0) && (tmp == '"' || tmp == '\'')) {
          quote_char = tmp;
          _current_index++;
          continue;
        }

        /* close the quote */
        if ((quote_char != 0) && tmp == quote_char) {
          quote_char = 0;
          _current_index++;
          continue;
        }

        if ((quote_char == 0) && (isPunct(tmp) || (std::isspace(tmp) != 0))) {
          break;
        }
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
      return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
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

  [[nodiscard]] auto describePosition() const -> std::string {
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
  static auto isPunct(char c) -> bool { return c == '='; }

  auto consumeTokenPos() -> std::pair<ConfigLexemeType, size_t> {
    auto start_index = _current_index;
    auto token = consume();
    return {token, start_index};
  }

  void skipWhitespace() {
    while (_current_index < _input.size()) {
      if (std::isspace(_input[_current_index]) == 0) { break; }
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
  size_t _current_index{0};
};

auto parse_filename(ConfigLexer& lexer) -> std::string {
  lexer.expect(ConfigLexemeType::VALUE);
  auto filename_value = lexer.consumeValueAsString();
  return filename_value;
}

auto parse_list(ConfigLexer& lexer) -> std::vector<std::string> {
  std::vector<std::string> values;

  do {
    lexer.expect(ConfigLexemeType::VALUE);
    values.push_back(lexer.consumeValueAsString());
  } while (lexer.peak() != ConfigLexemeType::END);

  return values;
}

auto parse_set(ConfigLexer& lexer) -> std::unordered_set<std::string> {
  auto tmp = parse_list(lexer);
  return {tmp.begin(), tmp.end()};
}

auto parse_double(ConfigLexer& lexer) -> double {
  lexer.expect(ConfigLexemeType::VALUE);

  return lexer.consumeValueAsFloat();
}

auto parse_size_t(ConfigLexer& lexer) -> size_t {
  lexer.expect(ConfigLexemeType::VALUE);

  return lexer.consumeValueAsSizeT();
}

auto determine_fossil_type(const std::string& fossil_type_string)
    -> FossilType {
  if (fossil_type_string == "n" || fossil_type_string == "node") {
    return FossilType::NODE;
  }
  if (fossil_type_string == "b" || fossil_type_string == "branch") {
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

auto parse_fossil(ConfigLexer& lexer) -> Fossil {
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

auto determine_output_file_type(const std::string& type_string) -> OutputType {
  if (type_string == "csv") { return OutputType::CSV; }
  if (type_string == "json") { return OutputType::JSON; }
  return OutputType::UNKNOWN;
}

auto parse_opt_method(const std::string& method_string) -> OptimizationMethod {
  if (method_string == "nelder-mead") {
    return OptimizationMethod::NELDER_MEAD;
  }
  if (method_string == "bobyqa") { return OptimizationMethod::BOBYQA; }
  if (method_string == "cobyla") { return OptimizationMethod::COBYLA; }
  if (method_string == "bfgs") { return OptimizationMethod::BFGS; }
  if (method_string == "direct") { return OptimizationMethod::DIRECT; }
  if (method_string == "stogo") { return OptimizationMethod::STOGO; }

  return OptimizationMethod::UNKNOWN;
}

auto ConfigFile::parse_config_file(std::istream& instream) -> ConfigFile {
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
        config._mrcas[mrca_name] = std::make_shared<MRCAEntry>(
            MRCAEntry{.clade = std::move(mrca_entries)});
      } else if (config_value == "fossil") {
        config._fossils.push_back(parse_fossil(lexer));
      } else if (config_value == "calctype") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto calc_type_value = lexer.consumeValueAsString();
        config._marginal = std::tolower(calc_type_value[0]) == 'm';
      } else if (config_value == "report") {
        lexer.expect(ConfigLexemeType::VALUE);
        auto report_value = lexer.consumeValueAsString();
        if (report_value != "split") { config._all_splits = false; }
      } else if (config_value == "splits") {
        if (lexer.peak() == ConfigLexemeType::VALUE) {
          config._split_nodes = parse_set(lexer);
        } else {
          config._all_splits = true;
        }
      } else if (config_value == "states") {
        if (lexer.peak() == ConfigLexemeType::VALUE) {
          config._state_nodes = parse_set(lexer);
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

        config._lh_epsilon = parse_double(lexer);
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
      } else if (config_value == "opt-method") {
        lexer.expect(ConfigLexemeType::EQUALS_SIGN);
        lexer.expect(ConfigLexemeType::VALUE);

        auto value = lexer.consumeValueAsLowerString();
        config.opt_method(parse_opt_method(value));
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

auto ConfigFile::tree_filename() const -> const std::filesystem::path& {
  return _tree_filename;
}

void ConfigFile::tree_filename(const std::filesystem::path& path) {
  _tree_filename = path;
}

auto ConfigFile::data_filename() const -> const std::filesystem::path& {
  return _data_filename;
}

void ConfigFile::data_filename(const std::filesystem::path& path) {
  _data_filename = path;
}

auto ConfigFile::log_filename() const -> const std::filesystem::path& {
  return _log_filename;
}

void ConfigFile::log_filename(const std::filesystem::path& path) {
  _log_filename = path;
}

auto ConfigFile::prefix() const -> const std::filesystem::path& {
  return _prefix;
}

void ConfigFile::prefix(const std::filesystem::path& path) { _prefix = path; }

auto ConfigFile::has_rate_matrix_filename() const -> bool {
  return _rate_matrix_filename.hasValue();
}

auto ConfigFile::output_file_type() const -> OutputType {
  if (_output_file_type.hasValue()) { return _output_file_type.get(); }
  return OutputType::JSON;
}

void ConfigFile::output_file_type(OutputType type) { _output_file_type = type; }

auto ConfigFile::rate_matrix_filename() const -> const std::filesystem::path& {
  return _rate_matrix_filename.get();
}

void ConfigFile::rate_matrix_filename(const std::filesystem::path& path) {
  _rate_matrix_filename = path;
}

auto ConfigFile::has_max_areas() const -> bool { return _max_areas.hasValue(); }

auto ConfigFile::max_areas() const -> size_t { return _max_areas.get(); }

void ConfigFile::max_areas(size_t m) { _max_areas = m; }

auto ConfigFile::area_names() const -> const std::vector<std::string>& {
  return _area_names;
}

void ConfigFile::area_names(const std::vector<std::string>& area_names) {
  _area_names = area_names;
}

auto ConfigFile::periods() const -> const Periods& { return _periods; }

void ConfigFile::periods(const Periods& periods) { _periods = periods; }

auto ConfigFile::mrca(const MRCALabel& label) const
    -> const std::shared_ptr<MRCAEntry>& {
  return _mrcas.at(label);
}

auto ConfigFile::mrcas() const -> const MRCAMap& { return _mrcas; }

void ConfigFile::add_mrca(const MRCALabel& label,
                          const std::shared_ptr<MRCAEntry>& entry) {
  _mrcas[label] = entry;
}

auto ConfigFile::fossils() const -> const std::vector<Fossil>& {
  return _fossils;
}

void ConfigFile::add_fossil(const Fossil& f) { _fossils.push_back(f); }

auto ConfigFile::marginal() const -> bool { return _marginal; }

void ConfigFile::marginal(bool m) { _marginal = m; }

auto ConfigFile::compute_all_splits() const -> bool { return _all_splits; }

void ConfigFile::compute_all_splits(bool b) { _all_splits = b; }

auto ConfigFile::compute_all_states() const -> bool { return _all_states; }

void ConfigFile::compute_all_states(bool b) { _all_states = b; }

auto ConfigFile::state_nodes() const -> const std::unordered_set<MRCALabel>& {
  return _state_nodes;
}

void ConfigFile::state_nodes(const std::unordered_set<MRCALabel>& labels) {
  _state_nodes = labels;
}

auto ConfigFile::split_nodes() const -> const std::unordered_set<MRCALabel>& {
  return _split_nodes;
}

void ConfigFile::split_nodes(const std::unordered_set<MRCALabel>& labels) {
  _split_nodes = labels;
}

auto ConfigFile::period_params() const -> PeriodParams {
  return {.dispersion_rate = _dispersion, .extinction_rate = _extinction};
}

void ConfigFile::period_params(double d, double e) {
  _dispersion = d;
  _extinction = e;
}

auto ConfigFile::lh_epsilon() const -> double { return _lh_epsilon; }

void ConfigFile::lh_epsilon(double e) { _lh_epsilon = e; }

auto ConfigFile::has_expm_mode() const -> bool { return _expm_mode.hasValue(); }

auto ConfigFile::expm_mode() const -> const LagrangeEXPMComputationMode& {
  return _expm_mode.get();
}

void ConfigFile::expm_mode(const LagrangeEXPMComputationMode& mode) {
  _expm_mode = mode;
}

auto ConfigFile::alignment_file_type() const -> AlignmentFileType {
  if (_alignment_file_type.hasValue()) { return _alignment_file_type.get(); }

  auto extension = data_filename().extension();
  if (extension == ".fasta" || extension == ".fas") {
    return AlignmentFileType::FASTA;
  }
  if (extension == ".phylip" || extension == ".phy") {
    return AlignmentFileType::PHYLIP;
  }
  throw AlignmentFiletypeError{"Failed to recognize alignment file type"};
}

void ConfigFile::alignment_file_type(const AlignmentFileType& type) {
  _alignment_file_type = type;
}

auto ConfigFile::region_count() const -> size_t { return _region_count.get(); }

void ConfigFile::region_count(size_t r) { _region_count = r; }

auto ConfigFile::workers() const -> size_t {
  if (_workers.hasValue()) { return _workers.get(); }
  return 1;
}

void ConfigFile::workers(size_t w) { _workers = w; }

auto ConfigFile::threads_per_worker() const -> size_t {
  return static_cast<size_t>(_threads_per_worker);
}

void ConfigFile::threads_per_worker(size_t t) { _threads_per_worker = t; }

auto ConfigFile::run_mode() const -> const LagrangeOperationMode& {
  return _run_mode.get();
}

void ConfigFile::run_mode(const LagrangeOperationMode& mode) {
  _run_mode = mode;
}

auto ConfigFile::opt_method() const -> OptimizationMethod {
  return _opt_method.get();
}

void ConfigFile::opt_method(const OptimizationMethod& om) { _opt_method = om; }

auto ConfigFile::jsonResultsFilename() const -> std::filesystem::path {
  auto results_filename = _prefix;
  results_filename += ".results.json";
  return results_filename;
}

auto ConfigFile::nodeTreeFilename() const -> std::filesystem::path {
  auto node_tree_filename = _prefix;
  node_tree_filename += ".nodes.tre";
  return node_tree_filename;
}

auto ConfigFile::scaledTreeFilename() const -> std::filesystem::path {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".scaled.tre";
  return scaled_tree_filename;
}

auto ConfigFile::splitsTreeFilename() const -> std::filesystem::path {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".splits.tre";
  return scaled_tree_filename;
}

auto ConfigFile::statesTreeFilename() const -> std::filesystem::path {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".states.tre";
  return scaled_tree_filename;
}

auto ConfigFile::splitsCSVResultsFilename() const -> std::filesystem::path {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".splits.csv";
  return scaled_tree_filename;
}

auto ConfigFile::statesCSVResultsFilename() const -> std::filesystem::path {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".states.csv";
  return scaled_tree_filename;
}

auto ConfigFile::periodsCSVResultsFilename() const -> std::filesystem::path {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".periods.csv";
  return scaled_tree_filename;
}

auto ConfigFile::distributionsCSVResultsFilename() const
    -> std::filesystem::path {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".distributions.csv";
  return scaled_tree_filename;
}

auto ConfigFile::nodeInfoCSVResultsFilename() const -> std::filesystem::path {
  auto scaled_tree_filename = _prefix;
  scaled_tree_filename += ".node-info.csv";
  return scaled_tree_filename;
}

auto ConfigFile::computeStatesStrict() const -> bool {
  return _all_states || !_state_nodes.empty();
}

auto ConfigFile::computeStates() const -> bool {
  return _all_states || !_state_nodes.empty() || computeSplits();
}

auto ConfigFile::computeSplits() const -> bool {
  return _all_splits || !_split_nodes.empty();
}
}  // namespace lagrange
