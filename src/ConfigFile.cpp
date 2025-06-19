#include "ConfigFile.hpp"

#include <algorithm>
#include <cctype>
#include <expected>
#include <filesystem>
#include <format>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>

#include "Fossil.hpp"
#include "logger.hpp"

namespace lagrange {
using namespace std::string_literals;

template <typename T>
struct AssignmentValue {
  std::string id;
  T value;
};

enum class ConfigLexemeType : uint8_t { VALUE, EQUALS_SIGN, END };

enum class LexerError {
  expected_wrong_token,
  unknown_token,
  value_conversion_failed,
};

enum class ToLowerOption {
  lower,
  nolower,
};

template <typename T>
using LexerResult = std::expected<T, LexerError>;

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

  template <typename T>
  auto consume() noexcept -> LexerResult<T> {
    std::string tmp;
    std::swap(tmp, _value);
    return tmp;
  }

  template <typename T>
  auto consume(const std::function<LexerResult<T>(const std::string&)>&
                   conversion_func) noexcept -> LexerResult<T> {
    auto res = consume<std::string>();
    if (res) { return conversion_func(*res); }
    return std::unexpected{res.error()};
  }

  template <std::ranges::range T, typename U>
  auto consume(const std::function<U(U)>& trans_func) noexcept
      -> LexerResult<T> {
    auto res = consume<T>();
    if (!res) { return std::unexpected{res.error()}; }
    std::transform(res->cbegin(), res->cend(), res->begin(), trans_func);
    return res;
  }

  auto consumeAsRange() noexcept -> LexerResult<Range> {
    auto res = consume<std::string>();
    if (res) { return convert_dist_binary_string_to_dist(*res); }
    return std::unexpected{res.error()};
  }

  auto consumeValueAsLowerString() noexcept -> LexerResult<std::string> {
    return consume<std::string, char>([](char c) -> char {
      return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    });
  }

  [[nodiscard]] auto describePosition() const -> std::string {
    std::stringstream builder;
    builder << "position " << _current_index;
    return builder.str();
  }

  auto expect(ConfigLexemeType token_type) -> LexerResult<void> {
    auto ret = consumeTokenPos();
    if (ret.first != token_type) {
      LOG_DEBUG("Got the wrong token type at position {} was expecting {}",
                std::to_string(ret.second + 1),
                describeToken(token_type));
      return std::unexpected{LexerError::expected_wrong_token};
    }
    return {};
  }

  template <typename T>
  auto expectAndConsume(ConfigLexemeType token_type) -> LexerResult<T> {
    if (auto r = expect(token_type); !r) { return std::unexpected{r.error()}; }
    return consume<T>();
  }

  auto expectAndConsume(ConfigLexemeType token_type,
                        ToLowerOption l = ToLowerOption::nolower)
      -> LexerResult<std::string> {
    if (auto r = expect(token_type); !r) { return std::unexpected{r.error()}; }
    if (l == ToLowerOption::nolower) {
      return consume<std::string>();
    } else if (l == ToLowerOption::lower) {
      return consumeValueAsLowerString();
    } else {
      return std::unexpected{LexerError::value_conversion_failed};
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

template <>
auto ConfigLexer::consume<double>() noexcept -> LexerResult<double> {
  return consume<std::string>().and_then(
      [this](const std::string& f_str) -> LexerResult<double> {
        size_t pos = 0;
        double val = std::stod(f_str, &pos);
        if (pos != f_str.size()) {
          LOG_ERROR("Float conversion failed around {}", describePosition());
          return std::unexpected{LexerError::value_conversion_failed};
        }
        return val;
      });
}

template <>
auto ConfigLexer::consume<size_t>() noexcept -> LexerResult<size_t> {
  return consume<std::string>().and_then(
      [this](const std::string& d_str) -> LexerResult<size_t> {
        size_t pos = 0;
        size_t val = std::stoull(d_str, &pos);
        if (pos != d_str.size()) {
          LOG_ERROR("size_t conversion failed around {}", describePosition());
          return std::unexpected{LexerError::value_conversion_failed};
        }
        return val;
      });
}

template <>
auto ConfigLexer::consume<bool>() noexcept -> LexerResult<bool> {
  auto b_str = consumeValueAsLowerString();
  if (!b_str) { return std::unexpected{b_str.error()}; }

  if (b_str == "true" || b_str == "1") { return true; }
  if (b_str == "false" || b_str == "0") { return false; }

  return std::unexpected{LexerError::value_conversion_failed};
}

template <typename T>
concept VectorLike = requires(T& t, T::value_type i) {
  typename T::value_type;
  t.push_back(i);
} && !std::convertible_to<T, std::string>;

template <typename T>
concept SetLike = requires(T& t, T::value_type i) {
  typename T::value_type;
  t.insert(i);
};

template <typename T>
auto parse(ConfigLexer& lexer) -> ParsingResult<T> {
  if (auto r = lexer.expectAndConsume<T>(ConfigLexemeType::VALUE)) {
    return *r;
  } else if (r == std::unexpected{LexerError::value_conversion_failed}) {
    return std::unexpected{ParsingError::conversion_error};
  } else if (r == std::unexpected{LexerError::expected_wrong_token}) {
    return std::unexpected{ParsingError::expected_value};
  } else {
    return std::unexpected{ParsingError::lexing_error};
  }
}

auto parse(ConfigLexer& lexer, ToLowerOption l = ToLowerOption::nolower)
    -> ParsingResult<std::string> {
  if (auto r = lexer.expectAndConsume(ConfigLexemeType::VALUE, l)) {
    return *r;
  } else if (r == std::unexpected{LexerError::value_conversion_failed}) {
    return std::unexpected{ParsingError::conversion_error};
  } else if (r == std::unexpected{LexerError::expected_wrong_token}) {
    return std::unexpected{ParsingError::expected_value};
  } else {
    return std::unexpected{ParsingError::lexing_error};
  }
}

template <VectorLike T>
auto parse(ConfigLexer& lexer) -> ParsingResult<T> {
  T values;
  using U = T::value_type;

  do {
    if (auto r = parse<U>(lexer)) {
      values.push_back(*r);
    } else {
      return std::unexpected{r.error()};
    }
  } while (lexer.peak() != ConfigLexemeType::END);

  return values;
}

template <SetLike T>
auto parse(ConfigLexer& lexer) -> ParsingResult<T> {
  T values;
  using U = T::value_type;

  do {
    if (auto r = parse<U>(lexer)) {
      values.insert(*r);
    } else {
      return std::unexpected{r.error()};
    }
  } while (lexer.peak() != ConfigLexemeType::END);

  return values;
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

auto parse_fossil(ConfigLexer& lexer) -> ParsingResult<Fossil> {
  if (!lexer.expect(ConfigLexemeType::VALUE)) {
    return std::unexpected{ParsingError::expected_value};
  }
  auto fossil_type_string = *lexer.consumeValueAsLowerString();

  FossilType ft = determine_fossil_type(fossil_type_string);

  auto mrca_name = parse<std::string>(lexer);
  if (!mrca_name) { return std::unexpected{mrca_name.error()}; }

  if (!lexer.expect(ConfigLexemeType::EQUALS_SIGN)) {
    return std::unexpected{ParsingError::expected_equals_sign};
  }

  if (!lexer.expect(ConfigLexemeType::VALUE)) {
    return std::unexpected{ParsingError::expected_value};
  }
  auto fossil_area = *lexer.consumeAsRange();

  return Fossil{.mrca_name = *mrca_name,
                .clade = {},
                .age = ft == FossilType::BRANCH ? *parse<double>(lexer) : 0.0,
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

template <typename T>
auto parse_assignment(ConfigLexer& lexer) -> ParsingResult<T> {
  if (auto r = lexer.expect(ConfigLexemeType::EQUALS_SIGN); !r) {
    return std::unexpected{ParsingError::expected_equals_sign};
  }

  return parse<T>(lexer);
};

auto parse_assignment(ConfigLexer& lexer,
                      ToLowerOption l = ToLowerOption::nolower)
    -> ParsingResult<std::string> {
  if (auto r = lexer.expect(ConfigLexemeType::EQUALS_SIGN); !r) {
    return std::unexpected{ParsingError::expected_equals_sign};
  }

  return parse(lexer, l);
};

template <typename T>
auto parse_id_assignment(ConfigLexer& lexer)
    -> ParsingResult<AssignmentValue<T>> {
  auto id_name = parse<std::string>(lexer);
  if (!id_name) { return std::unexpected{id_name.error()}; }

  auto val = parse_assignment<T>(lexer);
  if (!val) { return std::unexpected{val.error()}; }

  return AssignmentValue{.id = *id_name, .value = *val};
};

template <typename T>
ParsingResult<AssignmentValue<T>> parse_and_check_period_map(
    ConfigLexer& lexer, const PeriodMap& map) {
  auto res = parse_id_assignment<T>(lexer);
  if (!res) { return res; }
  if (!map.contains(res->id)) {
    LOG_ERROR("The period '{}' was not declared before use at {}",
              res->id,
              lexer.describePosition());
    return std::unexpected{ParsingError::id_not_found};
  }
  return res;
}

ParsingResult<void> parse_period_include_statement(ConfigLexer& lexer,
                                                   PeriodMap& period_map) {
  auto res =
      parse_and_check_period_map<std::vector<std::string>>(lexer, period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).include_areas = list;
  return {};
}

ParsingResult<void> parse_period_exclude_statement(ConfigLexer& lexer,
                                                   PeriodMap& period_map) {
  auto res =
      parse_and_check_period_map<std::vector<std::string>>(lexer, period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).exclude_areas = list;
  return {};
}

ParsingResult<void> parse_period_start_statement(ConfigLexer& lexer,
                                                 PeriodMap& period_map) {
  auto res = parse_and_check_period_map<double>(lexer, period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).start = list;
  return {};
}

ParsingResult<void> parse_period_end_statement(ConfigLexer& lexer,
                                               PeriodMap& period_map) {
  auto res = parse_and_check_period_map<double>(lexer, period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).end = list;
  return {};
}

ParsingResult<void> parse_period_matrix_statement(ConfigLexer& lexer,
                                                  PeriodMap& period_map) {
  auto res =
      parse_and_check_period_map<std::filesystem::path>(lexer, period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).adjustment_matrix_filename = list;
  return {};
}

ParsingResult<void> parse_period_statement(ConfigLexer& lexer,
                                           PeriodMap& period_map) {
  auto value = parse<std::string>(lexer);
  if (!value) { return std::unexpected{value.error()}; }

  if (value == "include") {
    return parse_period_include_statement(lexer, period_map);
  } else if (value == "exclude") {
    return parse_period_exclude_statement(lexer, period_map);
  } else if (value == "start") {
    return parse_period_start_statement(lexer, period_map);
  } else if (value == "end") {
    return parse_period_end_statement(lexer, period_map);
  } else if (value == "matrix") {
    return parse_period_matrix_statement(lexer, period_map);
  } else {
    if (period_map.contains(*value)) {
      LOG_ERROR("Duplicate period declaration for period id '{}' at {}",
                *value,
                lexer.describePosition());
      return std::unexpected{ParsingError::duplicate_id};
    } else {
      period_map.insert({*value, {}});
      return {};
    }
  }

  if (!lexer.expect(ConfigLexemeType::END)) {
    return std::unexpected{ParsingError::expected_end};
  }
  return {};
}

template <typename T>
auto parse_and_assign(T& dst, ConfigLexer& lexer) -> ParsingResult<void> {
  if (auto r = lexer.expect(ConfigLexemeType::EQUALS_SIGN); !r) {
    return std::unexpected{ParsingError::expected_end};
  }
  if (auto r = parse<T>(lexer)) {
    dst = *r;
    return {};
  } else {
    return std::unexpected{r.error()};
  }
}

template <typename T>
auto parse_and_assign(std::optional<T>& dst, ConfigLexer& lexer)
    -> ParsingResult<void> {
  if (auto r = lexer.expect(ConfigLexemeType::EQUALS_SIGN); !r) {
    return std::unexpected{ParsingError::expected_end};
  }
  if (auto r = parse<T>(lexer)) {
    dst = *r;
    return {};
  } else {
    return std::unexpected{r.error()};
  }
}

auto ConfigFile::parse_line(ConfigLexer& lexer,
                            ConfigFile& config,
                            size_t line_number) -> ParsingResult<void> {
  if (lexer.peak() == ConfigLexemeType::END) { return {}; }

  auto config_value = parse<std::string>(lexer);

  if (config_value == "treefile") {
    if (auto r = parse_and_assign(config._tree_filename, lexer); !r) {
      return r;
    }
  } else if (config_value == "datafile") {
    if (auto r = parse_and_assign(config._data_filename, lexer); !r) {
      return r;
    }
  } else if (config_value == "areanames") {
    if (auto r = parse_and_assign(config._area_names, lexer); !r) { return r; }
  } else if (config_value == "prefix") {
    if (auto r = parse_and_assign(config._prefix, lexer); !r) { return r; }
  } else if (config_value == "periods") {
    std::vector<double> period_times;
    if (auto r = parse_and_assign(period_times, lexer)) {
      config._periods = Periods(period_times);
    } else {
      return r;
    }
  } else if (config_value == "mrca") {
    if (auto r = parse_id_assignment<std::vector<std::string>>(lexer)) {
      auto [mrca_name, mrca_entries] = *r;
      config._mrcas[mrca_name] = std::make_shared<MRCAEntry>(
          MRCAEntry{.clade = std::move(mrca_entries)});
    } else {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "fossil") {
    if (auto r = parse_fossil(lexer)) {
      config._fossils.push_back(*r);
    } else {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "splits") {
    if (lexer.peak() == ConfigLexemeType::VALUE) {
      if (auto r = parse<std::unordered_set<MRCALabel>>(lexer)) {
        config._split_nodes = *r;
      } else {
        return std::unexpected{r.error()};
      }
    } else {
      config._all_splits = true;
    }
  } else if (config_value == "states") {
    if (lexer.peak() == ConfigLexemeType::VALUE) {
      if (auto r = parse<std::unordered_set<MRCALabel>>(lexer)) {
        config._state_nodes = *r;
      } else {
        return std::unexpected{r.error()};
      }
    } else {
      config._all_states = true;
    }
  } else if (config_value == "dispersion") {
    if (auto r = parse_and_assign(config._dispersion, lexer); !r) { return r; }
  } else if (config_value == "extinction") {
    if (auto r = parse_and_assign(config._extinction, lexer); !r) { return r; }
  } else if (config_value == "lh-epsilon") {
    if (auto r = parse_and_assign(config._lh_epsilon, lexer); !r) { return r; }
  } else if (config_value == "workers") {
    if (auto r = parse_and_assign(config._workers, lexer); !r) { return r; }
  } else if (config_value == "threads-per-worker") {
    if (auto r = parse_and_assign(config._threads_per_worker, lexer); !r) {
      return r;
    }
  } else if (config_value == "maxareas") {
    if (auto r = parse_and_assign(config._max_areas, lexer); !r) { return r; }
  } else if (config_value == "expm-mode") {
    auto expm_type = parse_assignment<std::string>(lexer);

    if (expm_type == "krylov") {
      config._expm_mode = LagrangeEXPMComputationMode::KRYLOV;
    } else if (expm_type == "pade") {
      config._expm_mode = LagrangeEXPMComputationMode::PADE;
    } else if (expm_type == "adaptive") {
      config._expm_mode = LagrangeEXPMComputationMode::ADAPTIVE;
    } else if (expm_type) {
      LOG_ERROR("Unknown expm type '{}'", *expm_type);
      return std::unexpected{ParsingError::unknown_expm_type};
    } else {
      return std::unexpected{expm_type.error()};
    }
  } else if (config_value == "mode") {
    auto mode_type = parse_assignment<std::string>(lexer);
    if (mode_type == "optimize") {
      config._run_mode = LagrangeOperationMode::OPTIMIZE;
    } else if (mode_type == "evaluate") {
      config._run_mode = LagrangeOperationMode::EVALUATE;
    } else if (mode_type) {
      LOG_ERROR("Unknown run mode type '{}'", *mode_type);
      return std::unexpected{ParsingError::unknown_expm_type};
    } else {
      return std::unexpected{mode_type.error()};
    }
  } else if (config_value == "alignment-type") {
    auto align_type_str = parse_assignment<std::string>(lexer);
    if (align_type_str == "fasta") {
      config._alignment_file_type = AlignmentFileType::FASTA;
    } else if (align_type_str == "phylip") {
      config._alignment_file_type = AlignmentFileType::PHYLIP;
    } else if (align_type_str) {
      LOG_ERROR("Unknown alignment type '{}'", *align_type_str);
    } else {
      return std::unexpected{align_type_str.error()};
    }
  } else if (config_value == "logfile") {
    if (auto r = parse_and_assign(config._log_filename, lexer); r) {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "output-type") {
    if (auto r = parse_assignment(lexer, ToLowerOption::lower)) {
      config.output_file_type(determine_output_file_type(*r));
    } else {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "opt-method") {
    if (auto r = parse_assignment(lexer, ToLowerOption::lower)) {
      config.opt_method(parse_opt_method(*r));
    } else {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "allow-ambiguous") {
    if (lexer.peak() == ConfigLexemeType::EQUALS_SIGN) {
      if (auto r = parse_assignment<bool>(lexer)) {
        config.allow_ambigious(*r);
      } else {
        return std::unexpected{r.error()};
      }
    } else {
      config.allow_ambigious(true);
    }
  } else if (config_value == "dump-graph") {
    if (lexer.peak() == ConfigLexemeType::EQUALS_SIGN) {
      if (auto r = parse_assignment<bool>(lexer)) {
        config.dump_graph(*r);
      } else {
        return std::unexpected{r.error()};
      }
    } else {
      config.dump_graph(true);
    }
  } else if (config_value == "lwr-threshold") {
    if (auto r = parse_and_assign(config._lh_epsilon, lexer); !r) {
      return std::unexpected{r.error()};
    }
  } else if (config_value) {
    LOG_ERROR("Option '{}' on line {} was not recognized",
              *config_value,
              line_number);
  } else {
    LOG_ERROR("Failed to parse the option on line {}", line_number);
  }

  if (auto r = lexer.expect(ConfigLexemeType::END); !r) {
    return std::unexpected{ParsingError::expected_end};
  }

  return {};
}

auto ConfigFile::parse_config_file(std::istream& instream) -> ConfigFile {
  ConfigFile config;
  std::string line;
  size_t line_number = 1;
  while (getline(instream, line)) {
    ConfigLexer lexer(line);
    try {
      if (auto r = parse_line(lexer, config, line_number); !r) {
        LOG_ERROR("Failed to parse config on line '{}'", line_number);
      }
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
  return _adjustment_matrix_filename.has_value();
}

auto ConfigFile::output_file_type() const -> OutputType {
  if (_output_file_type.has_value()) { return _output_file_type.value(); }
  return OutputType::JSON;
}

void ConfigFile::output_file_type(OutputType type) { _output_file_type = type; }

auto ConfigFile::rate_matrix_filename() const -> const std::filesystem::path& {
  return _adjustment_matrix_filename.value();
}

void ConfigFile::rate_matrix_filename(const std::filesystem::path& path) {
  _adjustment_matrix_filename = path;
}

auto ConfigFile::has_max_areas() const -> bool {
  return _max_areas.has_value();
}

auto ConfigFile::max_areas() const -> size_t { return _max_areas.value(); }

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
  return {.dispersion_rate = _dispersion,
          .extinction_rate = _extinction,
          .adjustment_matrix = nullptr,
          .regions = *_region_count};
}

void ConfigFile::period_params(double d, double e) {
  _dispersion = d;
  _extinction = e;
}

auto ConfigFile::lh_epsilon() const -> double { return _lh_epsilon; }

void ConfigFile::lh_epsilon(double e) { _lh_epsilon = e; }

auto ConfigFile::has_expm_mode() const -> bool {
  return _expm_mode.has_value();
}

auto ConfigFile::expm_mode() const -> const LagrangeEXPMComputationMode& {
  return _expm_mode.value();
}

void ConfigFile::expm_mode(const LagrangeEXPMComputationMode& mode) {
  _expm_mode = mode;
}

auto ConfigFile::alignment_file_type() const -> AlignmentFileType {
  if (_alignment_file_type.has_value()) { return _alignment_file_type.value(); }

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

auto ConfigFile::region_count() const -> size_t {
  return _region_count.value();
}

void ConfigFile::region_count(size_t r) {
  _region_count = r;
  if (!has_max_areas()) { max_areas(_region_count.value()); }
}

auto ConfigFile::workers() const -> size_t {
  if (_workers.has_value()) { return _workers.value(); }
  return 1;
}

void ConfigFile::workers(size_t w) { _workers = w; }

auto ConfigFile::threads_per_worker() const -> size_t {
  return static_cast<size_t>(*_threads_per_worker);
}

void ConfigFile::threads_per_worker(size_t t) { _threads_per_worker = t; }

auto ConfigFile::run_mode() const -> const LagrangeOperationMode& {
  return _run_mode.value();
}

void ConfigFile::run_mode(const LagrangeOperationMode& mode) {
  _run_mode = mode;
}

auto ConfigFile::opt_method() const -> OptimizationMethod {
  return _opt_method.value();
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

auto ConfigFile::cleanTreeFilename() const -> std::filesystem::path {
  auto node_tree_filename = _prefix;
  node_tree_filename += ".clean.tre";
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

auto ConfigFile::forwardGraphFilename() const -> std::filesystem::path {
  auto graphvis_filename = _prefix;
  graphvis_filename += ".forward-graph.gv";
  return graphvis_filename;
}

auto ConfigFile::reverseGraphFilename() const -> std::filesystem::path {
  auto graphvis_filename = _prefix;
  graphvis_filename += ".reverse-graph.gv";
  return graphvis_filename;
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

auto ConfigFile::allow_ambigious() const -> bool {
  return _allow_ambigious.value_or(true);
}

void ConfigFile::allow_ambigious(bool b) { _allow_ambigious = b; }

auto ConfigFile::dump_graph() const -> bool {
  return _dump_graph.value_or(false);
}

void ConfigFile::dump_graph(bool b) { _dump_graph = b; }

auto ConfigFile::lwrOutputThreshold() const -> double {
  return _output_threshold;
}

void ConfigFile::lwrOutputThreshold(double t) { _output_threshold = t; }
}  // namespace lagrange
