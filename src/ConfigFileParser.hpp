#pragma once

#include <concepts>
#include <cstdint>
#include <expected>
#include <logger.hpp>
#include <map>
#include <string>
#include <string_view>

#include "AdjustmentMatrix.hpp"
#include "ConfigFileLexer.hpp"
#include "Fossil.hpp"

namespace lagrange {
using namespace std::string_literals;

template <typename T>
struct AssignmentValue {
  std::string id;
  T value;
};

constexpr size_t KRYLOV_RANGE_COUNT_THRESHOLD = 5;

enum class OutputType : uint8_t { JSON, CSV, UNKNOWN };

enum class OptimizationMethod : uint8_t {
  NELDER_MEAD,
  BOBYQA,
  COBYLA,
  BFGS,
  DIRECT,
  STOGO,
  UNKNOWN,
};

enum class ParsingError {
  id_not_found,
  duplicate_id,
  lexing_error,
  expected_value,
  expected_equals_sign,
  expected_end,
  conversion_error,
  unknown_expm_type,
  invalid_option,
  unknown_error,
};

struct PeriodConfig {
  std::optional<double> dispersion;
  std::optional<double> extinction;

  std::optional<std::filesystem::path> adjustment_matrix_filename;
  std::optional<AdjustmentMatrix> adjustment_matrix;

  double start = 0.0;
  double end = std::numeric_limits<double>::infinity();

  std::optional<std::string> include_areas;
  std::optional<std::string> exclude_areas;
};

using PeriodConfigMap = std::unordered_map<std::string, PeriodConfig>;

template <typename T>
using ParsingResult = std::expected<T, ParsingError>;

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
concept IdMapLike = SetLike<T> && requires(T& t) {
  typename T::key_type;
  requires std::same_as<typename T::key_type, std::string>;
};

class ConfigFileParsingError : public std::runtime_error {
 public:
  ConfigFileParsingError(const std::string& msg) : std::runtime_error{msg} {}
};

class ConfigFile;

class ConfigFileParser {
 public:
  ConfigFileParser() = delete;

  ConfigFileParser(ConfigLexer&& lexer) : _lexer{lexer} {}

  ConfigFileParser(const std::string& line, size_t line_number) :
      _lexer{line, line_number} {}

  auto parse_line(ConfigFile& config) -> ParsingResult<void>;

  template <typename T>
  auto parse() -> ParsingResult<T> {
    if (auto r = _lexer.expectAndConsume<T>(ConfigLexemeType::VALUE)) {
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
  auto parse() -> ParsingResult<T> {
    T values;
    using U = T::value_type;

    do {
      if (auto r = parse<U>()) {
        values.push_back(*r);
      } else {
        return std::unexpected{r.error()};
      }
    } while (_lexer.peak() != ConfigLexemeType::END);

    return values;
  }

  template <SetLike T>
  auto parse() -> ParsingResult<T> {
    T values;
    using U = T::value_type;

    do {
      if (auto r = parse<U>()) {
        values.insert(*r);
      } else {
        return std::unexpected{r.error()};
      }
    } while (_lexer.peak() != ConfigLexemeType::END);

    return values;
  }

  auto parse(ToLowerOption l = ToLowerOption::nolower)
      -> ParsingResult<std::string>;

  auto parse_fossil() -> ParsingResult<Fossil>;

  auto parse_assignment(ToLowerOption l = ToLowerOption::nolower)
      -> ParsingResult<std::string>;

  template <typename T>
  auto parse_assignment() -> ParsingResult<T> {
    if (auto r = _lexer.expect(ConfigLexemeType::EQUALS_SIGN); !r) {
      return std::unexpected{ParsingError::expected_equals_sign};
    }

    return parse<T>();
  };

  template <typename T>
  auto parse_id_assignment() -> ParsingResult<AssignmentValue<T>> {
    auto id_name = parse<std::string>();
    if (!id_name) { return std::unexpected{id_name.error()}; }

    auto val = parse_assignment<T>();
    if (!val) { return std::unexpected{val.error()}; }

    return AssignmentValue{.id = *id_name, .value = *val};
  };

  template <typename T>
  ParsingResult<AssignmentValue<T>> parse_and_check_id_map(
      const IdMapLike auto& map) {
    auto res = parse_id_assignment<T>();
    if (!res) { return res; }
    if (!map.contains(res->id)) {
      LOG_ERROR("The id '{}' was not declared before use at {}",
                res->id,
                _lexer.describePosition());
      return std::unexpected{ParsingError::id_not_found};
    }
    return res;
  }

  template <typename T>
  auto parse_and_assign(T& dst) -> ParsingResult<void> {
    if (auto r = parse_assignment<T>()) {
      dst = *r;
      return {};
    } else {
      return std::unexpected{r.error()};
    }
  }

  template <typename T>
  auto parse_and_assign(std::optional<T>& dst) -> ParsingResult<void> {
    if (auto r = parse_assignment<T>()) {
      dst = *r;
      return {};
    } else {
      return std::unexpected{r.error()};
    }
  }

  ParsingResult<void> parse_period_include_statement(PeriodConfig& period);

  ParsingResult<void> parse_period_exclude_statement(PeriodConfig& period);

  ParsingResult<void> parse_period_start_statement(PeriodConfig& period);

  ParsingResult<void> parse_period_end_statement(PeriodConfig& period);

  ParsingResult<void> parse_period_matrix_statement(PeriodConfig& period);

  ParsingResult<void> parse_period_sub_statement(PeriodConfig& period);

  ParsingResult<void> parse_period_statement(PeriodConfigMap& period_map);

  auto describePosition() const -> std::string;

  static void print_help_long();
  static void print_help_short();

 private:
  static auto determine_fossil_type(const std::string& fossil_type_string)
      -> FossilType;
  static auto determine_output_file_type(const std::string& type_string)
      -> OutputType;
  static auto determine_opt_method(const std::string& method_string)
      -> OptimizationMethod;

  struct ConfigActionMapItem {
    std::function<ParsingResult<void>(ConfigFileParser&, ConfigFile&)> _action;
    std::string_view _name;
    std::string_view _help;
    std::string_view _usage;
    std::vector<std::string_view> _examples;
  };

  using ActionMapType = std::map<std::string_view, ConfigActionMapItem>;

  static ActionMapType _config_action_map;

  ConfigLexer _lexer;
};

}  // namespace lagrange
