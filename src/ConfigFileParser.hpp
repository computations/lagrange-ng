#pragma once

#include <cstdint>
#include <expected>
#include <logger.hpp>
#include <string>

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

  std::optional<double> start;
  std::optional<double> end;

  std::optional<std::vector<std::string>> include_areas;
  std::optional<std::vector<std::string>> exclude_areas;
};

using PeriodMap = std::unordered_map<std::string, PeriodConfig>;

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
    -> ParsingResult<std::string>;

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

auto determine_fossil_type(const std::string& fossil_type_string) -> FossilType;

auto parse_fossil(ConfigLexer& lexer) -> ParsingResult<Fossil>;

auto determine_output_file_type(const std::string& type_string) -> OutputType;

auto parse_opt_method(const std::string& method_string) -> OptimizationMethod;

template <typename T>
auto parse_assignment(ConfigLexer& lexer) -> ParsingResult<T> {
  if (auto r = lexer.expect(ConfigLexemeType::EQUALS_SIGN); !r) {
    return std::unexpected{ParsingError::expected_equals_sign};
  }

  return parse<T>(lexer);
};

auto parse_assignment(ConfigLexer& lexer,
                      ToLowerOption l = ToLowerOption::nolower)
    -> ParsingResult<std::string>;

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
                                                   PeriodMap& period_map);

ParsingResult<void> parse_period_exclude_statement(ConfigLexer& lexer,
                                                   PeriodMap& period_map);

ParsingResult<void> parse_period_start_statement(ConfigLexer& lexer,
                                                 PeriodMap& period_map);

ParsingResult<void> parse_period_end_statement(ConfigLexer& lexer,
                                               PeriodMap& period_map);

ParsingResult<void> parse_period_matrix_statement(ConfigLexer& lexer,
                                                  PeriodMap& period_map);

ParsingResult<void> parse_period_statement(ConfigLexer& lexer,
                                           PeriodMap& period_map);

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
};
}  // namespace lagrange
