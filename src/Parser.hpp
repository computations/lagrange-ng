#pragma once

#include <cstdint>
#include <expected>
#include <functional>
#include <logger.hpp>
#include <sstream>
#include <string>

#include "Fossil.hpp"
#include "Utils.hpp"

namespace lagrange {
using namespace std::string_literals;

template <typename T>
struct AssignmentValue {
  std::string id;
  T value;
};

constexpr size_t KRYLOV_RANGE_COUNT_THRESHOLD = 5;

enum class ConfigLexemeType : uint8_t { VALUE, EQUALS_SIGN, END };

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

enum class LexerError {
  expected_wrong_token,
  unknown_token,
  value_conversion_failed,
};

enum class ToLowerOption {
  lower,
  nolower,
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

template <>
auto ConfigLexer::consume<double>() noexcept -> LexerResult<double>;

template <>
auto ConfigLexer::consume<size_t>() noexcept -> LexerResult<size_t>;

template <>
auto ConfigLexer::consume<bool>() noexcept -> LexerResult<bool>;

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
