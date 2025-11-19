#pragma once

#include <cstdint>
#include <expected>
#include <functional>
#include <sstream>
#include <string>

#include "Utils.hpp"

namespace lagrange {
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

enum class ConfigLexemeType : uint8_t { VALUE, EQUALS_SIGN, COMMENT, END };

class ConfigLexer {
 public:
  explicit ConfigLexer(std::string input, size_t line_number = 0) :
      _input{std::move(input)},
      _line_number{line_number} {};

  auto readToken() -> ConfigLexemeType {
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
    if (_current_index >= _input.size()) { return ConfigLexemeType::END; }

    size_t tmp_index = _current_index;
    char current_char = _input[tmp_index++];

    if (current_char == '=') { return ConfigLexemeType::EQUALS_SIGN; }
    if (current_char == '#') { return ConfigLexemeType::COMMENT; }
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
    return std::format("line {}, position {}", _line_number, _current_index);
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
    while (token_type != readToken()) {}
  }

  auto atEnd() -> bool { return _input.size() == _current_index; }

 private:
  static auto isPunct(char c) -> bool { return c == '='; }

  auto consumeTokenPos() -> std::pair<ConfigLexemeType, size_t> {
    auto start_index = _current_index;
    auto token = readToken();
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
  size_t _line_number;
};

template <>
auto ConfigLexer::consume<double>() noexcept -> LexerResult<double>;

template <>
auto ConfigLexer::consume<size_t>() noexcept -> LexerResult<size_t>;

template <>
auto ConfigLexer::consume<bool>() noexcept -> LexerResult<bool>;
};  // namespace lagrange
