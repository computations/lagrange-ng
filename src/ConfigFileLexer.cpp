#include "ConfigFileLexer.hpp"

namespace lagrange{
template <>
auto ConfigLexer::consume<double>() noexcept -> LexerResult<double> {
  return consume<std::string>().and_then(
      [](const std::string& f_str) -> LexerResult<double> {
        try {
          size_t pos = 0;
          double val = std::stod(f_str, &pos);
          if (pos != f_str.size()) {
            LOG_ERROR("Float conversion failed");
            return std::unexpected{LexerError::value_conversion_failed};
          }
          return val;
        } catch (std::invalid_argument& err) {
          LOG_ERROR("Float conversion failed");
          return std::unexpected{LexerError::value_conversion_failed};
        }
      });
}

template <>
auto ConfigLexer::consume<size_t>() noexcept -> LexerResult<size_t> {
  return consume<std::string>().and_then(
      [](const std::string& d_str) -> LexerResult<size_t> {
        try {
          size_t pos = 0;
          size_t val = std::stoull(d_str, &pos);
          if (pos != d_str.size()) {
            LOG_ERROR("size_t conversion failed");
            return std::unexpected{LexerError::value_conversion_failed};
          }
          return val;
        } catch (std::invalid_argument& err) {
          LOG_ERROR("size_t conversion failed");
          return std::unexpected{LexerError::value_conversion_failed};
        }
      });
}

template <>
auto ConfigLexer::consume<bool>() noexcept -> LexerResult<bool> {
  auto b_str = consumeValueAsLowerString();
  if (!b_str) { return std::unexpected{b_str.error()}; }

  if (b_str == "true" || b_str == "1") { return true; }
  if (b_str == "false" || b_str == "0") { return false; }

  LOG_ERROR("Bool conversion failed");
  return std::unexpected{LexerError::value_conversion_failed};
}
};
