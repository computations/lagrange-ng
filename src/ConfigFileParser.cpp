#include "ConfigFileParser.hpp"

namespace lagrange {
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

auto parse(ConfigLexer& lexer, ToLowerOption l) -> ParsingResult<std::string> {
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

auto parse_assignment(ConfigLexer& lexer, ToLowerOption l)
    -> ParsingResult<std::string> {
  if (auto r = lexer.expect(ConfigLexemeType::EQUALS_SIGN); !r) {
    return std::unexpected{ParsingError::expected_equals_sign};
  }

  return parse(lexer, l);
};

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
}  // namespace lagrange
