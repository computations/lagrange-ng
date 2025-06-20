#include "ConfigFileParser.hpp"

#include "ConfigFile.hpp"
#include "Periods.hpp"

namespace lagrange {

auto ConfigFileParser::parse(ToLowerOption l) -> ParsingResult<std::string> {
  if (auto r = _lexer.expectAndConsume(ConfigLexemeType::VALUE, l)) {
    return *r;
  } else if (r == std::unexpected{LexerError::value_conversion_failed}) {
    return std::unexpected{ParsingError::conversion_error};
  } else if (r == std::unexpected{LexerError::expected_wrong_token}) {
    return std::unexpected{ParsingError::expected_value};
  } else {
    return std::unexpected{ParsingError::lexing_error};
  }
}

auto ConfigFileParser::determine_fossil_type(
    const std::string& fossil_type_string) -> FossilType {
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

auto ConfigFileParser::parse_fossil() -> ParsingResult<Fossil> {
  if (!_lexer.expect(ConfigLexemeType::VALUE)) {
    return std::unexpected{ParsingError::expected_value};
  }
  auto fossil_type_string = *_lexer.consumeValueAsLowerString();

  FossilType ft = determine_fossil_type(fossil_type_string);

  auto mrca_name = parse<std::string>();
  if (!mrca_name) { return std::unexpected{mrca_name.error()}; }

  if (!_lexer.expect(ConfigLexemeType::EQUALS_SIGN)) {
    return std::unexpected{ParsingError::expected_equals_sign};
  }

  if (!_lexer.expect(ConfigLexemeType::VALUE)) {
    return std::unexpected{ParsingError::expected_value};
  }
  auto fossil_area = *_lexer.consumeAsRange();

  return Fossil{.mrca_name = *mrca_name,
                .clade = {},
                .age = ft == FossilType::BRANCH ? *parse<double>() : 0.0,
                .area = fossil_area,
                .type = ft};
}

auto ConfigFileParser::determine_output_file_type(
    const std::string& type_string) -> OutputType {
  if (type_string == "csv") { return OutputType::CSV; }
  if (type_string == "json") { return OutputType::JSON; }
  return OutputType::UNKNOWN;
}

auto ConfigFileParser::parse_opt_method(const std::string& method_string)
    -> OptimizationMethod {
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

auto ConfigFileParser::parse_assignment(ToLowerOption l)
    -> ParsingResult<std::string> {
  if (auto r = _lexer.expect(ConfigLexemeType::EQUALS_SIGN); !r) {
    return std::unexpected{ParsingError::expected_equals_sign};
  }

  return parse(l);
};

ParsingResult<void> ConfigFileParser::parse_period_include_statement(
    PeriodMap& period_map) {
  auto res = parse_and_check_id_map<std::vector<std::string>>(period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).include_areas = list;
  return {};
}

ParsingResult<void> ConfigFileParser::parse_period_exclude_statement(
    PeriodMap& period_map) {
  auto res = parse_and_check_id_map<std::vector<std::string>>(period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).exclude_areas = list;
  return {};
}

ParsingResult<void> ConfigFileParser::parse_period_start_statement(
    PeriodMap& period_map) {
  auto res = parse_and_check_id_map<double>(period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).start = list;
  return {};
}

ParsingResult<void> ConfigFileParser::parse_period_end_statement(
    PeriodMap& period_map) {
  auto res = parse_and_check_id_map<double>(period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).end = list;
  return {};
}

ParsingResult<void> ConfigFileParser::parse_period_matrix_statement(
    PeriodMap& period_map) {
  auto res = parse_and_check_id_map<std::filesystem::path>(period_map);

  if (!res) { return std::unexpected{res.error()}; }

  auto [period_name, list] = *res;

  period_map.at(period_name).adjustment_matrix_filename = list;
  return {};
}

ParsingResult<void> ConfigFileParser::parse_period_statement(
    PeriodMap& period_map) {
  auto value = parse<std::string>();
  if (!value) { return std::unexpected{value.error()}; }

  if (value == "include") {
    return parse_period_include_statement(period_map);
  } else if (value == "exclude") {
    return parse_period_exclude_statement(period_map);
  } else if (value == "start") {
    return parse_period_start_statement(period_map);
  } else if (value == "end") {
    return parse_period_end_statement(period_map);
  } else if (value == "matrix") {
    return parse_period_matrix_statement(period_map);
  } else {
    if (period_map.contains(*value)) {
      LOG_ERROR("Duplicate period declaration for period id '{}' at {}",
                *value,
                _lexer.describePosition());
      return std::unexpected{ParsingError::duplicate_id};
    } else {
      period_map.insert({*value, {}});
      return {};
    }
  }

  if (!_lexer.expect(ConfigLexemeType::END)) {
    return std::unexpected{ParsingError::expected_end};
  }
  return {};
}

auto ConfigFileParser::parse_line(ConfigFile& config) -> ParsingResult<void> {
  if (_lexer.peak() == ConfigLexemeType::END) { return {}; }

  auto config_value = parse<std::string>();

  if (config_value == "treefile") {
    if (auto r = parse_and_assign(config._tree_filename); !r) { return r; }
  } else if (config_value == "datafile") {
    if (auto r = parse_and_assign(config._data_filename); !r) { return r; }
  } else if (config_value == "areanames") {
    if (auto r = parse_and_assign(config._area_names); !r) { return r; }
  } else if (config_value == "prefix") {
    if (auto r = parse_and_assign(config._prefix); !r) { return r; }
  } else if (config_value == "periods") {
    std::vector<double> period_times;
    if (auto r = parse_and_assign(period_times)) {
      config._periods = Periods(period_times);
    } else {
      return r;
    }
  } else if (config_value == "mrca") {
    if (auto r = parse_id_assignment<std::vector<std::string>>()) {
      auto [mrca_name, mrca_entries] = *r;
      config._mrcas[mrca_name] = std::make_shared<MRCAEntry>(
          MRCAEntry{.clade = std::move(mrca_entries)});
    } else {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "fossil") {
    if (auto r = parse_fossil()) {
      config._fossils.push_back(*r);
    } else {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "splits") {
    if (_lexer.peak() == ConfigLexemeType::VALUE) {
      if (auto r = parse<std::unordered_set<MRCALabel>>()) {
        config._split_nodes = *r;
      } else {
        return std::unexpected{r.error()};
      }
    } else {
      config._all_splits = true;
    }
  }

  else if (config_value == "states") {
    if (_lexer.peak() == ConfigLexemeType::VALUE) {
      if (auto r = parse<std::unordered_set<MRCALabel>>()) {
        config._state_nodes = *r;
      } else {
        return std::unexpected{r.error()};
      }
    } else {
      config._all_states = true;
    }
  } else if (config_value == "dispersion") {
    if (auto r = parse_and_assign(config._dispersion); !r) { return r; }
  } else if (config_value == "extinction") {
    if (auto r = parse_and_assign(config._extinction); !r) { return r; }
  } else if (config_value == "lh-epsilon") {
    if (auto r = parse_and_assign(config._lh_epsilon); !r) { return r; }
  } else if (config_value == "workers") {
    if (auto r = parse_and_assign(config._workers); !r) { return r; }
  } else if (config_value == "threads-per-worker") {
    if (auto r = parse_and_assign(config._threads_per_worker); !r) { return r; }
  } else if (config_value == "maxareas") {
    if (auto r = parse_and_assign(config._max_areas); !r) { return r; }
  } else if (config_value == "expm-mode") {
    auto expm_type = parse_assignment<std::string>();

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
    auto mode_type = parse_assignment<std::string>();
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
    auto align_type_str = parse_assignment<std::string>();
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
    if (auto r = parse_and_assign(config._log_filename); !r) {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "output-type") {
    if (auto r = parse_assignment(ToLowerOption::lower)) {
      config.output_file_type(determine_output_file_type(*r));
    } else {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "opt-method") {
    if (auto r = parse_assignment(ToLowerOption::lower)) {
      config.opt_method(parse_opt_method(*r));
    } else {
      return std::unexpected{r.error()};
    }
  } else if (config_value == "allow-ambiguous") {
    if (_lexer.peak() == ConfigLexemeType::EQUALS_SIGN) {
      if (auto r = parse_assignment<bool>()) {
        config.allow_ambigious(*r);
      } else {
        return std::unexpected{r.error()};
      }
    } else {
      config.allow_ambigious(true);
    }
  } else if (config_value == "dump-graph") {
    if (_lexer.peak() == ConfigLexemeType::EQUALS_SIGN) {
      if (auto r = parse_assignment<bool>()) {
        config.dump_graph(*r);
      } else {
        return std::unexpected{r.error()};
      }
    } else {
      config.dump_graph(true);
    }
  } else if (config_value == "lwr-threshold") {
    if (auto r = parse_and_assign(config._lh_epsilon); !r) {
      return std::unexpected{r.error()};
    }
  } else if (config_value) {
    LOG_ERROR("Option '{}' at {} was not recognized",
              *config_value,
              _lexer.describePosition());
    return std::unexpected{ParsingError::invalid_option};
  } else {
    LOG_ERROR("Failed to parse the option at {}", _lexer.describePosition());
    return std::unexpected{ParsingError::unknown_error};
  }

  if (auto r = _lexer.expect(ConfigLexemeType::END); !r) {
    return std::unexpected{ParsingError::expected_end};
  }

  return {};
}

auto ConfigFileParser::describePosition() const -> std::string {
  return _lexer.describePosition();
}

}  // namespace lagrange
