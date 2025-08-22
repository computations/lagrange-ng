#include "ConfigFileParser.hpp"

#include <cctype>
#include <functional>
#include <string_view>

#ifdef __cpp_lib_generator
  #include <generator>
#endif

#include "ConfigFile.hpp"
#include "logger.hpp"

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

auto ConfigFileParser::determine_opt_method(const std::string& method_string)
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
    LOG_ERROR("Expected an equal sign at {}", _lexer.describePosition());
    return std::unexpected{ParsingError::expected_equals_sign};
  }

  return parse(l);
};

ParsingResult<void> ConfigFileParser::parse_period_sub_statement(
    PeriodConfig& period) {
  auto command = parse<std::string>();

  if (!command) {
    LOG_ERROR(
        "Expected one of [start | end | matrix | exclude | include] at {}",
        _lexer.describePosition());
    return std::unexpected{ParsingError::invalid_option};
  }

  if (command == "start") { return parse_and_assign(period.start); }
  if (command == "end") { return parse_and_assign(period.end); }
  if (command == "matrix") {
    return parse_and_assign(period.adjustment_matrix_filename);
  }
  if (command == "include") { return parse_and_assign(period.include_areas); }
  if (command == "exclude") { return parse_and_assign(period.exclude_areas); }
  return std::unexpected{ParsingError::invalid_option};
}

ParsingResult<void> ConfigFileParser::parse_period_statement(
    PeriodConfigMap& period_map) {
  auto value = parse<std::string>();
  if (!value) {
    LOG_ERROR("Expected a value after a period statement at {}",
              _lexer.describePosition());
    return std::unexpected{value.error()};
  }

  if (_lexer.peak() == ConfigLexemeType::END) {
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

  if (!period_map.contains(*value)) {
    LOG_ERROR("Period {} was not declared before use at {}",
              *value,
              _lexer.describePosition());
    return std::unexpected{ParsingError::invalid_option};
  }

  if (auto r = parse_period_sub_statement(period_map.at(*value)); !r) {
    return r;
  }

  if (!_lexer.expect(ConfigLexemeType::END)) {
    LOG_ERROR(
        "Expected end token at {}. There probably are some extra characters "
        "on this line.",
        _lexer.describePosition());
    return std::unexpected{ParsingError::expected_end};
  }
  return {};
}

ConfigFileParser::ActionMapType ConfigFileParser::_config_action_map{
    {
        "treefile",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._tree_filename);
            }},
            ._name{"Tree file"},
            ._help{"File containing a rooted phlyogenetic tree in newick "
                   "format."},
            ._usage{"treefile = <FILENAME>"},
            ._examples{"treefile = example.nwk"},
        },
    },
    {
        "datafile",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._data_filename);
            }},
            ._name{"Range data file"},
            ._help{"File containing the tip and associated ranges as "
                   "either a phylip for fasta file. Ranges should be specified "
                   "as binary strings (e.g. 01011). "
                   "Areas are specified in the order of area names below. "
                   "For example, if the area names are E, W, and N, then the "
                   "string 011 corrisponds to a range where the species is "
                   "present in areas W and N, but not in E."},
            ._usage{"datafile = <FILENAME>"},
            ._examples{"datafile = example.phy"},
        },
    },
    {
        "areanames",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._area_names);
            }},
            ._name{"Ares name list"},
            ._help{
                "List of the area names, which will be used for output. "
                "There "
                "must be as many area names as there are areas in the data. "
                "Area names can contain spaces, so long as the whole name is "
                "quoted (e.g. 'North America')"},
            ._usage{"areanames = <NAME1> [<NAME2> ...]"},
            ._examples{
                "areanames = West South North",
                "areanames = 'North America' 'South America'",
                "areanames = \"West Europe\" \"South Europe\"",
            },
        },
    },
    {
        "prefix",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._prefix);
            }},
            ._name{"Output file prefix"},
            ._help{"Specifies the output file prefix. This is a path that "
                   "will be "
                   "prepended to the output files. The prefix can included "
                   "directories, and if they do not exist, they will be "
                   "created"},
            ._usage{"prefix = <PATH>"},
            ._examples{
                "prefix = results/test.1 (makes results/test.1.results.json)",
            },
        },
    },
    {
        "period",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_period_statement(config._period_map);
            }},
            ._name{"Create and manipulate periods."},
            ._help{"Create and manipulate periods. Each period is declared by "
                   "giving it a name (see the first example). From there, the "
                   "period can be modified via a series of commands. None of "
                   "the options are strictly required, with some caveats. In "
                   "the case that there is no start time, the start time is "
                   "assumed to be 0.0, and if there is no end time, then it is "
                   "assumed to be infinity. Periods cannot be overlapping, so "
                   "if more than one period does not have start or end times, "
                   "then an error will be thrown. For more information about "
                   "periods, please see README.md"},
            ._usage{"period <PERIOD> [(start|end|matrix|include|exclude) = "
                    "<VALUE>]"},
            ._examples{
                "period foo",
                "period foo start = 1.0",
                "period foo end = 2.0",
                "period foo include = 101",
                "period foo exclude = 011",
                "period foo matrix = matrix_file.csv",
                "period foo matrix = 'matrix file with spaces.csv'",
            },
        },
    },
    {
        "mrca",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              if (auto r = p.parse_id_assignment<std::vector<std::string>>()) {
                auto [mrca_name, mrca_entries] = *r;
                config._mrcas[mrca_name] = std::make_shared<MRCAEntry>(
                    MRCAEntry{.clade = std::move(mrca_entries)});
                return {};
              } else {
                return std::unexpected{r.error()};
              }
            }},
            ._name{"Inner node specification (MRCA)"},
            ._help{"Assign a label to an inner node for later use in other "
                   "options. The node is specified by listing taxa, and label "
                   "will corrispond to the node which is the MRCA of all taxa "
                   "in the list. The entire clade does not need to be "
                   "specified, only the \"bracketing\" taxa for the clade."},
            ._usage{"mrca <MRCA LABEL> = <TAXON> [<TAXON> ...]"},
            ._examples{"mrca crown = Penguin Pidgeon (crown is the label for "
                       "the MRCA)"},
        },
    },
    {
        "fossil",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              auto r = p.parse_fossil();
              if (!r) { return std::unexpected{r.error()}; }
              config._fossils.push_back(*r);
              return {};
            }},
            ._name{"Fossil specification"},
            ._help{"Specify a fossil constraint for an ancestral node. Can "
                   "either be an include, exclude, or fixed constraint. "
                   "Nodes constrained with an include fossil will only "
                   "consider ancestral states with areas that are set in the "
                   "given range, and vice versa for an exclude constraint. "
                   "Fixed "
                   "constraints require that the inner node's ancestral "
                   "state is "
                   "the given range."},
            ._usage{"fossil [include|exclude|fixed] <MRCA LABEL> = <RANGE>"},
            ._examples{
                "fossil include crown = 0010 (will include the 3rd area)",
                "fossil exclude crown = 0010 (will exclude the 3rd area)",
                "fossil fixed crown = 0010 (will be 0010)",
            },
        },
    },
    {
        "splits",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              if (p._lexer.peak() == ConfigLexemeType::VALUE) {
                if (auto r = p.parse<std::unordered_set<MRCALabel>>()) {
                  config._split_nodes = *r;
                } else {
                  return std::unexpected{r.error()};
                }
              } else {
                config._all_splits = true;
              }
              return {};
            }},
            ._name{"Split result specification"},
            ._help{
                "Specify which ancestral splits to compute. Can be a list "
                "of nodes, each specified by the mrca command. If no "
                "nodes are provided, then ancestral splits for all nodes are "
                "computed."},
            ._usage{"splits [<MRCA> ...]"},
            ._examples{
                "splits (compute splits for all nodes)",
                "splits crown stem (compute splits for nodes crown and stem)",
            },
        },
    },
    {
        "states",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              if (p._lexer.peak() == ConfigLexemeType::VALUE) {
                if (auto r = p.parse<std::unordered_set<MRCALabel>>()) {
                  config._state_nodes = *r;
                } else {
                  return std::unexpected{r.error()};
                }
              } else {
                config._all_states = true;
              }
              return {};
            }},
            ._name{"State result specification"},
            ._help{"Specify which ancestral states to compute. Can be a list "
                   "of nodes, each specified by the mrca command. If no "
                   "nodes are "
                   "provided, then ancestral states for all nodes are "
                   "computed."},
            ._usage{"splits [<MRCA> ...]"},
            ._examples{
                "states (compute states for all nodes)",
                "states crown stem (compute states for nodes crown and stem)",
            },
        },
    },
    {
        "dispersion",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._dispersion);
            }},
            ._name{"Dispersion rate"},
            ._help{"Fix the dispersion rate in evaluation mode. Has no "
                   "effect in "
                   "optimize mode."},
            ._usage{"dispersion = <FLOAT>"},
            ._examples{"dispersion = 0.1"},
        },
    },
    {
        "extinction",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._extinction);
            }},
            ._name{"Extinction rate"},
            ._help{"Fix the extinction rate in evaluation mode. Has no "
                   "effect in "
                   "optimize mode."},
            ._usage{"extinction = <FLOAT>"},
            ._examples{"extinction = 0.1"},
        },
    },
    {
        "lh-epsilon",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._lh_epsilon);
            }},
            ._name{"Log-likelihood threshold"},
            ._help{"Set the log-likelihood epsilon, which is the threshold "
                   "before aborting the optimization routine."},
            ._usage{"lh-epsilon = <FLOAT>"},
            ._examples{"lh-epsilon = 0.001"},
        },
    },
    {
        "workers",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._workers);
            }},
            ._name{"Worker thread specification"},
            ._help{"Set the number of workers. A good number to use for the "
                   "number of workers is the number of physical cores on this "
                   "machine. Defaults to 1."},
            ._usage{"workers = <INTEGER>"},
            ._examples{"workers = 4"},
        },
    },
    {
        "threads-per-worker",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._threads_per_worker);
            }},
            ._name{"Threads per worker specification"},
            ._help{"Set the number of threads per workers. In general, this "
                   "should be left to the default unless you have a large "
                   "number of areas (8 and over). In that case, it can be "
                   "faster to halve the number of workers, and set "
                   "threads-per-worker to 2. Defaults to 1."},
            ._usage{"threads-per-worker = <INTEGER>"},
            ._examples{"threads-per-workers = 2"},
        },
    },
    {
        "maxareas",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._max_areas);
            }},
            ._name{"Maximum area limit"},
            ._help{"Set the maximum areas allowed in the analysis. Ranges "
                   "which occupy more than the number of max areas will "
                   "not be considered when computing the likelihood or when "
                   "computing ancestral states/splits."},
            ._usage{"maxareas = <INTEGER>"},
            ._examples{"maxareas = 3 (ranges with more than 3 areas are not "
                       "considred)"},
        },
    },
    {
        "expm-mode",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              auto expm_type = p.parse_assignment<std::string>();

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
              return {};
            }},
            ._name{"Exponential method specification"},
            ._help{"Sets the strategy to compute the matrix exponential. Avoid "
                   "setting unless you are having issues with numerical "
                   "stability. In that case, the pade method is more stable, "
                   "but much slower. Defaults to adaptive mode, which "
                   "automatically switches between krylov and pade when errors "
                   "are detected."},
            ._usage{"expm-mode = [krylov|pade|adaptive]"},
            ._examples{},
        },
    },
    {
        "mode",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              auto mode_type = p.parse_assignment<std::string>();
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
              return {};
            }},
            ._name{"Run mode specification"},
            ._help{"Sets the run mode for results. Setting this to optimize "
                   "will optimize all free parameters before computing "
                   "ancestral states.  On the other hand, evaluate mode will "
                   "skip "
                   "the optimization step, and only compute results. Defaults "
                   "to optimize."},
            ._usage{"mode = [optimize|evaluate]"},
            ._examples{},
        },
    },
    {
        "alignment-format",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              auto align_type_str = p.parse_assignment<std::string>();
              if (align_type_str == "fasta") {
                config._alignment_file_type = AlignmentFileType::FASTA;
              } else if (align_type_str == "phylip") {
                config._alignment_file_type = AlignmentFileType::PHYLIP;
              } else if (align_type_str) {
                LOG_ERROR("Unknown alignment type '{}'", *align_type_str);
              } else {
                return std::unexpected{align_type_str.error()};
              }
              return {};
            }},
            ._name{"Input data format"},
            ._help{"Sets the alignment format for the data file. Set to "
                   "phylip by default."},
            ._usage{"alignment-format = [phylip|fasta]"},
            ._examples{},
        },
    },
    {
        "logfile",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              return p.parse_and_assign(config._log_filename);
            }},
            ._name{"Log output filename"},
            ._help{"Sets the filename for the logfile which will contain the "
                   "messages printed to screen. By default, screen output is "
                   "not written to a file. There are no results contained in "
                   "the logfile that are not written to a results file. Only "
                   "specify this for debugging purposes."},
            ._usage{"logfile = <FILENAME>"},
            ._examples{"logfile = results/example.log"},
        },
    },
    {
        "output-type",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              if (auto r = p.parse_assignment(ToLowerOption::lower)) {
                config.output_file_type(determine_output_file_type(*r));
                return {};
              } else {
                return std::unexpected{r.error()};
              }
            }},
            ._name{"Output filetype"},
            ._help{"Sets the output format for the result files. By "
                   "default, results are saved as a single JSON file. However, "
                   "for very large datasets it will be faster to write the "
                   "results as CSV files. For a full description of what the "
                   "output files of Lagrange-NG look like, please read the "
                   "README.md"},
            ._usage{"output-type = [json|csv]"},
            ._examples{},
        },
    },
    {
        "opt-method",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              if (auto r = p.parse_assignment(ToLowerOption::lower)) {
                config.opt_method(p.determine_opt_method(*r));
              } else {
                return std::unexpected{r.error()};
              }
              return {};
            }},
            ._name{"Optimization Method"},
            ._help{"Sets the optimization algorithm which will be used to find "
                   "the optimal model parameters. In general, all methods are "
                   "fast enough for datasets with one period. However, if "
                   "there are several periods, it might be faster to set this "
                   "option to bfgs instead. Set to BOBYQA by default."},
            ._usage{"opt-method = "
                    "[nelder-mead|bobyqa|cobyla|bfgs|stogo|direct]"},
            ._examples{},
        },
    },
    {
        "allow-ambiguous",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              if (p._lexer.peak() == ConfigLexemeType::EQUALS_SIGN) {
                if (auto r = p.parse_assignment<bool>()) {
                  config.allow_ambigious(*r);
                } else {
                  return std::unexpected{r.error()};
                }
              } else {
                config.dump_graph(true);
              }
              return {};
            }},
            ._name{"Ambiguous characters"},
            ._help{"In the case where maxareas is given, but a tip state has a "
                   "number of occupied areas which exceeds that limit, then "
                   "this option will allow for those tip states to be "
                   "converted into ambiguous states. For example, if maxareas "
                   "= 2, and a tip has the range 111, the tip state will be "
                   "converted to the ambiguous state combining the ranges 011, "
                   "101, and 110. If this option is not set, then tips with "
                   "ranges that exceed the maxareas limit will cause the "
                   "program to quit with an error. If set, a warning will be "
                   "printed every time a state is converted. Set on by "
                   "default."},
            ._usage{"allow-ambiguous [= true|on|false|off]"},
            ._examples{},
        },
    },
    {
        "dump-graph",
        {
            ._action{[](ConfigFileParser& p,
                        ConfigFile& config) -> ParsingResult<void> {
              if (p._lexer.peak() == ConfigLexemeType::EQUALS_SIGN) {
                if (auto r = p.parse_assignment<bool>()) {
                  config.dump_graph(*r);
                } else {
                  return std::unexpected{r.error()};
                }
              } else {
                config.dump_graph(true);
              }
              return {};
            }},
            ._name{"Computation graph output"},
            ._help{"Produce a graphviz file specifying the operation graph for "
                   "both the forward and backward passes. Most users should "
                   "not need to set this option, as it just helps with "
                   "debugging the internals of Lagrange-NG's computation. Set "
                   "to off by default."},
            ._usage{"dump-graph [= true|on|false|off]"},
            ._examples{},
        },
    },
};

std::string_view break_line(std::string_view& str, size_t line_length) {
  size_t start_index = 0;
  while (start_index < str.length()) {
    if (start_index + line_length > str.length()) {
      auto ret = str.substr(start_index, line_length);
      str = {str.begin() + line_length, str.end()};
      return ret;
    }

    size_t cur_line_length = line_length;

    if (cur_line_length + start_index >= str.size()) {
      cur_line_length = str.size() - start_index - 1;
    }

    while (cur_line_length != 0
           && std::isspace(str[start_index + cur_line_length]) == 0) {
      cur_line_length -= 1;
    }
    auto ret = str.substr(start_index, cur_line_length);
    start_index += cur_line_length + 1;
    str = {str.begin() + start_index, str.end()};
    return ret;
  }
  return {};
}

#ifdef __cpp_lib_generator
std::generator<std::string_view> break_line(const std::string_view& str,
                                            size_t line_length) {
  size_t start_index = 0;
  while (start_index < str.length()) {
    if (start_index + line_length > str.length()) {
      co_yield str.substr(start_index, line_length);
      break;
    }

    size_t cur_line_length = line_length;

    if (cur_line_length + start_index >= str.size()) {
      cur_line_length = str.size() - start_index - 1;
    }

    while (cur_line_length != 0
           && std::isspace(str[start_index + cur_line_length]) == 0) {
      cur_line_length -= 1;
    }
    co_yield str.substr(start_index, cur_line_length);
    start_index += cur_line_length + 1;
  }
}
#endif

void ConfigFileParser::print_help_long() {
  MESSAGE_INFO(COLORIZE(ANSI_COLOR_BLUE, "Configuration file options:\n"));
  for (auto [option, config] : _config_action_map) {
    MESSAGE_INFO(COLORIZE(ANSI_COLOR_BLUE, "  {}"),
                 !config._usage.empty() ? config._usage : config._name);
    if (!config._help.empty()) {
#ifndef __cpp_lib_generator
      for (auto line : break_line(config._help, 78)) {
#else
      std::string_view str{config._help};
      while (true) {
        auto line = break_line(str, 78);
        if (line.empty()) { break; }
#endif
        MESSAGE_INFO("    {}", line);
      }
    }

    if (config._examples.size() != 0) {
      MESSAGE_INFO("\n    Examples:")
      for (auto e : config._examples) { MESSAGE_INFO("       {}", e); }
    }
    MESSAGE_INFO("");
  }
}

void ConfigFileParser::print_help_short() {
  MESSAGE_INFO("usage: lagrange-ng [help|<CONFIG FILE>]");
  MESSAGE_INFO(COLORIZE(ANSI_COLOR_BLUE, "Configuration file options:"));
  for (auto [option, config] : _config_action_map) {
    MESSAGE_INFO("  {}", !config._usage.empty() ? config._usage : config._name);
  }
}

auto ConfigFileParser::parse_line(ConfigFile& config) -> ParsingResult<void> {
  if (_lexer.peak() == ConfigLexemeType::END) { return {}; }

  auto parse_result = parse<std::string>();
  if (!parse_result) {
    LOG_ERROR("Failed to parse the option at {}", _lexer.describePosition());
    return std::unexpected{parse_result.error()};
  }

  auto config_value = *parse_result;

  if (_config_action_map.contains(config_value)) {
    auto r = _config_action_map[config_value]._action(*this, config);
    if (!r) {
      LOG_ERROR("Failed to parse option {} at {}",
                config_value,
                _lexer.describePosition());
      return r;
    }
  } else {
    LOG_ERROR("Option '{}' at {} is invalid\n",
              config_value,
              _lexer.describePosition());
    print_help_short();
    return std::unexpected{ParsingError::invalid_option};
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
