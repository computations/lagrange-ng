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
    if (auto r = parse_and_assign(config._log_filename, lexer); !r) {
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
    return std::unexpected{ParsingError::invalid_option};
  } else {
    LOG_ERROR("Failed to parse the option on line {}", line_number);
    return std::unexpected{ParsingError::unknown_error};
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
        LOG_ERROR("Failed to parse config on line {}, {}",
                  line_number,
                  lexer.describePosition());
        throw std::runtime_error{"failed to parse the config file"};
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
