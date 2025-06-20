#include "ConfigFile.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <unordered_set>

#include "Fossil.hpp"
#include "logger.hpp"

namespace lagrange {

auto ConfigFile::parse_config_file(std::istream& instream) -> ConfigFile {
  ConfigFile config;
  std::string line;
  size_t line_number = 1;
  while (getline(instream, line)) {
    ConfigFileParser parser(line, line_number);
    if (auto r = parser.parse_line(config); !r) {
      LOG_ERROR("Failed to parse config file at {}", parser.describePosition())
      throw ConfigFileParsingError{"Failed to parse config file, giving up"};
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
