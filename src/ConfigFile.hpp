#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <filesystem>
#include <logger.hpp>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "Alignment.hpp"
#include "Common.hpp"
#include "Fossil.hpp"
#include "MRCA.hpp"
#include "Periods.hpp"
#include "Utils.hpp"

namespace lagrange {
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

class ConfigFile;

auto parse_config_file(std::istream& instream) -> ConfigFile;

class ConfigFile {
 public:
  ConfigFile() = default;

  ConfigFile(std::istream& instream) : ConfigFile{parse_config_file(instream)} {
    finalize();
  }

  ConfigFile(const ConfigFile&) = default;
  ConfigFile(ConfigFile&&) = default;

  auto operator=(const ConfigFile&) -> ConfigFile& = default;
  auto operator=(ConfigFile&&) -> ConfigFile& = default;

  auto tree_filename() const -> const std::filesystem::path&;
  void tree_filename(const std::filesystem::path&);

  auto data_filename() const -> const std::filesystem::path&;
  void data_filename(const std::filesystem::path&);

  auto log_filename() const -> const std::filesystem::path&;
  void log_filename(const std::filesystem::path&);

  auto prefix() const -> const std::filesystem::path&;
  void prefix(const std::filesystem::path&);

  auto output_file_type() const -> OutputType;
  void output_file_type(OutputType);

  auto has_rate_matrix_filename() const -> bool;
  auto rate_matrix_filename() const -> const std::filesystem::path&;
  void rate_matrix_filename(const std::filesystem::path&);

  auto has_max_areas() const -> bool;
  auto max_areas() const -> size_t;
  void max_areas(size_t);

  auto area_names() const -> const std::vector<std::string>&;
  void area_names(const std::vector<std::string>&);

  auto periods() const -> const Periods&;
  void periods(const Periods&);

  void finalize_periods() {
    finalize_periods(_region_count.get(), _max_areas.get());
  }

  void finalize_periods(size_t regions, size_t max_areas) {
    _periods.setRegionCount(regions);
    _periods.setMaxAreas(max_areas);
  }

  auto mrca(const MRCALabel&) const -> const std::shared_ptr<MRCAEntry>&;
  auto mrcas() const -> const MRCAMap&;
  void add_mrca(const MRCALabel&, const std::shared_ptr<MRCAEntry>&);

  auto fossils() const -> const std::vector<Fossil>&;
  void add_fossil(const Fossil& fossil);

  auto marginal() const -> bool;
  void marginal(bool);

  auto compute_all_splits() const -> bool;
  void compute_all_splits(bool);

  auto compute_all_states() const -> bool;
  void compute_all_states(bool);

  auto state_nodes() const -> const std::unordered_set<MRCALabel>&;
  void state_nodes(const std::unordered_set<MRCALabel>&);

  auto split_nodes() const -> const std::unordered_set<MRCALabel>&;
  void split_nodes(const std::unordered_set<MRCALabel>&);

  auto period_params() const -> PeriodParams;
  void period_params(double, double);

  auto lh_epsilon() const -> double;
  void lh_epsilon(double);

  auto has_expm_mode() const -> bool;
  auto expm_mode() const -> const LagrangeEXPMComputationMode&;
  void expm_mode(const LagrangeEXPMComputationMode&);

  auto alignment_file_type() const -> AlignmentFileType;
  void alignment_file_type(const AlignmentFileType&);

  auto region_count() const -> size_t;
  void region_count(size_t);

  auto workers() const -> size_t;
  void workers(size_t);

  auto threads_per_worker() const -> size_t;
  void threads_per_worker(size_t);

  auto run_mode() const -> const LagrangeOperationMode&;
  void run_mode(const LagrangeOperationMode&);

  auto opt_method() const -> OptimizationMethod;
  void opt_method(const OptimizationMethod&);

  auto allow_ambigious() const -> bool;
  void allow_ambigious(bool);

  auto dump_graph() const -> bool;
  void dump_graph(bool);

  auto jsonResultsFilename() const -> std::filesystem::path;
  auto nodeTreeFilename() const -> std::filesystem::path;
  auto scaledTreeFilename() const -> std::filesystem::path;
  auto splitsTreeFilename() const -> std::filesystem::path;
  auto statesTreeFilename() const -> std::filesystem::path;
  auto splitsCSVResultsFilename() const -> std::filesystem::path;
  auto statesCSVResultsFilename() const -> std::filesystem::path;
  auto periodsCSVResultsFilename() const -> std::filesystem::path;
  auto distributionsCSVResultsFilename() const -> std::filesystem::path;
  auto nodeInfoCSVResultsFilename() const -> std::filesystem::path;

  auto forwardGraphFilename() const -> std::filesystem::path;
  auto reverseGraphFilename() const -> std::filesystem::path;

  auto computeStates() const -> bool;
  auto computeSplits() const -> bool;
  auto computeStatesStrict() const -> bool;

 private:
  static auto parse_config_file(std::istream& instream) -> ConfigFile;

  auto validate_and_make_prefix() -> bool {
    bool OllKorrect = true;
    try {
      if (!_prefix.parent_path().empty()) {
        std::filesystem::create_directories(_prefix.parent_path());
      }
    } catch (const std::filesystem::filesystem_error& err) {
      LOG(ERROR, "Failed to create prefix directory: {}", err.what());
      OllKorrect = false;
    }

    return OllKorrect;
  }

  void check_prefix() {
    if (prefix().empty()) { prefix(tree_filename()); }
    if (prefix().filename().string()[0] == '.') {
      LOG(WARNING,
          "The current prefix starts with a dot. The results will be hidden "
          "on unix-like systems");
    }
  }

  void set_mrcas_for_fossils() {
    for (auto& f : _fossils) { f.clade = mrca(f.mrca_name); }
  }

  void set_threads() {
    if (!_workers.hasValue()) { workers(1); }
    if (!_threads_per_worker.hasValue()) { _threads_per_worker = 1; }
  }

  void setup_log() const {
    if (!log_filename().empty()) {
      logger::get_log_states().add_file_stream(
          log_filename().c_str(),
          INFO | IMPORTANT | WARNING | ERROR | PROGRESS | DEBUG);
    }
  }

  void finalize() {
    setup_log();
    set_threads();
    set_mrcas_for_fossils();
    check_prefix();
    validate_and_make_prefix();
  }

  std::filesystem::path _tree_filename;
  std::filesystem::path _data_filename;
  std::filesystem::path _log_filename;
  std::filesystem::path _prefix;

  Option<OutputType> _output_file_type;

  Option<std::filesystem::path> _rate_matrix_filename;

  Option<size_t> _max_areas;

  std::vector<std::string> _area_names;

  Periods _periods;
  MRCAMap _mrcas;

  std::vector<Fossil> _fossils;

  bool _marginal = true;
  bool _all_splits = false;
  bool _all_states = false;

  std::unordered_set<MRCALabel> _state_nodes;
  std::unordered_set<MRCALabel> _split_nodes;

  double _dispersion = 0.01;
  double _extinction = 0.01;

  double _lh_epsilon = 1e-9;

  Option<LagrangeEXPMComputationMode> _expm_mode;

  Option<AlignmentFileType> _alignment_file_type;

  Option<size_t> _region_count;
  Option<size_t> _workers;
  Option<size_t> _threads_per_worker;
  Option<LagrangeOperationMode> _run_mode{LagrangeOperationMode::OPTIMIZE};
  Option<OptimizationMethod> _opt_method{OptimizationMethod::BOBYQA};

  Option<bool> _allow_ambigious;
  Option<bool> _dump_graph;
};

class ConfigFileLexingError : public std::runtime_error {
 public:
  ConfigFileLexingError(const std::string& msg) : std::runtime_error{msg} {}
};

class ConfigFileParsingError : public std::runtime_error {
 public:
  ConfigFileParsingError(const std::string& msg) : std::runtime_error{msg} {}
};

}  // namespace lagrange
#endif
