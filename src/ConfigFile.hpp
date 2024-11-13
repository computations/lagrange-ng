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

enum class OutputType { JSON, CSV, UNKNOWN };

enum class OptimizationMethod {
  NELDER_MEAD,
  BOBYQA,
  COBYLA,
  BFGS,
  DIRECT,
  STOGO,
  UNKNOWN,
};

class ConfigFile;

ConfigFile parse_config_file(std::istream& instream);

class ConfigFile {
 public:
  ConfigFile() = default;

  ConfigFile(std::istream& instream) : ConfigFile{parse_config_file(instream)} {
    finalize();
  }

  ConfigFile(const ConfigFile&) = default;
  ConfigFile(ConfigFile&&) = default;

  ConfigFile& operator=(const ConfigFile&) = default;
  ConfigFile& operator=(ConfigFile&&) = default;

  const std::filesystem::path& tree_filename() const;
  void tree_filename(const std::filesystem::path&);

  const std::filesystem::path& data_filename() const;
  void data_filename(const std::filesystem::path&);

  const std::filesystem::path& log_filename() const;
  void log_filename(const std::filesystem::path&);

  const std::filesystem::path& prefix() const;
  void prefix(const std::filesystem::path&);

  OutputType output_file_type() const;
  void output_file_type(OutputType);

  bool has_rate_matrix_filename() const;
  const std::filesystem::path& rate_matrix_filename() const;
  void rate_matrix_filename(const std::filesystem::path&);

  bool has_max_areas() const;
  RangeSize max_areas() const;
  void max_areas(RangeSize);

  const std::vector<std::string>& area_names() const;
  void area_names(const std::vector<std::string>&);

  const Periods& periods() const;
  void periods(const Periods&);

  void finalize_periods() {
    finalize_periods(_region_count.get(), _max_areas.get());
  }

  void finalize_periods(size_t regions, size_t max_areas) {
    _periods.setRegionCount(regions);
    _periods.setMaxAreas(max_areas);
  }

  const std::shared_ptr<MRCAEntry>& mrca(const MRCALabel&) const;
  const MRCAMap& mrcas() const;
  void add_mrca(const MRCALabel&, const std::shared_ptr<MRCAEntry>&);

  const std::vector<Fossil>& fossils() const;
  void add_fossil(const Fossil& f);

  bool marginal() const;
  void marginal(bool);

  bool compute_all_splits() const;
  void compute_all_splits(bool);

  bool compute_all_states() const;
  void compute_all_states(bool);

  const std::unordered_set<MRCALabel>& state_nodes() const;
  void state_nodes(const std::unordered_set<MRCALabel>&);

  const std::unordered_set<MRCALabel>& split_nodes() const;
  void split_nodes(const std::unordered_set<MRCALabel>&);

  PeriodParams period_params() const;
  void period_params(double, double);

  double lh_epsilon() const;
  void lh_epsilon(double);

  bool has_expm_mode() const;
  const LagrangeEXPMComputationMode& expm_mode() const;
  void expm_mode(const LagrangeEXPMComputationMode&);

  AlignmentFileType alignment_file_type() const;
  void alignment_file_type(const AlignmentFileType&);

  RangeSize region_count() const;
  void region_count(RangeSize);

  size_t workers() const;
  void workers(size_t);

  size_t threads_per_worker() const;
  void threads_per_worker(size_t);

  const LagrangeOperationMode& run_mode() const;
  void run_mode(const LagrangeOperationMode&);

  OptimizationMethod opt_method() const;
  void opt_method(const OptimizationMethod&);

  std::filesystem::path jsonResultsFilename() const;
  std::filesystem::path nodeTreeFilename() const;
  std::filesystem::path scaledTreeFilename() const;
  std::filesystem::path splitsTreeFilename() const;
  std::filesystem::path statesTreeFilename() const;
  std::filesystem::path splitsCSVResultsFilename() const;
  std::filesystem::path statesCSVResultsFilename() const;
  std::filesystem::path periodsCSVResultsFilename() const;
  std::filesystem::path distributionsCSVResultsFilename() const;
  std::filesystem::path nodeInfoCSVResultsFilename() const;

  bool computeStates() const;
  bool computeSplits() const;
  bool computeStatesStrict() const;

 private:
  static ConfigFile parse_config_file(std::istream& instream);

  bool validate_and_make_prefix() {
    bool ok = true;
    try {
      if (!_prefix.parent_path().empty()) {
        std::filesystem::create_directories(_prefix.parent_path());
      }
    } catch (const std::filesystem::filesystem_error& err) {
      MESSAGE(ERROR, err.what());
      ok = false;
    }

    return ok;
  }

  void check_prefix() {
    if (prefix().empty()) { prefix(tree_filename()); }
    if (prefix().filename().string()[0] == '.') {
      MESSAGE(
          WARNING,
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

  void setup_log() {
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

  Option<RangeSize> _region_count;
  Option<RangeSize> _max_areas;
  Option<size_t> _workers;
  Option<size_t> _threads_per_worker;
  Option<LagrangeOperationMode> _run_mode{LagrangeOperationMode::OPTIMIZE};
  Option<OptimizationMethod> _opt_method{OptimizationMethod::BOBYQA};
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
