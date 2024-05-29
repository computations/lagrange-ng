#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <filesystem>
#include <logger.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#include "Alignment.hpp"
#include "Common.hpp"
#include "Fossil.hpp"
#include "MRCA.hpp"
#include "Periods.hpp"
#include "Utils.hpp"

namespace lagrange {
constexpr size_t KRYLOV_RANGE_COUNT_THRESHOLD = 5;

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

  bool has_rate_matrix_filename() const;
  const std::filesystem::path& rate_matrix_filename() const;
  void rate_matrix_filename(const std::filesystem::path&);

  bool has_max_areas() const;
  size_t max_areas() const;
  void max_areas(size_t);

  const std::vector<std::string>& area_names() const;
  void area_names(const std::vector<std::string>&);

  const Periods& periods() const;
  void periods(const Periods&);

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

  const std::vector<MRCALabel>& state_nodes() const;
  void state_nodes(const std::vector<MRCALabel>&);

  const std::vector<MRCALabel>& split_nodes() const;
  void split_nodes(const std::vector<MRCALabel>&);

  PeriodParams period_params() const;
  void period_params(double, double);

  double lh_epsilon() const;
  void lh_epsilon(double);

  bool has_expm_mode() const;
  const LagrangeEXPMComputationMode& expm_mode() const;
  void expm_mode(const LagrangeEXPMComputationMode&);

  AlignmentFileType alignment_file_type() const;
  void alignment_file_type(const AlignmentFileType&);

  size_t region_count() const;
  void region_count(size_t);

  size_t workers() const;
  void workers(size_t);

  size_t threads_per_worker() const;
  void threads_per_worker(size_t);

  const LagrangeOperationMode& run_mode() const;
  void run_mode(const LagrangeOperationMode&);

  std::filesystem::path resultsFilename() const;
  std::filesystem::path NodeTreeFilename() const;
  std::filesystem::path scaledTreeFilename() const;

  bool computeStates() const;
  bool computeSplits() const;

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

  Option<std::filesystem::path> _rate_matrix_filename;

  Option<size_t> _max_areas;

  std::vector<std::string> _area_names;

  Periods _periods;
  MRCAMap _mrcas;

  std::vector<Fossil> _fossils;

  bool _marginal = true;
  bool _all_splits = false;
  bool _all_states = false;

  std::vector<MRCALabel> _state_nodes;
  std::vector<MRCALabel> _split_nodes;

  double _dispersion = 0.01;
  double _extinction = 0.01;

  double _lh_epsilon = 1e-9;

  Option<LagrangeEXPMComputationMode> _expm_mode;

  Option<AlignmentFileType> _alignment_file_type;

  Option<size_t> _region_count;
  Option<size_t> _workers;
  Option<size_t> _threads_per_worker;
  Option<LagrangeOperationMode> _run_mode{LagrangeOperationMode::OPTIMIZE};
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
