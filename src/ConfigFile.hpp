#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <logger.hpp>
#include <optional>
#include <ranges>
#include <string>
#include <unordered_set>
#include <vector>

#include "AdjustmentMatrix.hpp"
#include "Alignment.hpp"
#include "Common.hpp"
#include "ConfigFileParser.hpp"
#include "Fossil.hpp"
#include "MRCA.hpp"
#include "Periods.hpp"
#include "Utils.hpp"

namespace lagrange {

class ConfigFile;

class ConfigLexer;

// auto parse_config_file(std::istream& instream) -> ConfigFile;

class ConfigFile {
 public:
  ConfigFile() = default;

  ConfigFile(std::istream& instream, bool testing = false) :
      ConfigFile{parse_config_file(instream)} {
    finalize(testing);
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

  auto periods() const -> const PeriodTimes&;
  void periods(const PeriodTimes&);

  void read_period_matrix_files() {
    for (auto& [_, p] : _period_map) {
      if (p.adjustment_matrix_filename) {
        p.adjustment_matrix =
            AdjustmentMatrix{p.adjustment_matrix_filename.value(), _area_names};
      }
    }
  }

  [[nodiscard]] bool validate_periods() const {
    bool OllKorrect = true;
    for (auto p : _period_params) {
      if (lagrange_popcount(p.include_area_mask) > 1) {
        LOG_ERROR(
            "Include area mask for period '{}' has more than one set region. "
            "This is invalid.",
            p.name);
        OllKorrect = false;
      }
    }
    return OllKorrect;
  }

  [[nodiscard]] bool finalize_periods() {
    bool OllKorrect = true;
    read_period_matrix_files();
    OllKorrect &= setup_periods();
    finalize_periods(region_count(), max_areas());
    OllKorrect &= validate_periods();
    return OllKorrect;
  }

  void finalize_periods(size_t regions, size_t max_areas) {
    _periods_times.setRegionCount(regions);
    _periods_times.setMaxAreas(max_areas);
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

  auto period_count() const -> size_t;
  auto period_params() -> std::vector<PeriodParams>;
  auto period_params() const -> std::vector<PeriodParams>;
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
  auto cleanTreeFilename() const -> std::filesystem::path;
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

  auto lwrOutputThreshold() const -> double;
  void lwrOutputThreshold(double);

  friend ConfigFileParser;

 private:
  static auto parse_line(ConfigFileParser& lexer,
                         ConfigFile& config,
                         size_t line_number) -> ParsingResult<void>;
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
    if (!_workers.has_value()) { workers(1); }
    if (!_threads_per_worker.has_value()) { _threads_per_worker = 1; }
  }

  void setup_log() const {
    if (!log_filename().empty()) {
      logger::get_log_states().add_file_stream(
          log_filename().c_str(),
          INFO | IMPORTANT | WARNING | ERROR | PROGRESS | DEBUG);
    }
  }

  /*
   * This takes the period map, and turns it into a list of periods. The period
   * map is from the period name, as given by the user, to a "period map entry".
   * This period map entry has the following information:
   *
   * - start and end times,
   * - an adjustment matrix, and
   * - range constraints.
   *
   * What I need to check here is that there are no "gaps" in the periods given
   * by the user. To do this I will need to check that, after ordering by
   * _start_ time:
   *
   * - every end time is _exactly equal_ to an end time;
   * - that no periods are "contianed" within another period;
   * - that there is a start period, which has no specified start time; and
   * - that there is an end period, which has no specified end time.
   */
  bool setup_periods() {
    if (_period_map.empty()) { return true; }
    /* we can't just use value type here, because the key is const........ */
    using PeriodMapEntry =
        std::pair<PeriodConfigMap::key_type, PeriodConfigMap::mapped_type>;

    std::vector<PeriodMapEntry> period_buffer{_period_map.begin(),
                                              _period_map.end()};
    std::sort(period_buffer.begin(),
              period_buffer.end(),
              [](const auto& a, const auto& b) -> bool {
                return a.second.start < b.second.start;
              });

    bool OllKorrect = true;

    if (period_buffer.front().second.start != 0.0) { OllKorrect &= false; }
    if (!std::isinf(period_buffer.back().second.end)) { OllKorrect &= false; }

#if defined(__cpp_lib_ranges_zip) \
    && (!defined(__clang_major__) && __clang_major__ > 18)
    for (auto [a, b] : period_buffer | std::views::adjacent<2>) {
      auto pc_a = a.second;
      auto pc_b = b.second;
#else
    for (auto a_itr = period_buffer.begin(),
              b_itr = next(period_buffer.begin());
         b_itr != period_buffer.end();
         ++a_itr, ++b_itr) {
      auto a = *a_itr;
      auto b = *b_itr;
      auto pc_a = a_itr->second;
      auto pc_b = b_itr->second;
#endif
      if (pc_a.end != pc_b.start) {
        OllKorrect &= false;
        LOG_ERROR(
            "Failed to make periods, end of '{}' does not match start of '{}' "
            "({} vs {})",
            a.first,
            b.first,
            pc_a.end,
            pc_b.start);
      }

      if (std::isinf(pc_a.end)) {
        OllKorrect &= false;
        LOG_ERROR(
            "There are several end periods, please specify only one end "
            "period.");
      }

      if (pc_b.start == 0.0) {
        OllKorrect &= false;
        LOG_ERROR(
            "There are several start periods, please specify only one start "
            "period.");
      }
    }

    if (!OllKorrect) { return false; }

    _period_params =
        period_buffer
        | std::views::transform([&](const auto& a) -> PeriodParams {
            auto [key, p] = a;
            return {
                .dispersion_rate = p.dispersion.value_or(0.01),
                .extinction_rate = p.extinction.value_or(0.01),
                .distance_penalty = 1.0,
                .adjustment_matrix = p.adjustment_matrix
                                         ? p.adjustment_matrix->to_matrix()
                                         : nullptr,
                .name = key,
                .regions = _region_count.value_or(0),
                .include_area_mask =
                    p.include_areas
                        ? convert_dist_binary_string_to_dist(*p.include_areas)
                        : 0,
                .exclude_area_mask =
                    p.exclude_areas
                        ? convert_dist_binary_string_to_dist(*p.exclude_areas)
                        : 0,
            };
          })
        | std::ranges::to<std::vector<PeriodParams>>();

    auto period_times = period_buffer
                        | std::views::take(period_buffer.size() - 1)
                        | std::views::transform([](const auto& a) -> double {
                            return a.second.end;
                          });

    _periods_times = PeriodTimes(period_times);

    return OllKorrect;
  }

  void setup_regions() { _region_count = _area_names.size(); }

  void finalize(bool testing = false) {
    if (!testing) { setup_log(); }
    setup_regions();
    set_threads();
    set_mrcas_for_fossils();
    check_prefix();
    validate_and_make_prefix();
  }

  std::filesystem::path _tree_filename;
  std::filesystem::path _data_filename;
  std::filesystem::path _log_filename;
  std::filesystem::path _prefix;

  std::optional<OutputType> _output_file_type;

  std::optional<std::filesystem::path> _adjustment_matrix_filename;

  std::optional<size_t> _max_areas;

  std::vector<std::string> _area_names;

  PeriodConfigMap _period_map;
  std::vector<PeriodParams> _period_params;
  PeriodTimes _periods_times;

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
  double _output_threshold = 1e-6;

  std::optional<LagrangeEXPMComputationMode> _expm_mode;

  std::optional<AlignmentFileType> _alignment_file_type;

  std::optional<size_t> _region_count;
  std::optional<size_t> _workers;
  std::optional<size_t> _threads_per_worker;
  std::optional<LagrangeOperationMode> _run_mode{
      LagrangeOperationMode::OPTIMIZE};
  std::optional<OptimizationMethod> _opt_method{OptimizationMethod::BOBYQA};

  std::optional<bool> _allow_ambigious;
  std::optional<bool> _dump_graph;
};

}  // namespace lagrange
#endif
