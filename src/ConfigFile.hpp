#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include "Alignment.hpp"
#include "Common.hpp"
#include "Fossil.hpp"
#include "Periods.hpp"
#include "Utils.hpp"

constexpr size_t KRYLOV_RANGE_COUNT_THRESHOLD = 5;

struct ConfigFile {
  std::string tree_filename;
  std::string data_filename;
  std::string rate_matrix_filename;
  std::string log_filename;
  std::filesystem::path prefix;

  size_t maxareas = 0;

  Periods periods;
  std::unordered_map<std::string, std::shared_ptr<MRCAEntry>> mrcas;
  std::vector<lagrange_dist_t> exclude_dists;
  std::vector<lagrange_dist_t> include_dists;
  std::vector<std::string> area_names;
  std::vector<std::string> anc_states;

  std::vector<Fossil> fossils;
  std::unordered_map<std::string, lagrange_dist_t> fix_node_with_mrca;

  bool marginal = true;  // false means joint
  bool splits = false;
  bool states = false;

  double dispersal = 0.01;
  double extinction = 0.01;
  bool estimate = true;

  double lh_epsilon = 1e-9;

  LagrangeOption<LagrangeEXPMComputationMode> expm_mode;

  LagrangeOption<AlignmentFileType> alignment_file_type;

  size_t region_count{};
  LagrangeOption<size_t> workers;
  LagrangeOption<size_t> threads_per_worker;
  LagrangeOption<LagrangeOperationMode> run_mode{
      LagrangeOperationMode::OPTIMIZE};
  LagrangeOption<std::pair<double, double>> rate_params;

  std::filesystem::path resultsFilename() const;
  std::filesystem::path NodeTreeFilename() const;
  std::filesystem::path scaledTreeFilename() const;
};

ConfigFile parse_config_file(std::istream& instream);

#endif
