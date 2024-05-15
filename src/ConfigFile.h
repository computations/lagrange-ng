#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include "Alignment.h"
#include "Common.h"
#include "Fossil.h"
#include "Periods.hpp"
#include "Utils.h"

constexpr size_t KRYLOV_RANGE_COUNT_THRESHOLD = 5;

struct ConfigFile {
  std::string treefile;
  std::string datafile;
  std::string ratematrixfile;
  std::string logfile;
  std::filesystem::path prefix;

  size_t maxareas = 0;

  Periods periods;
  std::unordered_map<std::string, std::shared_ptr<MRCAEntry>> mrcas;
  std::vector<lagrange_dist_t> excludedists;
  std::vector<lagrange_dist_t> includedists;
  std::vector<std::string> areaNames;
  std::unordered_map<std::string, lagrange_dist_t> areaNameToDistMap;
  std::vector<std::string> ancstates;
  std::vector<std::string> areacolors;

  std::vector<Fossil> fossils;
  std::unordered_map<std::string, lagrange_dist_t> fixNodeWithMrca;

  bool marginal = true;  // false means joint
  bool splits = false;
  bool states = false;

  double dispersal = 0.01;
  double extinction = 0.01;
  bool estimate = true;

  double lh_epsilon = 1e-9;

  lagrange_option_t<lagrange_expm_computation_mode> expm_mode;

  lagrange_option_t<AlignmentFileType> alignment_file_type;

  size_t region_count{};
  lagrange_option_t<size_t> workers;
  lagrange_option_t<size_t> threads_per_worker;
  lagrange_option_t<lagrange_operation_mode> mode{
      lagrange_operation_mode::OPTIMIZE};
  lagrange_option_t<std::pair<double, double>> params;

  std::filesystem::path get_results_filename() const;
  std::filesystem::path get_node_tree_filename() const;
  std::filesystem::path get_scaled_tree_filename() const;
};

ConfigFile parse_config_file(std::istream& instream);

#endif
