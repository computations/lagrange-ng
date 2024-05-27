#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <filesystem>
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

struct ConfigFile {
  std::string tree_filename;
  std::string data_filename;
  std::string rate_matrix_filename;
  std::string log_filename;
  std::filesystem::path prefix;

  size_t maxareas = 0;

  Periods periods;
  MRCAMap mrcas;
  std::vector<std::string> area_names;
  std::vector<std::string> anc_states;

  std::vector<Fossil> fossils;

  bool marginal = true;  // false means joint
  bool all_splits = false;
  bool all_states = false;

  std::vector<MRCALabel> state_nodes;
  std::vector<MRCALabel> split_nodes;

  double dispersal = 0.01;
  double extinction = 0.01;
  bool estimate = true;

  double lh_epsilon = 1e-9;

  Option<LagrangeEXPMComputationMode> expm_mode;

  Option<AlignmentFileType> alignment_file_type;

  size_t region_count{};
  Option<size_t> workers;
  Option<size_t> threads_per_worker;
  Option<LagrangeOperationMode> run_mode{LagrangeOperationMode::OPTIMIZE};
  Option<std::pair<double, double>> rate_params;

  std::filesystem::path resultsFilename() const;
  std::filesystem::path NodeTreeFilename() const;
  std::filesystem::path scaledTreeFilename() const;

  bool computeStates() const;
  bool computeSplits() const;
};

ConfigFile parse_config_file(std::istream& instream);

}  // namespace lagrange
#endif
