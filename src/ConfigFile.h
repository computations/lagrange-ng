#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "Common.h"
#include "Fossil.h"
#include "Utils.h"

constexpr size_t KRYLOV_RANGE_COUNT_THRESHOLD = 5;

struct ConfigFile {
  std::string treefile;
  std::string datafile;
  std::string ratematrixfile;
  std::string logfile;

  size_t maxareas = 0;

  std::vector<double> periods;
  std::unordered_map<std::string, std::vector<std::string>> mrcas;
  std::unordered_map<std::string, lagrange_dist_t> fixnodewithmrca;
  std::vector<lagrange_dist_t> excludedists;
  std::vector<lagrange_dist_t> includedists;
  std::vector<std::string> areaNames;
  std::unordered_map<std::string, lagrange_dist_t> areaNameToDistMap;
  std::vector<std::string> ancstates;
  std::vector<std::string> areacolors;

  std::vector<Fossil> fossils;

  std::vector<std::string> fossilmrca;
  std::vector<std::string> fossiltype;
  std::vector<std::string> fossilarea;
  std::vector<double> fossilage;

  bool marginal = true;  // false means joint
  bool splits = false;
  bool states = false;

  double dispersal = 0.01;
  double extinction = 0.01;
  bool estimate = true;

  double lh_epsilon = 1e-9;

  lagrange_option_t<lagrange_expm_computation_mode> expm_mode;

  size_t region_count{};
  lagrange_option_t<size_t> workers;
  lagrange_option_t<size_t> threads_per_worker;
  lagrange_option_t<lagrange_operation_mode> mode{
      lagrange_operation_mode::OPTIMIZE};
  lagrange_option_t<std::pair<double, double>> params;
};

ConfigFile parse_config_file(std::istream& instream);

#endif
