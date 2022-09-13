/*
 * main.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: Stephen A. Smith
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "Common.h"
#include "Context.h"
#include "InputReader.h"
#include "Utils.h"
#include "WorkerState.h"
#include "nlohmann/json.hpp"

struct config_options_t {
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

  double lh_epsilon = 1e-8;

  bool expm_approximate = false;

  size_t region_count{};
  lagrange_option_t<size_t> workers;
  lagrange_option_t<size_t> threads_per_worker;
  lagrange_option_t<lagrange_mode> mode{lagrange_mode::OPTIMIZE};
  lagrange_option_t<std::pair<double, double>> params;
};

static auto normalizeStateDistrubtionByLWR(
    const std::unique_ptr<lagrange_matrix_base_t[]> &states, size_t states_len)
    -> std::unique_ptr<lagrange_matrix_base_t[]> {
  std::unique_ptr<lagrange_matrix_base_t[]> normalized_states{
      new lagrange_matrix_base_t[states_len]};

  for (size_t i = 0; i < states_len; i++) {
    normalized_states.get()[i] = states.get()[i];
  }

  double max_llh = -std::numeric_limits<double>::infinity();

  assert(states_len != 1);

  for (size_t i = 0; i < states_len; ++i) {
    double tmp = normalized_states.get()[i];

    max_llh = std::max(max_llh, tmp);
  }

  double total_llh = 0.0;
  for (size_t i = 1; i < states_len; i++) {
    total_llh += std::exp(normalized_states.get()[i] - max_llh);
  }

  for (size_t i = 0; i < states_len; i++) {
    double tmp = std::exp(normalized_states.get()[i] - max_llh) / total_llh;

    normalized_states.get()[i] = tmp;
  }

  return normalized_states;
}

static auto grab_token(const std::string &token,
                       const std::string &deliminators)
    -> std::vector<std::string> {
  std::vector<std::string> searchtokens;
  Tokenize(token, searchtokens, deliminators);
  for (auto &searchtoken : searchtokens) { TrimSpaces(searchtoken); }
  return searchtokens;
}

static auto parse_config(const std::string &config_filename)
    -> config_options_t {
  config_options_t config;

  std::ifstream ifs(config_filename);
  if (!ifs) {
    throw std::runtime_error{"Could not open the config file for parsing"};
  }

  std::string line;
  while (getline(ifs, line)) {
    if (line.empty() || line[0] == '#') { continue; }
    /* Make the token */
    std::vector<std::string> tokens = grab_token(line, "=");

    /* Parse the option in the token */
    if (tokens[0] == "treefile") {
      config.treefile = tokens[1];
    } else if (tokens[0] == "datafile") {
      config.datafile = tokens[1];
    } else if (tokens[0] == "ratematrix") {
      config.ratematrixfile = tokens[1];
      if (config.ratematrixfile == "d" || config.ratematrixfile == "D") {
        config.ratematrixfile = "";
      }
    } else if (tokens[0] == "areanames") {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      config.areaNames = searchtokens;
    } else if (tokens[0] == "periods") {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      for (auto &searchtoken : searchtokens) {
        config.periods.push_back(atof(searchtoken.c_str()));
      }
    } else if (tokens[0] == "mrca") {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      std::vector<std::string> mrc;
      for (unsigned int j = 1; j < searchtokens.size(); j++) {
        mrc.push_back(searchtokens[j]);
      }
      config.mrcas[searchtokens[0]] = mrc;
    } else if (tokens[0] == "ancstate") {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      config.ancstates.push_back(searchtokens[0]);
    } else if (tokens[0] == "fossil") {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      config.fossiltype.push_back(searchtokens[0]);
      config.fossilmrca.push_back(searchtokens[1]);
      config.fossilarea.push_back(searchtokens[2]);
      if (searchtokens.size() > 3) {
        config.fossilage.push_back(atof(searchtokens[3].c_str()));
      } else {
        config.fossilage.push_back(0.0);
      }
    } else if (tokens[0] == "calctype") {
      std::string calctype = tokens[1];
      if (calctype != "m" && calctype != "M") { config.marginal = false; }
    } else if (tokens[0] == "report") {
      if (tokens[1] != "split") { config.splits = false; }
    } else if (tokens[0] == "splits") {
      config.splits = true;
    } else if (tokens[0] == "states") {
      config.states = true;
    } else if (tokens[0] == "dispersal") {
      config.dispersal = atof(tokens[1].c_str());
      std::cout << "setting dispersal: " << config.dispersal << std::endl;
      config.estimate = false;
    } else if (tokens[0] == "extinction") {
      config.extinction = atof(tokens[1].c_str());
      std::cout << "setting extinction: " << config.extinction << std::endl;
      config.estimate = false;
    } else if (tokens[0] == "lh-epsilon") {
      config.lh_epsilon = stof(tokens[1]);
    } else if (tokens[0] == "workers") {
      config.workers = lagrange_parse_size_t(tokens[1]);
    } else if (tokens[0] == "threads-per-worker") {
      config.threads_per_worker = lagrange_parse_size_t(tokens[1]);
    } else if (tokens[0] == "maxareas") {
      config.maxareas = lagrange_parse_size_t(tokens[1]);
    } else if (tokens[0] == "approximate") {
      config.expm_approximate = true;
    } else if (tokens[0] == "mode") {
      if (tokens[1] == "optimize") {
        config.mode = lagrange_mode::OPTIMIZE;
      } else if (tokens[1] == "evaluate") {
        config.mode = lagrange_mode::EVALUATE;
      }
    }
  }
  ifs.close();
  return config;
}

static auto makeStateJsonOutput(
    const std::vector<std::unique_ptr<lagrange_matrix_base_t[]>> &states,
    size_t states_len, const std::vector<size_t> &stateToIdMap,
    const std::vector<std::string> &regionNames, size_t maxareas)
    -> nlohmann::json {
  nlohmann::json states_json;
  for (size_t i = 0; i < states.size(); ++i) {
    nlohmann::json node_json;
    node_json["number"] = stateToIdMap[i];
    const auto &state_distribution = states[i];
    auto lwr_distribution =
        normalizeStateDistrubtionByLWR(state_distribution, states_len);
    for (size_t dist = 0, dist_index = 0; dist_index < states_len;
         ++dist_index, dist = next_dist(dist, maxareas)) {
      nlohmann::json tmp;
      tmp["distribution"] = dist;
      double llh = state_distribution.get()[dist_index];
      tmp["distribution-string"] =
          lagrange_convert_dist_string(dist, regionNames);
      tmp["llh"] = llh;
      double ratio = lwr_distribution.get()[dist_index];
      tmp["ratio"] = ratio;
      tmp["regions"] = lagrange_convert_dist_to_list(dist, regionNames);
      node_json["states"].push_back(tmp);
    }
    states_json.push_back(node_json);
  }
  return states_json;
}

static void writeJsonToFile(const config_options_t &config,
                            const nlohmann ::json &root_json) {
  std::string json_filename = config.treefile + ".results.json";
  std::ofstream outfile(json_filename);
  outfile << root_json.dump();
}

#if 0
static void writeResultTree(const config_options_t &config,
                            const std::shared_ptr<Tree> &tree) {
  std::string result_tree_filename = config.treefile + ".results.tree";
  std::ofstream outfile(result_tree_filename);
  outfile << tree->getNewickLambda([](const Node &n) -> std::string {
    if (n.isInternal()) {
      std::ostringstream oss;
      oss << "[&&NHX";
      if (!n.getStateString().empty()) {
        oss << ":state=" << n.getStateString();
      }
      if (!n.getSplitString().empty()) {
        oss << ":split=" << n.getSplitString();
      }
      oss << "]";
      return oss.str();
    }
    return n.getName();
  });
}
#endif

static void handle_tree(
    const std::shared_ptr<Tree> &intree,
    const std::unordered_map<std::string, lagrange_dist_t> &data,
    const config_options_t &config) {
  nlohmann::json root_json;
  nlohmann::json attributes_json;
  attributes_json["periods"] =
      !config.periods.empty() ? static_cast<int>(config.periods.size()) : 1;
  attributes_json["regions"] = config.region_count;
  attributes_json["taxa"] = intree->getExternalNodeCount();
  root_json["attributes"] = attributes_json;
  Context context(intree, config.region_count, config.maxareas);
  context.registerLHGoal();
  if (config.states) { context.registerStateLHGoal(); }
  context.init();
  context.updateRates({config.dispersal, config.extinction});
  context.registerTipClvs(data);
  context.set_lh_epsilon(config.lh_epsilon);
  if (config.expm_approximate) {
    context.useArnoldi();
    std::cout << "Enabling Expm Approximation" << std::endl;
  }

  std::vector<WorkerState> worker_states;
  worker_states.reserve(config.workers.get());
  std::vector<std::thread> threads;
  std::cout << "Starting Workers" << std::endl;
  for (size_t i = 0; i < config.workers.get(); i++) {
    std::cout << "Making Worker #" << i + 1 << std::endl;
    worker_states.emplace_back();
    worker_states.back().set_assigned_threads(config.threads_per_worker.get());
    threads.emplace_back(&Context::optimizeAndComputeValues, std::ref(context),
                         std::ref(worker_states[i]), config.states,
                         config.splits, true, std::cref(config.mode.get()));
  }

  std::cout << "Waiting for workers to finish" << std::endl;
  for (auto &t : threads) { t.join(); }

  nlohmann::json params_json;
  auto params = context.currentParams();
  params_json["dispersion"] = params.dispersion_rate;
  params_json["extinction"] = params.extinction_rate;
  root_json["params"] = params_json;

  auto stateGoalIndexToIdMap =
      intree->traversePreorderInternalNodesOnlyNumbers();

  // invert the map
  if (config.states) {
    auto states = context.getStateResults();
    root_json["node-results"] =
        makeStateJsonOutput(states, context.stateCount(), stateGoalIndexToIdMap,
                            config.areaNames, config.maxareas);
  }
  if (config.splits) { auto splits = context.getSplitResults(); }
  writeJsonToFile(config, root_json);
  std::ofstream node_tree(config.treefile + ".nodes.tre");

  node_tree << intree->getNewickLambda([](const Node &n) -> std::string {
    return std::to_string(n.getNumber()) + ":" + std::to_string(n.getBL());
  });
  std::ofstream anal_tree(config.treefile + ".scaled.tre");
  anal_tree << intree->getNewickLambda([](const Node &n) -> std::string {
    return n.getName() + ":" + std::to_string(n.getBL());
  }) << std::endl;
}

static void setThreads(config_options_t &config) {
  if (!config.workers.has_value()) { config.workers = 1; }
  if (!config.threads_per_worker.has_value()) { config.threads_per_worker = 1; }
}

auto main(int argc, char *argv[]) -> int {
#if MKL_ENABLED
  mkl_set_num_threads(1);
#else
  openblas_set_num_threads(1);
#endif
  auto start_time = std::chrono::high_resolution_clock::now();
  if (argc != 2) {
    std::cout << "you need more arguments." << std::endl;
    std::cout << "usage: lagrange configfile" << std::endl;
    exit(0);
  } else {
    std::string config_filename(argv[1]);
    auto config = parse_config(config_filename);
    setThreads(config);

    /*****************
     * finish reading the configuration file
     *****************/
    /*
     * after reading the input file
     */
    InputReader ir;
    std::cout << "reading tree..." << std::endl;
    std::vector<std::shared_ptr<Tree>> intrees =
        InputReader::readMultipleTreeFile(config.treefile);
    std::cout << "reading data..." << std::endl;
    std::unordered_map<std::string, size_t> data =
        ir.readStandardInputData(config.datafile, config.maxareas);
    std::cout << "checking data..." << std::endl;
    InputReader::checkData(data, intrees);

    config.region_count = ir.nareas;
    if (config.maxareas == 0) { config.maxareas = config.region_count; }

    std::cout << "running analysis..." << std::endl;
    for (auto &intree : intrees) { handle_tree(intree, data, config); }
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end_time - start_time;
  std::cout << "Analysis took: " << duration.count() << "s" << std::endl;
  return 0;
}
