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

constexpr size_t KRYLOV_RANGE_COUNT_THRESHOLD = 5;

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

  double lh_epsilon = 1e-9;

  lagrange_option_t<lagrange_expm_computation_mode> expm_mode;

  size_t region_count{};
  lagrange_option_t<size_t> workers;
  lagrange_option_t<size_t> threads_per_worker;
  lagrange_option_t<lagrange_operation_mode> mode{
      lagrange_operation_mode::OPTIMIZE};
  lagrange_option_t<std::pair<double, double>> params;
};

static auto normalizeSplitDistributionByLWR(lagrange_split_return_t &splits)
    -> void {
  double max_llh = -std::numeric_limits<double>::infinity();
  for (const auto &kv : splits) {
    for (const auto &sp : kv.second) {
      max_llh = std::max(max_llh, sp.getLikelihood());
    }
  }

  double total_llh = 0.0;
  for (const auto &kv : splits) {
    for (const auto &sp : kv.second) {
      total_llh += std::exp(sp.getLikelihood() - max_llh);
    }
  }

  for (auto &kv : splits) {
    for (auto &sp : kv.second) {
      sp.setLWR(std::exp(sp.getLikelihood() - max_llh) / total_llh);
    }
  }
}

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
    } else if (tokens[0] == "expm-mode") {
      if (tokens[1] == "krylov") {
        config.expm_mode = lagrange_expm_computation_mode::KRYLOV;
      } else if (tokens[1] == "pade") {
        config.expm_mode = lagrange_expm_computation_mode::PADE;
      } else if (tokens[1] == "adaptive") {
        config.expm_mode = lagrange_expm_computation_mode::ADAPTIVE;
      }
    } else if (tokens[0] == "mode") {
      if (tokens[1] == "optimize") {
        config.mode = lagrange_operation_mode::OPTIMIZE;
      } else if (tokens[1] == "evaluate") {
        config.mode = lagrange_operation_mode::EVALUATE;
      }
    }
  }

  ifs.close();
  return config;
}

static auto makeStateResultNode(
    const std::unique_ptr<lagrange_matrix_base_t[]> &state_distribution,
    const std::vector<std::string> &region_names, size_t states_len,
    size_t max_areas) -> nlohmann::json {
  nlohmann::json node_json;
  auto lwr_distribution =
      normalizeStateDistrubtionByLWR(state_distribution, states_len);
  for (size_t dist = 0, dist_index = 0; dist_index < states_len;
       ++dist_index, dist = next_dist(dist, max_areas)) {
    nlohmann::json tmp;
    tmp["distribution"] = dist;
    double llh = state_distribution.get()[dist_index];
    tmp["distribution-string"] =
        lagrange_convert_dist_string(dist, region_names);
    tmp["llh"] = llh;
    double ratio = lwr_distribution.get()[dist_index];
    tmp["ratio"] = ratio;
    tmp["regions"] = lagrange_convert_dist_to_list(dist, region_names);
    node_json.push_back(tmp);
  }
  return node_json;
}

static auto makeSplitResultNode(lagrange_split_return_t &splits,
                                const std::vector<std::string> &region_names,
                                size_t states_len, size_t max_areas)
    -> nlohmann::json {
  nlohmann::json node_json;
  normalizeSplitDistributionByLWR(splits);
  for (size_t dist = 0, dist_index = 0; dist_index < states_len;
       ++dist_index, dist = next_dist(dist, max_areas)) {
    for (const auto &sp : splits[dist]) {
      nlohmann::json anc_json;
      nlohmann::json left_json;
      nlohmann::json right_json;

      assert(dist == sp.anc_dist);
      anc_json["distribution"] = dist;
      anc_json["distribution-string"] =
          lagrange_convert_dist_string(dist, region_names);
      anc_json["regions"] = lagrange_convert_dist_to_list(dist, region_names);

      left_json["distribution"] = sp.l_dist;
      left_json["distribution-string"] =
          lagrange_convert_dist_string(sp.l_dist, region_names);
      left_json["regions"] =
          lagrange_convert_dist_to_list(sp.l_dist, region_names);

      right_json["distribution"] = sp.r_dist;
      right_json["distribution-string"] =
          lagrange_convert_dist_string(sp.r_dist, region_names);
      right_json["regions"] =
          lagrange_convert_dist_to_list(sp.r_dist, region_names);

      nlohmann::json tmp;
      tmp["anc-dist"] = anc_json;
      tmp["left-dist"] = left_json;
      tmp["right-dist"] = right_json;
      tmp["llh"] = sp.getLikelihood();
      tmp["ratio"] = sp.getLWR();
      node_json.push_back(tmp);
    }
  }

  return node_json;
}

static auto makeNodeReultsJSONOutput(
    const std::vector<std::unique_ptr<lagrange_matrix_base_t[]>> &states,
    lagrange_split_list_t &splits, const std::vector<size_t> &state_id_map,
    const std::vector<std::string> &region_names, size_t states_len,
    size_t max_areas) -> nlohmann::json {
  nlohmann::json root_json;

  size_t node_count = std::max(states.size(), splits.size());

  for (size_t i = 0; i < node_count; ++i) {
    nlohmann::json node_json;

    if (i < states.size()) {
      node_json["states"] =
          makeStateResultNode(states[i], region_names, states_len, max_areas);
    }
    if (i < splits.size()) {
      node_json["splits"] =
          makeSplitResultNode(splits[i], region_names, states_len, max_areas);
    }

    node_json["number"] = state_id_map[i];
    root_json.push_back(node_json);
  }

  return root_json;
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

static void set_expm_mode(Context &context, const config_options_t &config) {
  auto &expm_mode = config.expm_mode;
  if (expm_mode.has_value()) {
    switch (expm_mode.get()) {
      case lagrange_expm_computation_mode::ADAPTIVE:
        if (config.region_count > KRYLOV_RANGE_COUNT_THRESHOLD) {
          std::cout << "Enabling adaptive expm computation" << std::endl;
          context.useArnoldi();
        } else {
          std::cout << "Using Pade's method expm computation" << std::endl;
          context.useArnoldi(false, false);
        }
        break;
      case lagrange_expm_computation_mode::PADE:
        std::cout << "Using Pade's method expm computation" << std::endl;
        context.useArnoldi(false, false);
        break;
      case lagrange_expm_computation_mode::KRYLOV:
        std::cout << "Using Krylov subspace based method for expm computation"
                  << std::endl;
        context.useArnoldi(true, false);
        break;
      default:
        throw std::runtime_error{"Unknown Expm Mode"};
    }
  } else {
    if (config.region_count > KRYLOV_RANGE_COUNT_THRESHOLD) {
      std::cout << "Enabling adaptive expm computation" << std::endl;
      context.useArnoldi();
    } else {
      std::cout << "Using Pade's method expm computation" << std::endl;
    }
  }
}

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
  attributes_json["nodes-tree"] =
      intree->getNewickLambda([](const Node &n) -> std::string {
        if (n.isInternal()) {
          return std::to_string(n.getNumber()) + ":" +
                 std::to_string(n.getBL());
        } else {
          return n.getName() + ":" + std::to_string(n.getBL());
        }
      });
  root_json["attributes"] = attributes_json;
  Context context(intree, config.region_count, config.maxareas);
  context.registerLHGoal();
  if (config.states) { context.registerStateLHGoal(); }
  if (config.splits) { context.registerSplitLHGoal(); }
  context.init();
  context.updateRates({config.dispersal, config.extinction});
  context.registerTipClvs(data);
  context.set_lh_epsilon(config.lh_epsilon);
  set_expm_mode(context, config);

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

  auto states = context.getStateResults();
  auto splits = context.getSplitResults();
  root_json["node-results"] = makeNodeReultsJSONOutput(
      states, splits, stateGoalIndexToIdMap, config.areaNames,
      context.stateCount(), config.maxareas);
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
