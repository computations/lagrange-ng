/*
 * main.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: Stephen A. Smith
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include <cctype>
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
#include "ConfigFile.h"
#include "Context.h"
#include "Fossil.h"
#include "InputReader.h"
#include "Utils.h"
#include "WorkerState.h"
#include "nlohmann/json.hpp"

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

static void writeJsonToFile(const ConfigFile &config,
                            const nlohmann ::json &root_json) {
  std::string json_filename = config.treefile + ".results.json";
  std::ofstream outfile(json_filename);
  outfile << root_json.dump();
  // nlohmann::json::to_cbor(root_json, outfile);
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

static void set_expm_mode(Context &context, const ConfigFile &config) {
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
    const ConfigFile &config) {
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
  attributes_json["max-areas"] = config.maxareas;
  attributes_json["state-count"] = lagrange_compute_restricted_state_count(
      config.region_count, config.maxareas);
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

static void setThreads(ConfigFile &config) {
  if (!config.workers.has_value()) { config.workers = 1; }
  if (!config.threads_per_worker.has_value()) { config.threads_per_worker = 1; }
}

static ConfigFile read_config_file(const std::string &filename) {
  std::ifstream infile(filename);
  return parse_config_file(infile);
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
    auto config = read_config_file(config_filename);
    setThreads(config);

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
