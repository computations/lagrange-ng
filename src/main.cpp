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

  int maxareas = 0;

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
  bool sparse = false;

  double dispersal = 0.01;
  double extinction = 0.01;
  bool estimate = true;

  std::vector<std::vector<lagrange_dist_t>> stochastic_number_from_tos;
  std::vector<lagrange_dist_t> stochastic_time_dists;

  // estimating the dispersal mask
  bool estimate_dispersal_mask = false;

  size_t region_count;
  lagrange_option_t<size_t> workers;
  lagrange_option_t<size_t> threads_per_worker;
};

std::unique_ptr<lagrange_matrix_base_t[]> normalizeStateDistrubtionByLWR(
    const std::unique_ptr<lagrange_matrix_base_t[]> &states,
    size_t states_len) {
  std::unique_ptr<lagrange_matrix_base_t[]> normalized_states{
      new lagrange_matrix_base_t[states_len]};

  for (size_t i = 0; i < states_len; i++) {
    normalized_states.get()[i] = states.get()[i];
  }

  double max_llh = -std::numeric_limits<double>::infinity();

  if (states_len == 1) { throw std::runtime_error{"YOU FUCKED UP BEN"}; }

  for (size_t i = 0; i < states_len; ++i) {
    double tmp = normalized_states.get()[i];

    max_llh = std::max(max_llh, tmp);
  }

  double total_llh = 0.0;
  for (size_t i = 0; i < states_len; i++) {
    total_llh += std::exp(normalized_states.get()[i] - max_llh);
  }

  for (size_t i = 0; i < states_len; i++) {
    double tmp = std::exp(normalized_states.get()[i] - max_llh) / total_llh;

    normalized_states.get()[i] = tmp;
  }

  return normalized_states;
}

std::vector<std::string> grab_token(const std::string &token,
                                    const std::string &deliminators) {
  std::vector<std::string> searchtokens;
  Tokenize(token, searchtokens, deliminators);
  for (unsigned int j = 0; j < searchtokens.size(); j++) {
    TrimSpaces(searchtokens[j]);
  }
  return searchtokens;
}

config_options_t parse_config(const std::string &config_filename) {
  config_options_t config;

  std::ifstream ifs(config_filename);
  if (!ifs) {
    throw std::runtime_error{"Could not open the config file for parsing"};
  }

  std::string line;
  while (getline(ifs, line)) {
    if (line.size() == 0 || line[0] == '#') { continue; }
    /* Make the token */
    std::vector<std::string> tokens = grab_token(line, "=");

    /* Parse the option in the token */
    if (!strcmp(tokens[0].c_str(), "treefile")) {
      config.treefile = tokens[1];
    } else if (!strcmp(tokens[0].c_str(), "datafile")) {
      config.datafile = tokens[1];
    } else if (!strcmp(tokens[0].c_str(), "ratematrix")) {
      config.ratematrixfile = tokens[1];
      if (config.ratematrixfile == "d" || config.ratematrixfile == "D") {
        config.ratematrixfile = "";
      }
    } else if (!strcmp(tokens[0].c_str(), "areanames")) {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      config.areaNames = searchtokens;
    } else if (!strcmp(tokens[0].c_str(), "fixnode")) {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      std::vector<int> dist;
      for (unsigned int j = 0; j < searchtokens[1].size(); j++) {
        char c = (searchtokens[1].c_str())[j];
        dist.push_back(atoi(&c));
      }
      config.fixnodewithmrca[searchtokens[0]] =
          convert_vector_to_lagrange_dist(dist);
    } else if (!strcmp(tokens[0].c_str(), "excludedists")) {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      for (unsigned int j = 0; j < searchtokens.size(); j++) {
        std::vector<int> dist;
        for (unsigned int k = 0; k < searchtokens[j].size(); k++) {
          char c = (searchtokens[j].c_str())[k];
          dist.push_back(atoi(&c));
        }
        config.excludedists.push_back(convert_vector_to_lagrange_dist(dist));
      }
    } else if (!strcmp(tokens[0].c_str(), "includedists")) {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      if (searchtokens[0].size() == 1) {
        config.maxareas = atoi(searchtokens[0].c_str());
      } else {
        for (unsigned int j = 0; j < searchtokens.size(); j++) {
          std::vector<int> dist;
          for (unsigned int k = 0; k < searchtokens[j].size(); k++) {
            char c = (searchtokens[j].c_str())[k];
            dist.push_back(atoi(&c));
          }
          config.includedists.push_back(convert_vector_to_lagrange_dist(dist));
        }
      }
    } else if (!strcmp(tokens[0].c_str(), "areacolors")) {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      config.areacolors = searchtokens;
    } else if (!strcmp(tokens[0].c_str(), "periods")) {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      for (unsigned int j = 0; j < searchtokens.size(); j++) {
        config.periods.push_back(atof(searchtokens[j].c_str()));
      }
    } else if (!strcmp(tokens[0].c_str(), "mrca")) {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      std::vector<std::string> mrc;
      for (unsigned int j = 1; j < searchtokens.size(); j++) {
        mrc.push_back(searchtokens[j]);
      }
      config.mrcas[searchtokens[0]] = mrc;
    } else if (!strcmp(tokens[0].c_str(), "ancstate")) {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      config.ancstates.push_back(searchtokens[0]);
    } else if (!strcmp(tokens[0].c_str(), "fossil")) {
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      config.fossiltype.push_back(searchtokens[0]);
      config.fossilmrca.push_back(searchtokens[1]);
      config.fossilarea.push_back(searchtokens[2]);
      if (searchtokens.size() > 3) {
        config.fossilage.push_back(atof(searchtokens[3].c_str()));
      } else {
        config.fossilage.push_back(0.0);
      }
    } else if (!strcmp(tokens[0].c_str(), "calctype")) {
      std::string calctype = tokens[1];
      if (calctype.compare("m") != 0 && calctype.compare("M") != 0) {
        config.marginal = false;
      }
    } else if (!strcmp(tokens[0].c_str(), "report")) {
      if (tokens[1].compare("split") != 0) { config.splits = false; }
    } else if (!strcmp(tokens[0].c_str(), "sparse")) {
      config.sparse = true;
    } else if (!strcmp(tokens[0].c_str(), "splits")) {
      config.splits = true;
    } else if (!strcmp(tokens[0].c_str(), "states")) {
      config.states = true;
    } else if (!strcmp(tokens[0].c_str(), "estimate_dispersal_mask")) {
      config.estimate_dispersal_mask = true;
    } else if (!strcmp(tokens[0].c_str(), "stochastic_time")) {
      config.states = true;  // requires ancestral states
      if (config.ancstates.size() > 0) {
        config.ancstates[0] = "_all_";
      } else {
        config.ancstates.push_back("_all_");
      }
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      for (unsigned int j = 0; j < searchtokens.size(); j++) {
        std::vector<int> dist;
        for (unsigned int k = 0; k < searchtokens[j].size(); k++) {
          char c = (searchtokens[j].c_str())[k];
          dist.push_back(atoi(&c));
        }
        config.stochastic_time_dists.push_back(
            convert_vector_to_lagrange_dist(dist));
      }
    } else if (!strcmp(tokens[0].c_str(), "stochastic_number")) {
      config.states = true;  // requires ancestral states
      if (config.ancstates.size() > 0) {
        config.ancstates[0] = "_all_";
      } else {
        config.ancstates.push_back("_all_");
      }
      std::vector<std::string> searchtokens = grab_token(tokens[1], ",     ");
      if (searchtokens.size() != 2) {
        std::cout << "ERROR: distributions for stochastic_number need to be "
                     "in the form from_to"
                  << std::endl;
      } else {
        std::vector<lagrange_dist_t> dists;
        std::vector<int> dist0;
        for (unsigned int k = 0; k < searchtokens[0].size(); k++) {
          char c = (searchtokens[0].c_str())[k];
          dist0.push_back(atoi(&c));
        }
        std::vector<int> dist1;
        for (unsigned int k = 0; k < searchtokens[1].size(); k++) {
          char c = (searchtokens[1].c_str())[k];
          dist1.push_back(atoi(&c));
        }
        dists.push_back(convert_vector_to_lagrange_dist(dist0));
        dists.push_back(convert_vector_to_lagrange_dist(dist1));
        config.stochastic_number_from_tos.push_back(dists);
      }
    } else if (!strcmp(tokens[0].c_str(), "dispersal")) {
      config.dispersal = atof(tokens[1].c_str());
      std::cout << "setting dispersal: " << config.dispersal << std::endl;
      config.estimate = false;
    } else if (!strcmp(tokens[0].c_str(), "extinction")) {
      config.extinction = atof(tokens[1].c_str());
      std::cout << "setting extinction: " << config.extinction << std::endl;
      config.estimate = false;
    } else if (!strcmp(tokens[0].c_str(), "workers")) {
      config.workers = std::stoi(tokens[1]);
    } else if (!strcmp(tokens[0].c_str(), "threads-per-worker")) {
      config.threads_per_worker = std::stoi(tokens[1]);
    } else if (!strcmp(tokens[0].c_str(), "maxaraes")) {
      config.maxareas = std::stoi(tokens[1]);
    }
  }
  ifs.close();
  return config;
}

nlohmann::json makeStateJsonOutput(
    const std::vector<std::unique_ptr<lagrange_matrix_base_t[]>> &states,
    size_t states_len, const std::vector<size_t> &stateToIdMap) {
  nlohmann::json states_json;
  for (size_t i = 0; i < states.size(); ++i) {
    nlohmann::json node_json;
    node_json["number"] = stateToIdMap[i];
    auto &state_distribution = states[i];
    auto lwr_distribution =
        normalizeStateDistrubtionByLWR(state_distribution, states_len);
    for (size_t dist = 0; dist < states_len; ++dist) {
      nlohmann::json tmp;
      tmp["distribution"] = dist;
      double llh = state_distribution.get()[dist];
      tmp["llh"] = llh;
      double ratio = lwr_distribution.get()[dist];
      tmp["ratio"] = ratio;
      node_json["states"].push_back(tmp);
    }
    states_json.push_back(node_json);
  }
  return states_json;
}

void writeJsonToFile(const config_options_t &config,
                     const nlohmann ::json &root_json) {
  std::string json_filename = config.treefile + ".results.json";
  std::ofstream outfile(json_filename);
  outfile << root_json.dump();
}

void writeResultTree(const config_options_t &config,
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

void handle_tree(std::shared_ptr<Tree> intree,
                 const std::unordered_map<std::string, lagrange_dist_t> data,
                 const config_options_t &config) {
  nlohmann::json root_json;
  nlohmann::json attributes_json;
  attributes_json["periods"] =
      config.periods.size() ? config.periods.size() != 0 : 1;
  attributes_json["regions"] = config.region_count;
  attributes_json["taxa"] = intree->getExternalNodeCount();
  root_json["attributes"] = attributes_json;
  Context context(intree, config.region_count, config.maxareas);
  context.registerLHGoal();
  if (config.states) { context.registerStateLHGoal(); }
  context.init();
  context.updateRates({config.dispersal, config.extinction});
  context.registerTipClvs(data);

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
                         config.splits, true);
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
    root_json["node-results"] = makeStateJsonOutput(
        states, context.stateCount(), stateGoalIndexToIdMap);
  }
  if (config.splits) { auto splits = context.getStateResults(); }
  // std::cout << context.treeCLVStatus() << std::endl;
  writeJsonToFile(config, root_json);
}

void setThreads(config_options_t &config) {
  if (!config.workers.has_value()) { config.workers = 1; }
  if (!config.threads_per_worker.has_value()) { config.threads_per_worker = 1; }
}

int main(int argc, char *argv[]) {
  openblas_set_num_threads(1);
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
    std::vector<std::shared_ptr<Tree>> intrees;
    ir.readMultipleTreeFile(config.treefile, intrees);
    std::cout << "reading data..." << std::endl;
    std::unordered_map<std::string, lagrange_dist_t> data =
        ir.readStandardInputData(config.datafile);
    std::cout << "checking data..." << std::endl;
    ir.checkData(data, intrees);

    config.region_count = ir.nareas;
    if (config.maxareas == 0) { config.maxareas = config.region_count; }

    for (unsigned int i = 0; i < intrees.size(); i++) {
      handle_tree(intrees[i], data, config);
    }
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end_time - start_time;
  std::cout << "Analysis took: " << duration.count() << "s" << std::endl;
  return 0;
}
