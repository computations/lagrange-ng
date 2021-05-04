/*
 * main.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: Stephen A. Smith
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>

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
#include "ThreadState.h"
#include "Utils.h"
#include "blis/blis.h"
#include "nlohmann/json.hpp"

struct config_options_t {
  std::string treefile;
  std::string datafile;
  std::string ratematrixfile;
  std::string logfile;

  int maxareas = 0;

  vector<double> periods;
  unordered_map<string, vector<string>> mrcas;
  unordered_map<string, lagrange_dist_t> fixnodewithmrca;
  vector<lagrange_dist_t> excludedists;
  vector<lagrange_dist_t> includedists;
  vector<string> areaNames;
  unordered_map<string, lagrange_dist_t> areaNameToDistMap;
  vector<string> ancstates;
  vector<string> areacolors;
  vector<string> fossilmrca;
  vector<string> fossiltype;
  vector<string> fossilarea;
  vector<double> fossilage;

  bool marginal = true;  // false means joint
  bool splits = false;
  bool states = false;
  bool sparse = false;

  double dispersal = 0.01;
  double extinction = 0.01;
  bool estimate = true;

  vector<vector<lagrange_dist_t>> stochastic_number_from_tos;
  vector<lagrange_dist_t> stochastic_time_dists;

  // estimating the dispersal mask
  bool estimate_dispersal_mask = false;

  size_t region_count;
  lagrange_option_t<size_t> threads;
};

std::unique_ptr<lagrange_matrix_base_t> normalizeStateDistrubtionByLWR(
    const std::unique_ptr<lagrange_matrix_base_t> &states) {
  std::unique_ptr<lagrange_matrix_base_t> normalized_states{
      new lagrange_matrix_base_t};

  // bli_printv("States prenormal:", states.get(), "%1.5f", "");
  bli_obj_create_conf_to(states.get(), normalized_states.get());
  bli_copyv(states.get(), normalized_states.get());

  double max_llh = -std::numeric_limits<double>::infinity();
  size_t states_size = bli_obj_length(normalized_states.get());

  if (states_size == 1) { throw std::runtime_error{"YOU FUCKED UP BEN"}; }

  for (size_t i = 0; i < states_size; ++i) {
    double tmp = bli_getiv_real(normalized_states.get(), i);
    max_llh = std::max(max_llh, tmp);
  }

  double total_llh = 0.0;
  for (size_t i = 0; i < states_size; i++) {
    total_llh += std::exp(bli_getiv_real(normalized_states.get(), i) - max_llh);
  }

  for (size_t i = 0; i < states_size; i++) {
    double tmp =
        std::exp(bli_getiv_real(normalized_states.get(), i) - max_llh) /
        total_llh;

    bli_setiv_real(normalized_states.get(), tmp, i);
  }
  // bli_printv("States postnormal:", normalized_states.get(), "%1.5f", "");

  return normalized_states;
}

std::vector<std::string> grab_token(const std::string &token,
                                    const std::string &deliminators) {
  vector<string> searchtokens;
  Tokenize(token, searchtokens, deliminators);
  for (unsigned int j = 0; j < searchtokens.size(); j++) {
    TrimSpaces(searchtokens[j]);
  }
  return searchtokens;
}

config_options_t parse_config(const std::string &config_filename) {
  config_options_t config;

  ifstream ifs(config_filename);
  if (!ifs) {
    throw std::runtime_error{"Could not open the config file for parsing"};
  }

  std::string line;
  while (getline(ifs, line)) {
    if (line.size() == 0 || line[0] == '#') { continue; }
    /* Make the token */
    vector<string> tokens = grab_token(line, "=");

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
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      config.areaNames = searchtokens;
    } else if (!strcmp(tokens[0].c_str(), "fixnode")) {
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      vector<int> dist;
      for (unsigned int j = 0; j < searchtokens[1].size(); j++) {
        char c = (searchtokens[1].c_str())[j];
        dist.push_back(atoi(&c));
      }
      config.fixnodewithmrca[searchtokens[0]] =
          convert_vector_to_lagrange_dist(dist);
    } else if (!strcmp(tokens[0].c_str(), "excludedists")) {
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      for (unsigned int j = 0; j < searchtokens.size(); j++) {
        vector<int> dist;
        for (unsigned int k = 0; k < searchtokens[j].size(); k++) {
          char c = (searchtokens[j].c_str())[k];
          dist.push_back(atoi(&c));
        }
        config.excludedists.push_back(convert_vector_to_lagrange_dist(dist));
      }
    } else if (!strcmp(tokens[0].c_str(), "includedists")) {
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      if (searchtokens[0].size() == 1) {
        config.maxareas = atoi(searchtokens[0].c_str());
      } else {
        for (unsigned int j = 0; j < searchtokens.size(); j++) {
          vector<int> dist;
          for (unsigned int k = 0; k < searchtokens[j].size(); k++) {
            char c = (searchtokens[j].c_str())[k];
            dist.push_back(atoi(&c));
          }
          config.includedists.push_back(convert_vector_to_lagrange_dist(dist));
        }
      }
    } else if (!strcmp(tokens[0].c_str(), "areacolors")) {
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      config.areacolors = searchtokens;
    } else if (!strcmp(tokens[0].c_str(), "periods")) {
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      for (unsigned int j = 0; j < searchtokens.size(); j++) {
        config.periods.push_back(atof(searchtokens[j].c_str()));
      }
    } else if (!strcmp(tokens[0].c_str(), "mrca")) {
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      vector<string> mrc;
      for (unsigned int j = 1; j < searchtokens.size(); j++) {
        mrc.push_back(searchtokens[j]);
      }
      config.mrcas[searchtokens[0]] = mrc;
    } else if (!strcmp(tokens[0].c_str(), "ancstate")) {
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      config.ancstates.push_back(searchtokens[0]);
    } else if (!strcmp(tokens[0].c_str(), "fossil")) {
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      config.fossiltype.push_back(searchtokens[0]);
      config.fossilmrca.push_back(searchtokens[1]);
      config.fossilarea.push_back(searchtokens[2]);
      if (searchtokens.size() > 3) {
        config.fossilage.push_back(atof(searchtokens[3].c_str()));
      } else {
        config.fossilage.push_back(0.0);
      }
    } else if (!strcmp(tokens[0].c_str(), "calctype")) {
      string calctype = tokens[1];
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
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      for (unsigned int j = 0; j < searchtokens.size(); j++) {
        vector<int> dist;
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
      vector<string> searchtokens = grab_token(tokens[1], ",     ");
      if (searchtokens.size() != 2) {
        cout << "ERROR: distributions for stochastic_number need to be "
                "in the form from_to"
             << endl;
      } else {
        vector<lagrange_dist_t> dists;
        vector<int> dist0;
        for (unsigned int k = 0; k < searchtokens[0].size(); k++) {
          char c = (searchtokens[0].c_str())[k];
          dist0.push_back(atoi(&c));
        }
        vector<int> dist1;
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
      cout << "setting dispersal: " << config.dispersal << endl;
      config.estimate = false;
    } else if (!strcmp(tokens[0].c_str(), "extinction")) {
      config.extinction = atof(tokens[1].c_str());
      cout << "setting extinction: " << config.extinction << endl;
      config.estimate = false;
    } else if (!strcmp(tokens[0].c_str(), "threads")) {
      config.threads = std::stoi(tokens[1]);
    }
  }
  ifs.close();
  return config;
}

nlohmann::json makeStateJsonOutput(
    const std::vector<std::unique_ptr<lagrange_matrix_base_t>> &states,
    const std::vector<size_t> &stateToIdMap) {
  nlohmann::json states_json;
  for (size_t i = 0; i < states.size(); ++i) {
    nlohmann::json node_json;
    node_json["number"] = stateToIdMap[i];
    auto &state_distribution = states[i];
    auto lwr_distribution = normalizeStateDistrubtionByLWR(state_distribution);
    for (size_t dist = 0;
         dist < static_cast<size_t>(bli_obj_length(state_distribution.get()));
         ++dist) {
      nlohmann::json tmp;
      tmp["distribution"] = dist;
      double llh = bli_getiv_real(state_distribution.get(), dist);
      tmp["llh"] = llh;
      double ratio = bli_getiv_real(lwr_distribution.get(), dist);
      tmp["ratio"] = ratio;
      node_json["states"].push_back(tmp);
    }
    states_json.push_back(node_json);
    bli_obj_free(lwr_distribution.get());
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
                 const std::unordered_map<string, lagrange_dist_t> data,
                 const config_options_t &config) {
  nlohmann::json root_json;
  nlohmann::json attributes_json;
  attributes_json["periods"] =
      config.periods.size() ? config.periods.size() != 0 : 1;
  attributes_json["regions"] = config.region_count;
  attributes_json["taxa"] = intree->getExternalNodeCount();
  root_json["attributes"] = attributes_json;
  Context context(intree, config.region_count);
  context.registerLHGoal();
  if (config.states) { context.registerStateLHGoal(); }
  context.init();
  context.updateRates({config.dispersal, config.extinction});
  context.registerTipClvs(data);

  std::vector<ThreadState> thread_states;
  thread_states.reserve(config.threads.get());
  std::vector<std::thread> threads;
  std::cout << "Starting Threads" << std::endl;
  for (size_t i = 0; i < config.threads.get(); i++) {
    std::cout << "Making Thread #" << i << std::endl;
    thread_states.emplace_back();
    threads.emplace_back(&Context::optimizeAndComputeValues, std::ref(context),
                         std::ref(thread_states[i]), config.states,
                         config.splits, true);
  }

  std::cout << "Joining Threads" << std::endl;
  for (auto &t : threads) { t.join(); }
  std::cout << "Joined Threads" << std::endl;

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
        makeStateJsonOutput(states, stateGoalIndexToIdMap);
    for (auto &s : states) { bli_obj_free(s.get()); }
  }
  if (config.splits) { auto splits = context.getStateResults(); }
  // std::cout << context.treeCLVStatus() << std::endl;
  writeJsonToFile(config, root_json);
}

void validateConfig(config_options_t &config) {
  if (!config.threads.has_value()) { config.threads = 1; }
}

int main(int argc, char *argv[]) {
  auto start_time = chrono::high_resolution_clock::now();
  bli_thread_set_num_threads(1);
  if (argc != 2) {
    cout << "you need more arguments." << endl;
    cout << "usage: lagrange configfile" << endl;
    exit(0);
  } else {
    std::string config_filename(argv[1]);
    auto config = parse_config(config_filename);
    validateConfig(config);

    /*****************
     * finish reading the configuration file
     *****************/
    /*
     * after reading the input file
     */
    InputReader ir;
    cout << "reading tree..." << endl;
    vector<std::shared_ptr<Tree>> intrees;
    ir.readMultipleTreeFile(config.treefile, intrees);
    cout << "reading data..." << endl;
    unordered_map<string, lagrange_dist_t> data =
        ir.readStandardInputData(config.datafile);
    cout << "checking data..." << endl;
    ir.checkData(data, intrees);

    config.region_count = ir.nareas;

    for (unsigned int i = 0; i < intrees.size(); i++) {
      handle_tree(intrees[i], data, config);
    }
  }
  auto end_time = chrono::high_resolution_clock::now();
  chrono::duration<double> duration = end_time - start_time;
  cout << "Analysis took: " << duration.count() << "s" << std::endl;
  return 0;
}
