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
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

#include "BioGeoTree.h"
#include "BioGeoTreeTools.h"
#include "InputReader.h"
#include "OptimizeBioGeo.h"
#include "OptimizeBioGeoAllDispersal_nlopt.h"
#include "RateMatrixUtils.h"
#include "RateModel.h"
#include "Utils.h"
#include "superdouble.h"
#include "vector_node_object.h"

#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif

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
  vector<string> areanames;
  unordered_map<string, int> areanamemap;
  vector<string> ancstates;
  vector<string> areacolors;
  vector<string> fossilmrca;
  vector<string> fossiltype;
  vector<string> fossilarea;
  vector<double> fossilage;

  bool marginal = true;  // false means joint
  bool splits = false;
  bool states = false;
  int numthreads = 0;
  bool sparse = false;

  double dispersal = 0.1;
  double extinction = 0.1;
  bool estimate = true;

  vector<vector<lagrange_dist_t>> stochastic_number_from_tos;
  vector<lagrange_dist_t> stochastic_time_dists;

  // estimating the dispersal mask
  bool estimate_dispersal_mask = false;
};

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
    if (line.size() == 0 || line[0] == '#') {
      continue;
    }
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
      config.areanames = searchtokens;
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
      if (tokens[1].compare("split") != 0) {
        config.splits = false;
      }
    } else if (!strcmp(tokens[0].c_str(), "sparse")) {
      config.sparse = true;
    } else if (!strcmp(tokens[0].c_str(), "splits")) {
      config.splits = true;
    } else if (!strcmp(tokens[0].c_str(), "states")) {
      config.states = true;
    } else if (!strcmp(tokens[0].c_str(), "estimate_dispersal_mask")) {
      config.estimate_dispersal_mask = true;
    } else if (!strcmp(tokens[0].c_str(), "numthreads")) {
      config.numthreads = atoi(tokens[1].c_str());
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
    }
  }
  ifs.close();
  return config;
}

void handle_ancsplits_ancstates(
    std::shared_ptr<BioGeoTree> bgt, std::shared_ptr<Tree> intree,
    std::shared_ptr<RateModel> rm, unordered_map<int, string> &areanamemaprev,
    const unordered_map<string, std::shared_ptr<Node>> &mrcanodeint,
    const config_options_t &config, Superdouble &outtotlike) {
  BioGeoTreeTools tt;

  if (config.ancstates.size() == 0) {
    return;
  }

  bgt->set_use_stored_matrices(true);

  bgt->prepare_ancstate_reverse();

  std::vector<std::shared_ptr<Node>> node_indicies;

  if (config.ancstates[0] == "_all_" || config.ancstates[0] == "_ALL_") {
    for (unsigned int j = 0; j < intree->getInternalNodeCount(); j++) {
      node_indicies.push_back(intree->getInternalNode(j));
    }
  } else {
    for (unsigned int j = 0; j < config.ancstates.size(); j++) {
      node_indicies.push_back(mrcanodeint.at(config.ancstates[j]));
    }
  }

  nlohmann::json root_json;
  for (auto nid : node_indicies) {
    nlohmann::json current_json;
    current_json["number"] = nid->getNumber();
    if (config.splits) {
      cout << "Ancestral splits for:\t" << nid->getNumber() << endl;
      unordered_map<lagrange_dist_t, vector<AncSplit>> ras =
          bgt->calculate_ancsplit_reverse(nid);
      tt.summarizeSplits(nid, ras, areanamemaprev, rm->get_int_dists_map(),
                         rm->get_num_areas());
      cout << endl;
    }
    if (config.states) {
      nlohmann::json output_json;
      cout << "Ancestral states for:\t" << nid->getNumber() << endl;
      vector<Superdouble> rast = bgt->calculate_ancstate_reverse(nid);
      outtotlike = calculate_vector_Superdouble_sum(rast);
      tt.summarizeAncState(nid, rast, areanamemaprev, rm->get_int_dists_map(),
                           rm->get_num_areas(), output_json);
      cout << endl;

      current_json["states"] = output_json;
    }
    root_json.push_back(current_json);
  }
  ofstream jsonfile((config.treefile + ".bgstates.json").c_str(), ios::app);
  jsonfile << root_json.dump() << std::endl;

  /*
   * key file output
   */
  ofstream outTreeKeyFile;
  outTreeKeyFile.open((config.treefile + ".bgkey.tre").c_str(), ios::app);
  // need to output numbers
  outTreeKeyFile << intree->getRoot()->getNewick(
                        true,
                        [](const Node &n) -> string {
                          return std::to_string(n.getNumber());
                        })
                 << ";" << endl;
  if (config.splits) {
    ofstream outTreeFile;
    outTreeFile.open((config.treefile + ".bgsplits.tre").c_str(), ios::app);
    // need to output object "split"
    outTreeFile << intree->getRoot()->getNewick(true,
                                                [](const Node &n) -> string {
                                                  return n.getSplitString();
                                                })
                << ";" << endl;
  }
  if (config.states) {
    ofstream outTreeFile;
    outTreeFile.open((config.treefile + ".bgstates.tre").c_str(), ios::app);
    // need to output object "state"
    outTreeFile << intree->getRoot()->getNewick(true,
                                                [](const Node &n) -> string {
                                                  return n.getStateString();
                                                })
                << ";" << endl;
  }
}

void handle_stoch(std::shared_ptr<BioGeoTree> bgt, std::shared_ptr<Tree> intree,
                  std::shared_ptr<RateModel> rm,
                  unordered_map<int, string> &areanamemaprev,
                  const config_options_t &config, Superdouble totlike) {
  BioGeoTreeTools tt;
  if (config.stochastic_time_dists.size() > 0) {
    cout << "calculating stochastic mapping time spent" << endl;
    for (unsigned int k = 0; k < config.stochastic_time_dists.size(); k++) {
      cout << tt.get_string_from_dist_int(
                  rm->get_dists_int_map().at(config.stochastic_time_dists[k]),
                  areanamemaprev, rm->get_int_dists_map(), rm->get_num_areas())
           << endl;
      bgt->prepare_stochmap_reverse_all_nodes(
          rm->get_dists_int_map().at(config.stochastic_time_dists[k]),
          rm->get_dists_int_map().at(config.stochastic_time_dists[k]));
      bgt->prepare_ancstate_reverse();
      ofstream outStochTimeFile;
      outStochTimeFile.open((config.treefile + ".bgstochtime.tre").c_str(),
                            ios::app);
      for (unsigned int j = 0; j < intree->getNodeCount(); j++) {
        if (intree->getNode(j) != intree->getRoot()) {
          vector<Superdouble> rsm =
              bgt->calculate_reverse_stochmap(intree->getNode(j), true);
          double stres = calculate_vector_Superdouble_sum(rsm) / totlike;
          intree->getNode(j)->setStochString(std::to_string(stres));
        }
      }
      // need to output object "stoch"
      outStochTimeFile
          << intree->getRoot()->getNewickLambda(
                 [](const Node &n) -> string { return n.getStochString(); })
          << tt.get_string_from_dist_int(
                 rm->get_dists_int_map().at(config.stochastic_time_dists[k]),
                 areanamemaprev, rm->get_int_dists_map(), rm->get_num_areas())
          << ";" << endl;
    }
  }
  if (config.stochastic_number_from_tos.size() > 0) {
    cout << "calculating stochastic mapping number of transitions" << endl;
    for (unsigned int k = 0; k < config.stochastic_number_from_tos.size();
         k++) {
      cout << tt.get_string_from_dist_int(
                  rm->get_dists_int_map().at(
                      config.stochastic_number_from_tos[k][0]),
                  areanamemaprev, rm->get_int_dists_map(), rm->get_num_areas())
           << " -> "
           << tt.get_string_from_dist_int(
                  rm->get_dists_int_map().at(
                      config.stochastic_number_from_tos[k][1]),
                  areanamemaprev, rm->get_int_dists_map(), rm->get_num_areas())
           << endl;
      bgt->prepare_stochmap_reverse_all_nodes(
          rm->get_dists_int_map().at(config.stochastic_number_from_tos[k][0]),
          rm->get_dists_int_map().at(config.stochastic_number_from_tos[k][1]));
      bgt->prepare_ancstate_reverse();
      ofstream outStochTimeFile;
      outStochTimeFile.open((config.treefile + ".bgstochnumber.tre").c_str(),
                            ios::app);
      for (unsigned int j = 0; j < intree->getNodeCount(); j++) {
        if (intree->getNode(j) != intree->getRoot()) {
          vector<Superdouble> rsm =
              bgt->calculate_reverse_stochmap(intree->getNode(j), false);
          // cout << calculate_vector_double_sum(rsm) / totlike << endl;
          double stres = calculate_vector_Superdouble_sum(rsm) / totlike;
          intree->getNode(j)->setStochString(std::to_string(stres));
        }
      }
      // need to output object "stoch"
      outStochTimeFile << intree->getRoot()->getNewickLambda(
                              [](const Node &n) -> string {
                                return n.getStochString();
                              })
                       << ";" << endl;
    }
  }
}

void handle_tree(std::shared_ptr<Tree> intree, std::shared_ptr<RateModel> rm,
                 const std::unordered_map<string, lagrange_dist_t> data,
                 unordered_map<int, string> &areanamemaprev,
                 const config_options_t &config) {
  auto bgt = std::make_shared<BioGeoTree>(intree, config.periods);

  /*
   * record the mrcas
   */
  unordered_map<string, std::shared_ptr<Node>> mrcanodeint;
  for (const auto &it : config.mrcas) {
    // records node by number, should maybe just point to node
    mrcanodeint[it.first] = intree->getMRCA(it.second);
    // tt.getLastCommonAncestor(*intree,nodeIds);
    cout << "Reading mrca: " << it.first << " = ";
    for (unsigned int k = 0; k < it.second.size(); k++) {
      cout << it.second[k] << " ";
    }
    cout << endl;
  }

  /*
   * set fixed nodes
   */
  for (auto &fnit : config.fixnodewithmrca) {
    lagrange_dist_t dista = fnit.second;
    for (unsigned int k = 0; k < rm->getDistsSize(); k++) {
      if (dista != rm->getDists()[k]) {
        bgt->set_excluded_dist(rm->getDists().at(k), mrcanodeint[fnit.first]);
      }
    }
    cout << "fixing " << fnit.first << " = ";
    print_lagrange_dist(fnit.second, rm->get_num_areas());
  }

  cout << "setting default model..." << endl;
  bgt->set_default_model(rm);
  cout << "setting up tips..." << endl;
  bgt->set_tip_conditionals(data);

  /*
   * setting up fossils
   */
  for (unsigned int k = 0; k < config.fossiltype.size(); k++) {
    if (config.fossiltype[k] == "n" || config.fossiltype[k] == "N") {
      bgt->setFossilatNodeByMRCA_id(
          mrcanodeint[config.fossilmrca[k]],
          config.areanamemap.at(config.fossilarea.at(k)));
      cout << "Setting node fossil at mrca: " << config.fossilmrca.at(k)
           << " at area: " << config.fossilarea.at(k) << endl;
    } else {
      bgt->setFossilatBranchByMRCA_id(
          mrcanodeint[config.fossilmrca.at(k)],
          config.areanamemap.at(config.fossilarea.at(k)),
          config.fossilage.at(k));
      cout << "Setting branch fossil at mrca: " << config.fossilmrca.at(k)
           << " at area: " << config.fossilarea.at(k)
           << " at age: " << config.fossilage.at(k) << endl;
    }
  }

  /*
   * initial likelihood calculation
   */
  cout << "starting likelihood calculations" << endl;
  cout << "initial -ln likelihood: "
       << double(bgt->eval_likelihood(config.marginal)) << endl;

  /*
   * optimize likelihood
   */
  Superdouble nlnlike = 0;
  if (config.estimate == true) {
    if (config.estimate_dispersal_mask == false) {
      cout << "Optimizing (simplex) -ln likelihood." << endl;
      OptimizeBioGeo opt(bgt, rm, config.marginal);
      vector<double> disext = opt.optimize_global_dispersal_extinction();

      cout << "dis: " << disext[0] << " ext: " << disext[1] << endl;

      rm->setup_D(disext[0]);
      rm->setup_E(disext[1]);

    } else {  // optimize all the dispersal matrix
      cout << "Optimizing (simplex) -ln likelihood with all dispersal "
              "parameters free."
           << endl;
      vector<double> disextrm =
          optimize_dispersal_extinction_all_nlopt(bgt, rm);
      cout << "dis: " << disextrm[0] << " ext: " << disextrm[1] << endl;
      vector<double> cols(rm->get_num_areas(), 0);
      vector<vector<double>> rows(rm->get_num_areas(), cols);
      vector<vector<vector<double>>> D_mask =
          vector<vector<vector<double>>>(config.periods.size(), rows);
      int count = 2;
      for (unsigned int ii = 0; ii < D_mask.size(); ii++) {
        for (unsigned int jj = 0; jj < D_mask[ii].size(); jj++) {
          D_mask[ii][jj][jj] = 0.0;
          for (unsigned int kk = 0; kk < D_mask[ii][jj].size(); kk++) {
            if (kk == jj) {
              continue;
            }
            D_mask[ii][jj][kk] = disextrm[count];
            count += 1;
          }
        }
      }
      cout << "D_mask" << endl;
      for (unsigned int ii = 0; ii < D_mask.size(); ii++) {
        cout << config.periods.at(ii) << endl;
        cout << "\t";
        for (unsigned int j = 0; j < D_mask[ii].size(); j++) {
          cout << config.areanames[j] << "\t";
        }
        cout << endl;
        for (unsigned int j = 0; j < D_mask[ii].size(); j++) {
          cout << config.areanames[j] << "\t";
          for (unsigned int k = 0; k < D_mask[ii][j].size(); k++) {
            cout << D_mask[ii][j][k] << "\t";
          }
          cout << endl;
        }
        cout << endl;
      }
      rm->setup_D_provided(disextrm[0], D_mask);
      rm->setup_E(disextrm[1]);
    }
  } else {
    rm->setup_D(config.dispersal);
    rm->setup_E(config.extinction);
  }

  rm->setup_Q();
  bgt->update_default_model(rm);
  bgt->set_store_p_matrices(true);
  cout << "final -ln likelihood: "
       << double(bgt->eval_likelihood(config.marginal)) << endl;
  bgt->set_store_p_matrices(false);

  /*
   * ancestral splits calculation
   */
  Superdouble totlike;
  handle_ancsplits_ancstates(bgt, intree, rm, areanamemaprev, mrcanodeint,
                             config, totlike);

  /*
   * stochastic mapping calculations
   * REQUIRES that ancestral calculation be done
   */
  handle_stoch(bgt, intree, rm, areanamemaprev, config, totlike);
  /*
   * end stochastic mapping
   */
}

void print_ratematrix(
    const std::vector<std::vector<std::vector<double>>> &dmconfig) {
  for (unsigned int i = 0; i < dmconfig.size(); i++) {
    for (unsigned int j = 0; j < dmconfig[i].size(); j++) {
      for (unsigned int k = 0; k < dmconfig[i][j].size(); k++) {
        if (dmconfig[i][j][k] != 1) {
          cout << dmconfig[i][j][k];
        } else {
          cout << " . ";
        }
        cout << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
}

void parse_ratematrix_file(const config_options_t &config,
                           const InputReader &ir,
                           std::shared_ptr<RateModel> rm) {
  cout << "Reading rate matrix file" << endl;
  vector<vector<vector<double>>> dmconfig = processRateMatrixConfigFile(
      config.ratematrixfile, ir.nareas, config.periods.size());
  print_ratematrix(dmconfig);
  for (unsigned int i = 0; i < dmconfig.size(); i++) {
    for (unsigned int j = 0; j < dmconfig[i].size(); j++) {
      for (unsigned int k = 0; k < dmconfig[i][j].size(); k++) {
        rm->set_Dmask_cell(i, j, k, dmconfig[i][j][k], false);
      }
    }
  }
}

int main(int argc, char *argv[]) {
  auto start_time = chrono::high_resolution_clock::now();
  if (argc != 2) {
    cout << "you need more arguments." << endl;
    cout << "usage: lagrange configfile" << endl;
    exit(0);
  } else {
    std::string config_filename(argv[1]);
    auto config = parse_config(config_filename);

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

    /*
     * read area names
     */
    unordered_map<int, string> areanamemaprev;
    if (config.areanames.size() > 0) {
      cout << "reading area names" << endl;
      for (unsigned int i = 0; i < config.areanames.size(); i++) {
        config.areanamemap[config.areanames[i]] = i;
        areanamemaprev[i] = config.areanames[i];
        cout << i << "=" << config.areanames[i] << endl;
      }
    } else {
      for (int i = 0; i < ir.nareas; i++) {
        std::ostringstream osstream;
        osstream << i;
        std::string string_x = osstream.str();
        config.areanamemap[string_x] = i;
        areanamemaprev[i] = string_x;
      }
    }
    /*
     * need to figure out how to work with multiple trees best
     */
    if (config.periods.size() < 1) {
      config.periods.push_back(10000);
    }
    auto rm = std::make_shared<RateModel>(ir.nareas, true, config.periods,
                                          config.sparse);
    if (config.numthreads != 0) {
      rm->set_nthreads(config.numthreads);
      cout << "Setting the number of threads: " << config.numthreads << endl;
    }
    rm->setup_Dmask();

    /*
     * if there is a ratematrixfile then it will be processed
     */
    if (!config.ratematrixfile.empty()) {
      parse_ratematrix_file(config, ir, rm);
    }
    /*
      need to add check to make sure that the tips are included in possible
      distributions
    */
    if (config.includedists.size() > 0 || config.excludedists.size() > 0 ||
        config.maxareas >= 2) {
      if (config.excludedists.size() > 0) {
        rm->setup_dists(config.excludedists, false);
      } else {
        if (config.maxareas >= 2) {
          config.includedists =
              generate_dists_from_num_max_areas(ir.nareas, config.maxareas);
        }
        rm->setup_dists(config.includedists, true);
      }
    } else {
      rm->setup_dists();
    }
    rm->setup_D(0.01);
    rm->setup_E(0.01);
    rm->setup_Q();

    /*
     * start calculating on all trees
     */
    for (unsigned int i = 0; i < intrees.size(); i++) {
      handle_tree(intrees[i], rm, data, areanamemaprev, config);
    }
    cout << "expm count: " << rm->get_expm_count() << endl;
  }
  auto end_time = chrono::high_resolution_clock::now();
  chrono::duration<double> duration = end_time - start_time;
  cout << "Analysis took: " << duration.count() << "s" << std::endl;
  return 0;
}
