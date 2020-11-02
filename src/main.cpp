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

int main(int argc, char *argv[]) {
  auto start_time = chrono::high_resolution_clock::now();
  if (argc != 2) {
    cout << "you need more arguments." << endl;
    cout << "usage: lagrange configfile" << endl;
    exit(0);
  } else {
    string treefile;
    string datafile;
    string ratematrixfile;
    string logfile;
    int maxareas = 0;
    vector<double> periods;
    unordered_map<string, vector<string>> mrcas;
    unordered_map<string, lagrange_dist_t> fixnodewithmrca;
    vector<lagrange_dist_t> excludedists;
    vector<lagrange_dist_t> includedists;
    vector<string> areanames;
    unordered_map<string, int> areanamemap;
    unordered_map<int, string> areanamemaprev;
    vector<string> ancstates;  // string should be the mrca or if it is just
                               // _all_
    // then everything will be computed
    vector<string> areacolors;
    vector<string> fossilmrca;
    vector<string> fossiltype;
    vector<string> fossilarea;
    vector<double> fossilage;  // 0's for N type and # for B type

    bool marginal = true;  // false means joint
    bool splits = false;
    bool states = false;
    int numthreads = 0;
    bool sparse = false;

    double dispersal = 0.1;
    double extinction = 0.1;
    bool estimate = true;

    /*
     * for stochastic mapping
     */
    vector<vector<lagrange_dist_t>> stochastic_number_from_tos;
    vector<lagrange_dist_t> stochastic_time_dists;

    // estimating the dispersal mask
    bool estimate_dispersal_mask = false;

    BioGeoTreeTools tt;

    /*************
     * read the configuration file
     **************/
    ifstream ifs(argv[1]);
    string line;
    while (getline(ifs, line)) {
      if (line.size() > 0) {
        if (line[0] != '#') {
          vector<string> tokens;
          string del("=");
          tokens.clear();
          Tokenize(line, tokens, del);
          for (unsigned int j = 0; j < tokens.size(); j++) {
            TrimSpaces(tokens[j]);
          }
          if (!strcmp(tokens[0].c_str(), "treefile")) {
            treefile = tokens[1];
          } else if (!strcmp(tokens[0].c_str(), "datafile")) {
            datafile = tokens[1];
          } else if (!strcmp(tokens[0].c_str(), "ratematrix")) {
            ratematrixfile = tokens[1];
            if (ratematrixfile == "d" || ratematrixfile == "D") {
              ratematrixfile = "";
            }
          } else if (!strcmp(tokens[0].c_str(), "areanames")) {
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
            areanames = searchtokens;
          } else if (!strcmp(tokens[0].c_str(), "fixnode")) {
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
            vector<int> dist;
            for (unsigned int j = 0; j < searchtokens[1].size(); j++) {
              char c = (searchtokens[1].c_str())[j];
              dist.push_back(atoi(&c));
            }
            fixnodewithmrca[searchtokens[0]] =
                convert_vector_to_lagrange_dist(dist);
          } else if (!strcmp(tokens[0].c_str(), "excludedists")) {
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              vector<int> dist;
              for (unsigned int k = 0; k < searchtokens[j].size(); k++) {
                char c = (searchtokens[j].c_str())[k];
                dist.push_back(atoi(&c));
              }
              excludedists.push_back(convert_vector_to_lagrange_dist(dist));
            }
          } else if (!strcmp(tokens[0].c_str(), "includedists")) {
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
            if (searchtokens[0].size() == 1) {
              maxareas = atoi(searchtokens[0].c_str());
            } else {
              for (unsigned int j = 0; j < searchtokens.size(); j++) {
                vector<int> dist;
                for (unsigned int k = 0; k < searchtokens[j].size(); k++) {
                  char c = (searchtokens[j].c_str())[k];
                  dist.push_back(atoi(&c));
                }
                includedists.push_back(convert_vector_to_lagrange_dist(dist));
              }
            }
          } else if (!strcmp(tokens[0].c_str(), "areacolors")) {
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
            areacolors = searchtokens;
          } else if (!strcmp(tokens[0].c_str(), "periods")) {
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
              periods.push_back(atof(searchtokens[j].c_str()));
            }
          } else if (!strcmp(tokens[0].c_str(), "mrca")) {
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
            vector<string> mrc;
            for (unsigned int j = 1; j < searchtokens.size(); j++) {
              mrc.push_back(searchtokens[j]);
            }
            mrcas[searchtokens[0]] = mrc;
          } else if (!strcmp(tokens[0].c_str(), "ancstate")) {
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
            ancstates.push_back(searchtokens[0]);
          } else if (!strcmp(tokens[0].c_str(), "fossil")) {
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
            fossiltype.push_back(searchtokens[0]);
            fossilmrca.push_back(searchtokens[1]);
            fossilarea.push_back(searchtokens[2]);
            if (searchtokens.size() > 3) {
              fossilage.push_back(atof(searchtokens[3].c_str()));
            } else {
              fossilage.push_back(0.0);
            }
          } else if (!strcmp(tokens[0].c_str(), "calctype")) {
            string calctype = tokens[1];
            if (calctype.compare("m") != 0 && calctype.compare("M") != 0) {
              marginal = false;
            }
          } else if (!strcmp(tokens[0].c_str(), "report")) {
            if (tokens[1].compare("split") != 0) {
              splits = false;
            }
          } else if (!strcmp(tokens[0].c_str(), "sparse")) {
            sparse = true;
          } else if (!strcmp(tokens[0].c_str(), "splits")) {
            splits = true;
          } else if (!strcmp(tokens[0].c_str(), "states")) {
            states = true;
          } else if (!strcmp(tokens[0].c_str(), "estimate_dispersal_mask")) {
            estimate_dispersal_mask = true;
          } else if (!strcmp(tokens[0].c_str(), "numthreads")) {
            numthreads = atoi(tokens[1].c_str());
          } else if (!strcmp(tokens[0].c_str(), "stochastic_time")) {
            states = true;  // requires ancestral states
            if (ancstates.size() > 0)
              ancstates[0] = "_all_";
            else
              ancstates.push_back("_all_");
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              vector<int> dist;
              for (unsigned int k = 0; k < searchtokens[j].size(); k++) {
                char c = (searchtokens[j].c_str())[k];
                dist.push_back(atoi(&c));
              }
              stochastic_time_dists.push_back(
                  convert_vector_to_lagrange_dist(dist));
            }
          } else if (!strcmp(tokens[0].c_str(), "stochastic_number")) {
            states = true;  // requires ancestral states
            if (ancstates.size() > 0)
              ancstates[0] = "_all_";
            else
              ancstates.push_back("_all_");
            vector<string> searchtokens;
            Tokenize(tokens[1], searchtokens, ",     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
              TrimSpaces(searchtokens[j]);
            }
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
              stochastic_number_from_tos.push_back(dists);
            }
          } else if (!strcmp(tokens[0].c_str(), "dispersal")) {
            dispersal = atof(tokens[1].c_str());
            cout << "setting dispersal: " << dispersal << endl;
            estimate = false;
          } else if (!strcmp(tokens[0].c_str(), "extinction")) {
            extinction = atof(tokens[1].c_str());
            cout << "setting extinction: " << extinction << endl;
            estimate = false;
          }
        }
      }
    }
    ifs.close();
    /*****************
     * finish reading the configuration file
     *****************/
    /*
     * after reading the input file
     */
    InputReader ir;
    cout << "reading tree..." << endl;
    vector<std::shared_ptr<Tree>> intrees;
    ir.readMultipleTreeFile(treefile, intrees);
    cout << "reading data..." << endl;
    unordered_map<string, lagrange_dist_t> data =
        ir.readStandardInputData(datafile);
    cout << "checking data..." << endl;
    ir.checkData(data, intrees);

    /*
     * read area names
     */
    if (areanames.size() > 0) {
      cout << "reading area names" << endl;
      for (unsigned int i = 0; i < areanames.size(); i++) {
        areanamemap[areanames[i]] = i;
        areanamemaprev[i] = areanames[i];
        cout << i << "=" << areanames[i] << endl;
      }
    } else {
      for (int i = 0; i < ir.nareas; i++) {
        std::ostringstream osstream;
        osstream << i;
        std::string string_x = osstream.str();
        areanamemap[string_x] = i;
        areanamemaprev[i] = string_x;
      }
    }
    /*
     * need to figure out how to work with multiple trees best
     */
    if (periods.size() < 1) {
      periods.push_back(10000);
    }
    auto rm = std::make_shared<RateModel>(ir.nareas, true, periods, sparse);
    if (numthreads != 0) {
      rm->set_nthreads(numthreads);
      cout << "Setting the number of threads: " << numthreads << endl;
    }
    rm->setup_Dmask();

    /*
     * if there is a ratematrixfile then it will be processed
     */
    if (ratematrixfile != "" && ratematrixfile.size() > 0) {
      cout << "Reading rate matrix file" << endl;
      vector<vector<vector<double>>> dmconfig = processRateMatrixConfigFile(
          ratematrixfile, ir.nareas, periods.size());
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
      for (unsigned int i = 0; i < dmconfig.size(); i++) {
        for (unsigned int j = 0; j < dmconfig[i].size(); j++) {
          for (unsigned int k = 0; k < dmconfig[i][j].size(); k++) {
            rm->set_Dmask_cell(i, j, k, dmconfig[i][j][k], false);
          }
        }
      }
    }
    /*
      need to add check to make sure that the tips are included in possible
      distributions
    */
    if (includedists.size() > 0 || excludedists.size() > 0 || maxareas >= 2) {
      if (excludedists.size() > 0) {
        rm->setup_dists(excludedists, false);
      } else {
        if (maxareas >= 2)
          includedists = generate_dists_from_num_max_areas(ir.nareas, maxareas);
        rm->setup_dists(includedists, true);
      }
    } else {
      rm->setup_dists();
    }
    rm->setup_D(0.01);
    rm->setup_E(0.01);
    rm->setup_Q();

    /*
     * outfile for tree reconstructed states
     */
    ofstream outTreeFile;
    ofstream outTreeKeyFile;

    /*
     * outfile for stochastic expectations
     */
    ofstream outStochTimeFile;
    ofstream outStochNumberFile;

    /*
     * start calculating on all trees
     */
    for (unsigned int i = 0; i < intrees.size(); i++) {
      auto bgt =
          std::make_shared<BioGeoTree>(intrees[i], periods, ir.nareas, data);
      /*
       * record the mrcas
       */
      unordered_map<string, std::shared_ptr<Node>> mrcanodeint;
      unordered_map<string, vector<string>>::iterator it;
      for (it = mrcas.begin(); it != mrcas.end(); it++) {
        // records node by number, should maybe just point to node
        mrcanodeint[(*it).first] = intrees[i]->getMRCA((*it).second);
        // tt.getLastCommonAncestor(*intrees[i],nodeIds);
        cout << "Reading mrca: " << (*it).first << " = ";
        for (unsigned int k = 0; k < (*it).second.size(); k++) {
          cout << (*it).second[k] << " ";
        }
        cout << endl;
      }

      /*
       * set fixed nodes
       */
      for (auto fnit = fixnodewithmrca.begin(); fnit != fixnodewithmrca.end();
           fnit++) {
        lagrange_dist_t dista = (*fnit).second;
        for (unsigned int k = 0; k < rm->getDistsSize(); k++) {
          if (dista != rm->getDists()[k]) {
            bgt->set_excluded_dist(rm->getDists().at(k),
                                   mrcanodeint[(*fnit).first]);
          }
        }
        cout << "fixing " << (*fnit).first << " = ";
        print_lagrange_dist((*fnit).second, rm->get_num_areas());
      }

      cout << "setting default model..." << endl;
      bgt->set_default_model(rm);
      cout << "setting up tips..." << endl;
      bgt->set_tip_conditionals(data);

      /*
       * setting up fossils
       */
      for (unsigned int k = 0; k < fossiltype.size(); k++) {
        if (fossiltype[k] == "n" || fossiltype[k] == "N") {
          bgt->setFossilatNodeByMRCA_id(mrcanodeint[fossilmrca[k]],
                                        areanamemap[fossilarea[k]]);
          cout << "Setting node fossil at mrca: " << fossilmrca[k]
               << " at area: " << fossilarea[k] << endl;
        } else {
          bgt->setFossilatBranchByMRCA_id(mrcanodeint[fossilmrca[k]],
                                          areanamemap[fossilarea[k]],
                                          fossilage[k]);
          cout << "Setting branch fossil at mrca: " << fossilmrca[k]
               << " at area: " << fossilarea[k] << " at age: " << fossilage[k]
               << endl;
        }
      }

      /*
       * initial likelihood calculation
       */
      cout << "starting likelihood calculations" << endl;
      cout << "initial -ln likelihood: "
           << double(bgt->eval_likelihood(marginal)) << endl;

      /*
       * optimize likelihood
       */
      Superdouble nlnlike = 0;
      if (estimate == true) {
        if (estimate_dispersal_mask == false) {
          cout << "Optimizing (simplex) -ln likelihood." << endl;
          OptimizeBioGeo opt(bgt, rm, marginal);
          vector<double> disext = opt.optimize_global_dispersal_extinction();
          cout << "dis: " << disext[0] << " ext: " << disext[1] << endl;
          rm->setup_D(disext[0]);
          rm->setup_E(disext[1]);
          rm->setup_Q();
          bgt->update_default_model(rm);
          bgt->set_store_p_matrices(true);
          cout << "final -ln likelihood: "
               << double(bgt->eval_likelihood(marginal)) << endl;
          bgt->set_store_p_matrices(false);
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
              vector<vector<vector<double>>>(periods.size(), rows);
          int count = 2;
          for (unsigned int ii = 0; ii < D_mask.size(); ii++) {
            for (unsigned int jj = 0; jj < D_mask[ii].size(); jj++) {
              D_mask[ii][jj][jj] = 0.0;
              for (unsigned int kk = 0; kk < D_mask[ii][jj].size(); kk++) {
                if (kk != jj) {
                  D_mask[ii][jj][kk] = disextrm[count];
                  count += 1;
                }
              }
            }
          }
          cout << "D_mask" << endl;
          for (unsigned int ii = 0; ii < D_mask.size(); ii++) {
            cout << periods.at(ii) << endl;
            cout << "\t";
            for (unsigned int j = 0; j < D_mask[ii].size(); j++) {
              cout << areanames[j] << "\t";
            }
            cout << endl;
            for (unsigned int j = 0; j < D_mask[ii].size(); j++) {
              cout << areanames[j] << "\t";
              for (unsigned int k = 0; k < D_mask[ii][j].size(); k++) {
                cout << D_mask[ii][j][k] << "\t";
              }
              cout << endl;
            }
            cout << endl;
          }
          rm->setup_D_provided(disextrm[0], D_mask);
          rm->setup_E(disextrm[1]);
          rm->setup_Q();
          bgt->update_default_model(rm);
          bgt->set_store_p_matrices(true);
          nlnlike = bgt->eval_likelihood(marginal);
          cout << "final -ln likelihood: " << double(nlnlike) << endl;
          bgt->set_store_p_matrices(false);
        }
      } else {
        rm->setup_D(dispersal);
        rm->setup_E(extinction);
        rm->setup_Q();
        bgt->update_default_model(rm);
        bgt->set_store_p_matrices(true);
        cout << "final -ln likelihood: "
             << double(bgt->eval_likelihood(marginal)) << endl;
        bgt->set_store_p_matrices(false);
      }

      /*
       * ancestral splits calculation
       */
      if (ancstates.size() > 0) {
        bgt->set_use_stored_matrices(true);

        bgt->prepare_ancstate_reverse();
        Superdouble totlike = 0;  // calculate_vector_double_sum(rast) , should
                                  // be the same for every node

        if (ancstates[0] == "_all_" || ancstates[0] == "_ALL_") {
          for (unsigned int j = 0; j < intrees[i]->getInternalNodeCount();
               j++) {
            if (splits) {
              cout << "Ancestral splits for:\t"
                   << intrees[i]->getInternalNode(j)->getNumber() << endl;
              unordered_map<lagrange_dist_t, vector<AncSplit>> ras =
                  bgt->calculate_ancsplit_reverse(
                      intrees[i]->getInternalNode(j));
              tt.summarizeSplits(intrees[i]->getInternalNode(j), ras,
                                 areanamemaprev, rm->get_int_dists_map(),
                                 rm->get_num_areas());
              cout << endl;
            }
            if (states) {
              cout << "Ancestral states for:\t"
                   << intrees[i]->getInternalNode(j)->getNumber() << endl;
              vector<Superdouble> rast = bgt->calculate_ancstate_reverse(
                  intrees[i]->getInternalNode(j));
              totlike = calculate_vector_Superdouble_sum(rast);
              tt.summarizeAncState(intrees[i]->getInternalNode(j), rast,
                                   areanamemaprev, rm->get_int_dists_map(),
                                   rm->get_num_areas());
              cout << endl;
            }
            // exit(0);
          }
          /*
           * key file output
           */
          outTreeKeyFile.open((treefile + ".bgkey.tre").c_str(), ios::app);
          // need to output numbers
          outTreeKeyFile << intrees[i]->getRoot()->getNewick(
                                true,
                                [](const Node &n) -> string {
                                  return std::to_string(n.getNumber());
                                })
                         << ";" << endl;
          outTreeKeyFile.close();
        } else {
          for (unsigned int j = 0; j < ancstates.size(); j++) {
            if (splits) {
              cout << "Ancestral splits for: " << ancstates[j] << endl;
              unordered_map<lagrange_dist_t, vector<AncSplit>> ras =
                  bgt->calculate_ancsplit_reverse(mrcanodeint[ancstates[j]]);
              tt.summarizeSplits(mrcanodeint[ancstates[j]], ras, areanamemaprev,
                                 rm->get_int_dists_map(), rm->get_num_areas());
            }
            if (states) {
              cout << "Ancestral states for: " << ancstates[j] << endl;
              vector<Superdouble> rast =
                  bgt->calculate_ancstate_reverse(mrcanodeint[ancstates[j]]);
              tt.summarizeAncState(mrcanodeint[ancstates[j]], rast,
                                   areanamemaprev, rm->get_int_dists_map(),
                                   rm->get_num_areas());
            }
          }
        }
        if (splits) {
          outTreeFile.open((treefile + ".bgsplits.tre").c_str(), ios::app);
          // need to output object "split"
          outTreeFile << intrees[i]->getRoot()->getNewick(
                             true,
                             [](const Node &n) -> string {
                               return n.getSplitString();
                             })
                      << ";" << endl;
          outTreeFile.close();
        }
        if (states) {
          outTreeFile.open((treefile + ".bgstates.tre").c_str(), ios::app);
          // need to output object "state"
          outTreeFile << intrees[i]->getRoot()->getNewick(
                             true,
                             [](const Node &n) -> string {
                               return n.getStateString();
                             })
                      << ";" << endl;
          outTreeFile.close();
        }
        /*
         * stochastic mapping calculations
         * REQUIRES that ancestral calculation be done
         */
        if (stochastic_time_dists.size() > 0) {
          cout << "calculating stochastic mapping time spent" << endl;
          for (unsigned int k = 0; k < stochastic_time_dists.size(); k++) {
            cout << tt.get_string_from_dist_int(
                        rm->get_dists_int_map().at(stochastic_time_dists[k]),
                        areanamemaprev, rm->get_int_dists_map(),
                        rm->get_num_areas())
                 << endl;
            bgt->prepare_stochmap_reverse_all_nodes(
                rm->get_dists_int_map().at(stochastic_time_dists[k]),
                rm->get_dists_int_map().at(stochastic_time_dists[k]));
            bgt->prepare_ancstate_reverse();
            outStochTimeFile.open((treefile + ".bgstochtime.tre").c_str(),
                                  ios::app);
            for (unsigned int j = 0; j < intrees[i]->getNodeCount(); j++) {
              if (intrees[i]->getNode(j) != intrees[i]->getRoot()) {
                vector<Superdouble> rsm = bgt->calculate_reverse_stochmap(
                    intrees[i]->getNode(j), true);
                double stres = calculate_vector_Superdouble_sum(rsm) / totlike;
                intrees[i]->getNode(j)->setStochString(std::to_string(stres));
              }
            }
            // need to output object "stoch"
            outStochTimeFile
                << intrees[i]->getRoot()->getNewickLambda(
                       [](const Node &n) -> string {
                         return n.getStochString();
                       })
                << tt.get_string_from_dist_int(
                       rm->get_dists_int_map().at(stochastic_time_dists[k]),
                       areanamemaprev, rm->get_int_dists_map(),
                       rm->get_num_areas())
                << ";" << endl;
            outStochTimeFile.close();
          }
        }
        if (stochastic_number_from_tos.size() > 0) {
          cout << "calculating stochastic mapping number of transitions"
               << endl;
          for (unsigned int k = 0; k < stochastic_number_from_tos.size(); k++) {
            cout << tt.get_string_from_dist_int(
                        rm->get_dists_int_map().at(
                            stochastic_number_from_tos[k][0]),
                        areanamemaprev, rm->get_int_dists_map(),
                        rm->get_num_areas())
                 << " -> "
                 << tt.get_string_from_dist_int(
                        rm->get_dists_int_map().at(
                            stochastic_number_from_tos[k][1]),
                        areanamemaprev, rm->get_int_dists_map(),
                        rm->get_num_areas())
                 << endl;
            bgt->prepare_stochmap_reverse_all_nodes(
                rm->get_dists_int_map().at(stochastic_number_from_tos[k][0]),
                rm->get_dists_int_map().at(stochastic_number_from_tos[k][1]));
            bgt->prepare_ancstate_reverse();
            outStochTimeFile.open((treefile + ".bgstochnumber.tre").c_str(),
                                  ios::app);
            for (unsigned int j = 0; j < intrees[i]->getNodeCount(); j++) {
              if (intrees[i]->getNode(j) != intrees[i]->getRoot()) {
                vector<Superdouble> rsm = bgt->calculate_reverse_stochmap(
                    intrees[i]->getNode(j), false);
                // cout << calculate_vector_double_sum(rsm) / totlike << endl;
                double stres = calculate_vector_Superdouble_sum(rsm) / totlike;
                intrees[i]->getNode(j)->setStochString(std::to_string(stres));
              }
            }
            // need to output object "stoch"
            outStochTimeFile << intrees[i]->getRoot()->getNewickLambda(
                                    [](const Node &n) -> string {
                                      return n.getStochString();
                                    })
                             << ";" << endl;
            outStochTimeFile.close();
          }
        }
        /*
         * end stochastic mapping
         */
      }
    }
    cout << "expm count: " << rm->get_expm_count() << endl;
  }
  auto end_time = chrono::high_resolution_clock::now();
  chrono::duration<double> duration = end_time - start_time;
  cout << "Analysis took: " << duration.count() << "s" << std::endl;
  return 0;
}
