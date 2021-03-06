/*
 * InputReader.cpp
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include <fstream>
#include <iostream>

#include "InputReader.h"
#include "TreeReader.h"
#include "Utils.h"

InputReader::InputReader() : nareas(0), nspecies(0) {}

void InputReader::readMultipleTreeFile(
    std::string filename, std::vector<std::shared_ptr<Tree>> &ret) {
  TreeReader tr;
  std::ifstream ifs(filename.c_str());
  std::string temp;
  int count = 1;
  while (getline(ifs, temp)) {
    if (temp.size() > 1) {
      auto intree = tr.readTree(temp);
      std::cout << "Tree " << count << " has " << intree->getExternalNodeCount()
                << " leaves." << std::endl;
      ret.push_back(intree);
      count++;
    }
  }
}

std::unordered_map<std::string, lagrange_dist_t>
InputReader::readStandardInputData(std::string filename) {
  std::ifstream ifs(filename.c_str());
  std::string temp;
  bool first = false;
  nareas = 0;
  nspecies = 0;
  std::unordered_map<std::string, lagrange_dist_t> data;
  std::string line;
  std::vector<std::string> tokens;
  std::string del("\t ");
  while (getline(ifs, line)) {
    if (first == false) {
      first = true;
      tokens.clear();
      Tokenize(line, tokens, del);
      for (unsigned int j = 0; j < tokens.size(); j++) {
        TrimSpaces(tokens[j]);
      }
      nspecies = atoi(tokens[0].c_str());
      nareas = atoi(tokens[1].c_str());
    } else {
      tokens.clear();
      Tokenize(line, tokens, del);
      for (unsigned int j = 0; j < tokens.size(); j++) {
        TrimSpaces(tokens[j]);
      }
      std::cout << "Reading species: " << tokens[0] << " ";
      std::vector<int> speciesdata(nareas, 0);
      for (int i = 0; i < nareas; i++) {
        char spot = tokens[1][i];
        if (spot == '1') speciesdata[i] = 1;
        std::cout << spot - '0';
      }
      std::cout << std::endl;
      data[tokens[0]] = convert_vector_to_lagrange_dist(speciesdata);
    }
  }
  ifs.close();
  return data;
}

void InputReader::checkData(
    const std::unordered_map<std::string, lagrange_dist_t> &data,
    const std::vector<std::shared_ptr<Tree>> &trees) {
  std::vector<std::string> dataspecies;
  for (auto itr = data.begin(); itr != data.end(); ++itr) {
    dataspecies.push_back(itr->first);
  }
  std::vector<std::string> treespecies;
  for (unsigned int j = 0; j < trees[0]->getExternalNodeCount(); j++) {
    treespecies.push_back(trees[0]->getExternalNode(j)->getName());
    int count = 0;
    for (unsigned int k = 0; k < dataspecies.size(); k++) {
      if (trees[0]->getExternalNode(j)->getName() == dataspecies[k]) count += 1;
    }
    if (count != 1) {
      std::cout << "Error: " << trees[0]->getExternalNode(j)->getName()
                << " found" << count << " times in data file." << std::endl;
      exit(0);
    }
  }
  for (size_t j = 0; j < dataspecies.size(); j++) {
    int count = 0;
    for (size_t k = 0; k < treespecies.size(); k++) {
      if (dataspecies[j] == treespecies[k]) { count += 1; }
    }
    if (count != 1) {
      std::cerr << "Error: " << dataspecies[j] << " found " << count
                << " times in tree file." << std::endl;
      exit(0);
    }
  }
}
