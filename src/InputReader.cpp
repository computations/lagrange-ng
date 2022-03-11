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
#include <memory>
#include <stdexcept>

#include "InputReader.h"
#include "TreeReader.h"
#include "Utils.h"

InputReader::InputReader() : nareas(0), nspecies(0) {}

auto InputReader::readMultipleTreeFile(const std::string &filename)
    -> std::vector<std::shared_ptr<Tree>> {
  std::vector<std::shared_ptr<Tree>> ret;
  std::ifstream ifs(filename.c_str());
  std::string temp;
  int count = 1;
  while (getline(ifs, temp)) {
    if (temp.size() > 1) {
      auto intree = TreeReader::readTree(temp);
      std::cout << "Tree " << count << " has " << intree->getExternalNodeCount()
                << " leaves." << std::endl;
      ret.push_back(intree);
      count++;
    }
  }
  return ret;
}

auto InputReader::readStandardInputData(const std::string &filename,
                                        size_t max_areas)
    -> std::unordered_map<std::string, size_t> {
  std::ifstream ifs(filename.c_str());
  nareas = 0;
  nspecies = 0;
  std::unordered_map<std::string, size_t> data;
  std::string line;
  std::vector<std::string> tokens;
  std::string del("\t ");

  getline(ifs, line);

  Tokenize(line, tokens, del);
  for (auto &token : tokens) { TrimSpaces(token); }

  nspecies = lagrange_parse_size_t(tokens[0]);
  nareas = lagrange_parse_size_t(tokens[1]);

  if (max_areas == 0) { max_areas = nareas; }

  while (getline(ifs, line)) {
    tokens.clear();

    Tokenize(line, tokens, del);

    for (auto &token : tokens) { TrimSpaces(token); }
    std::cout << "Reading species: " << tokens[0] << " ";

    std::vector<int> speciesdata(nareas, 0);

    for (size_t i = 0; i < nareas; i++) {
      char spot = tokens[1][i];
      if (spot == '1') { speciesdata[i] = 1; }
      std::cout << spot - '0';
    }
    std::cout << std::endl;

    lagrange_dist_t dist = convert_vector_to_lagrange_dist(speciesdata);

    try {
      data[tokens[0]] = compute_index_from_dist(dist, max_areas);
    } catch (std::lagrange_util_dist_index_conversion_exception &) {
      throw std::runtime_error(
          std::string(
              "found invalid dist when parsing the dist for species: ") +
          tokens[0]);
    }
  }
  ifs.close();
  return data;
}

void InputReader::checkData(
    const std::unordered_map<std::string, lagrange_dist_t> &data,
    const std::vector<std::shared_ptr<Tree>> &trees) {
  std::vector<std::string> dataspecies;
  dataspecies.reserve(data.size());
  for (const auto &itr : data) { dataspecies.push_back(itr.first); }
  std::vector<std::string> treespecies;
  for (unsigned int j = 0; j < trees[0]->getExternalNodeCount(); j++) {
    treespecies.push_back(trees[0]->getExternalNode(j)->getName());
    int count = 0;
    for (auto &ds : dataspecies) {
      if (trees[0]->getExternalNode(j)->getName() == ds) { count += 1; }
    }
    if (count != 1) {
      std::cout << "Error: " << trees[0]->getExternalNode(j)->getName()
                << " found " << count << " times in data file." << std::endl;
      exit(0);
    }
  }
  for (auto &dataspecie : dataspecies) {
    int count = 0;
    for (auto &treespecie : treespecies) {
      if (dataspecie == treespecie) { count += 1; }
    }
    if (count != 1) {
      std::cerr << "Error: " << dataspecie << " found " << count
                << " times in tree file." << std::endl;
      exit(0);
    }
  }
}
