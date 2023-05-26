/*
 * InputReader.cpp
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include "InputReader.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

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
