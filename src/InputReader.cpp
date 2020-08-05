/*
 * InputReader.cpp
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 */

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
using namespace std;

#include "BioGeoTreeTools.h"
#include "InputReader.h"
#include "RateMatrixUtils.h"
#include "Utils.h"

#include "node.h"
#include "tree.h"
#include "tree_reader.h"

InputReader::InputReader() : nareas(0), nspecies(0) {}

void InputReader::readMultipleTreeFile(string filename,
                                       vector<std::shared_ptr<Tree>> &ret) {
  TreeReader tr;
  ifstream ifs(filename.c_str());
  string temp;
  int count = 1;
  while (getline(ifs, temp)) {
    if (temp.size() > 1) {
      auto intree = tr.readTree(temp);
      cout << "Tree " << count << " has " << intree->getExternalNodeCount()
           << " leaves." << endl;
      ret.push_back(intree);
      count++;
    }
  }
}

unordered_map<string, vector<int>>
InputReader::readStandardInputData(string filename) {
  ifstream ifs(filename.c_str());
  string temp;
  bool first = false;
  nareas = 0;
  nspecies = 0;
  unordered_map<string, vector<int>> data;
  string line;
  vector<string> tokens;
  string del("\t ");
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
      cout << "Reading species: " << tokens[0] << " ";
      vector<int> speciesdata(nareas, 0);
      for (int i = 0; i < nareas; i++) {
        char spot = tokens[1][i];
        if (spot == '1')
          speciesdata[i] = 1;
        cout << spot - '0';
      }
      cout << endl;
      data[tokens[0]] = speciesdata;
    }
  }
  ifs.close();
  return data;
}

void InputReader::checkData(const unordered_map<string, vector<int>> &data,
                            const vector<std::shared_ptr<Tree>> &trees) {
  vector<string> dataspecies;
  unordered_map<string, vector<int>>::const_iterator itr;
  for (itr = data.begin(); itr != data.end(); ++itr) {
    dataspecies.push_back(itr->first);
  }
  vector<string> treespecies;
  for (int j = 0; j < trees[0]->getExternalNodeCount(); j++) {
    treespecies.push_back(trees[0]->getExternalNode(j)->getName());
    int count = 0;
    for (unsigned int k = 0; k < dataspecies.size(); k++) {
      if (trees[0]->getExternalNode(j)->getName() == dataspecies[k])
        count += 1;
    }
    if (count != 1) {
      cout << "Error: " << trees[0]->getExternalNode(j)->getName() << " found"
           << count << " times in data file." << endl;
      exit(0);
    }
  }
  for (size_t j = 0; j < dataspecies.size(); j++) {
    int count = 0;
    for (size_t k = 0; k < treespecies.size(); k++) {
      if (dataspecies[j] == treespecies[k]) {
        count += 1;
      }
    }
    if (count != 1) {
      cerr << "Error: " << dataspecies[j] << " found " << count
           << " times in tree file." << endl;
      exit(0);
    }
  }
}
