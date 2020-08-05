/*
 * InputReader.h
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 */

#ifndef INPUTREADER_COPPER_H_
#define INPUTREADER_COPPER_H_

#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

#include "tree.h"

class InputReader {
public:
  InputReader();
  void readMultipleTreeFile(string filename, vector<std::shared_ptr<Tree>> &);
  unordered_map<string, vector<int>> readStandardInputData(string filename);
  void checkData(const unordered_map<string, vector<int>> &,
                 const vector<std::shared_ptr<Tree>> &);
  int nareas;
  int nspecies;
};

#endif /* INPUTREADER_H_ */
