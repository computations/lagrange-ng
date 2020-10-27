/*
 * InputReader.h
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
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
  unordered_map<string, lagrange_dist_t> readStandardInputData(string filename);
  void checkData(const unordered_map<string, lagrange_dist_t> &,
                 const vector<std::shared_ptr<Tree>> &);
  int nareas;
  int nspecies;
};

#endif /* INPUTREADER_H_ */
