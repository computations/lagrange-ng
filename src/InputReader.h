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

#include "Tree.h"

class InputReader {
 public:
  InputReader();
  void readMultipleTreeFile(std::string filename,
                            std::vector<std::shared_ptr<Tree>> &);
  std::unordered_map<std::string, lagrange_dist_t> readStandardInputData(
      std::string filename);
  void checkData(const std::unordered_map<std::string, lagrange_dist_t> &,
                 const std::vector<std::shared_ptr<Tree>> &);
  int nareas;
  int nspecies;
};

#endif /* INPUTREADER_H_ */
