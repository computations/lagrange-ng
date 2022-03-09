/*
 * InputReader.h
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef INPUTREADER_COPPER_H
#define INPUTREADER_COPPER_H

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Tree.h"

class InputReader {
 public:
  InputReader();
  std::vector<std::shared_ptr<Tree>> readMultipleTreeFile(
      const std::string &filename);
  std::unordered_map<std::string, lagrange_dist_t> readStandardInputData(
      std::string filename, size_t max_areas);
  void checkData(const std::unordered_map<std::string, lagrange_dist_t> &,
                 const std::vector<std::shared_ptr<Tree>> &);
  size_t nareas;
  size_t nspecies;
};

#endif /* INPUTREADER_H_ */
