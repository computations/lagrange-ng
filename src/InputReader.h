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
  static auto readMultipleTreeFile(const std::string &filename)
      -> std::vector<std::shared_ptr<Tree>>;
  static void checkData(
      const std::unordered_map<std::string, lagrange_dist_t> &,
      const std::vector<std::shared_ptr<Tree>> &);
  size_t nareas;
  size_t nspecies;
};

#endif /* INPUTREADER_H_ */
