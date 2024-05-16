/*
 * tree_reader.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef TREE_READER_H
#define TREE_READER_H

#include <string>
 
#include "Tree.hpp"

class TreeReader {
 public:
  TreeReader() = default;
  static auto readTree(const std::string &tree) -> std::shared_ptr<Tree>;
};

#endif /* TREE_READER_H_ */
