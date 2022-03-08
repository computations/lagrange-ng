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
#include <vector>

#include "Node.h"
#include "Tree.h"

class TreeReader {
 public:
  TreeReader() {}
  std::shared_ptr<Tree> readTree(const std::string &tree);
};

#endif /* TREE_READER_H_ */
