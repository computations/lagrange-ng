/*
 * tree_reader.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef TREE_READER_H_
#define TREE_READER_H_

#include <string>
#include <vector>

using namespace std;

#include "node.h"
#include "tree.h"

class TreeReader {
public:
  TreeReader() {}
  std::shared_ptr<Tree> readTree(string tree);
};

#endif /* TREE_READER_H_ */
