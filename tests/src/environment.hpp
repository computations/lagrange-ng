#ifndef LAGRANGE_ENV_HPP
#define LAGRANGE_ENV_HPP

#include <fstream>
#include <iostream>

#include "Common.h"
#include "Utils.h"
#include "gtest/gtest.h"
#include "test_data_paths.h"

class LagrangeEnvironment : public ::testing::Environment {
 public:
  LagrangeEnvironment() {}

  std::ifstream get_datafile() { return std::ifstream(RANDOM_TREE_FILE); }

  std::ifstream get_pathological_data() {
    return std::ifstream(PATHOLOGICAL_TREE_FILE);
  }
};

extern LagrangeEnvironment *env;

#endif
