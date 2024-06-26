#ifndef LAGRANGE_ENV_HPP
#define LAGRANGE_ENV_HPP

#include <fstream>

#include "gtest/gtest.h"

#define STRING(s) #s
#define STRINGIFY(s) STRING(s)

class LagrangeEnvironment : public ::testing::Environment {
 public:
  LagrangeEnvironment() {}

  std::ifstream get_datafile() {
    return std::ifstream(STRINGIFY(TREEPATH/random_test_trees));
  }

  std::ifstream get_pathological_data() {
    return std::ifstream(STRINGIFY(TREEPATH/pathological_trees));
  }

  std::ifstream get_sloth_tree() {
    return std::ifstream(STRINGIFY(TREEPATH/sloths.tre));
  }

  std::ifstream get_test1_config() {
    return std::ifstream(STRINGIFY(TREEPATH/test1.conf));
  }
  std::ifstream get_test2_config() {
    return std::ifstream(STRINGIFY(TREEPATH/test2.conf));
  }
};

extern LagrangeEnvironment *env;

#endif
