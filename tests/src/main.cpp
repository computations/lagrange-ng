#include "environment.hpp"
#include "gtest/gtest.h"

LagrangeEnvironment *env = nullptr;

int main(int argc, char **argv) {
  env = new LagrangeEnvironment();
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(env);
  return RUN_ALL_TESTS();
}
