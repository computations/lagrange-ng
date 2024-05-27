#include "ConfigFile.hpp"

#include "Utils.hpp"
#include "environment.hpp"
#include "gtest/gtest.h"

using namespace lagrange;

class ConfigFileTest : public ::testing::Test {
 protected:
  void SetUp() override { _test1 = env->get_test1_config(); }

  std::ifstream _test1;
};

TEST_F(ConfigFileTest, basic) {
  auto config = parse_config_file(_test1);
  EXPECT_EQ(config.tree_filename, "foo.nwk");
  EXPECT_EQ(config.data_filename, "foo.phy");

  std::vector<std::string> expected_area_names = {"RA", "RB", "RC"};
  ASSERT_EQ(config.area_names.size(), expected_area_names.size());

  for (size_t i = 0; i < expected_area_names.size(); ++i) {
    EXPECT_EQ(expected_area_names[i], config.area_names[i]);
  }

  EXPECT_TRUE(config.all_states);

  EXPECT_TRUE(config.workers.hasValue());
  EXPECT_EQ(config.workers, 1);

  EXPECT_EQ(config.mrcas.size(), 1);
  EXPECT_TRUE(config.mrcas.at("bar"));
  EXPECT_THROW(config.mrcas.at("foo"), std::out_of_range);

  EXPECT_EQ(config.fossils.size(), 1);

  auto fossil = config.fossils.front();
  EXPECT_EQ(fossil.mrca_name, "bar");
  EXPECT_EQ(fossil.type, FossilType::INCLUDE);
  EXPECT_EQ(fossil.area, convert_dist_binary_string_to_dist("100"));

  EXPECT_EQ(config.prefix, "foo");
}
