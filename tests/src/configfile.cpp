#include "ConfigFile.hpp"

#include <ranges>
#include <sstream>

#include "Utils.hpp"
#include "environment.hpp"
#include "gtest/gtest.h"

using namespace lagrange;

class ConfigFileTest : public ::testing::Test {
 protected:
  void SetUp() override {
    _test1 = env->get_test1_config();
    _test2 = env->get_test2_config();
  }

  std::ifstream _test1;
  std::ifstream _test2;
};

TEST_F(ConfigFileTest, basic) {
  ConfigFile config{_test1};
  EXPECT_EQ(config.tree_filename(), "foo.nwk");
  EXPECT_EQ(config.data_filename(), "foo.phy");

  std::vector<std::string> expected_area_names = {"RA", "RB", "RC"};
  ASSERT_EQ(config.area_names().size(), expected_area_names.size());

#ifdef __cpp_lib_ranges_enumerate
  for (auto const [index, it] : expected_area_names | std::views::enumerate) {
    EXPECT_EQ(it, config.area_names()[index]);
  }
#else
  for (size_t i = 0; i < expected_area_names.size(); ++i) {
    EXPECT_EQ(expected_area_names[i], config.area_names()[i]);
  }
#endif

  EXPECT_TRUE(config.compute_all_states());

  EXPECT_EQ(config.workers(), 1);

  EXPECT_EQ(config.mrcas().size(), 1);
  EXPECT_TRUE(config.mrca("bar"));
  EXPECT_THROW(config.mrca("foo"), std::out_of_range);

  EXPECT_EQ(config.fossils().size(), 1);

  auto fossil = config.fossils().front();
  EXPECT_EQ(fossil.mrca_name, "bar");
  EXPECT_EQ(fossil.type, FossilType::INCLUDE);
  EXPECT_EQ(fossil.area, convert_dist_binary_string_to_dist("100"));

  EXPECT_EQ(config.prefix(), "foo");
}

TEST_F(ConfigFileTest, quoted_values) {
  ConfigFile config{_test2};
  EXPECT_EQ(config.tree_filename(), "bar foo.nwk");
  EXPECT_EQ(config.data_filename(), "bar foo.phy");

  std::vector<std::string> expected_area_names = {"R A", "R' B", "R C"};
  ASSERT_EQ(config.area_names().size(), expected_area_names.size());

  for (size_t i = 0; i < expected_area_names.size(); ++i) {
    EXPECT_EQ(expected_area_names[i], config.area_names()[i]);
  }

  EXPECT_TRUE(config.compute_all_states());

  EXPECT_EQ(config.workers(), 1);

  EXPECT_EQ(config.mrcas().size(), 1);
  EXPECT_TRUE(config.mrca("baz.bar"));
  EXPECT_THROW(config.mrca("foo"), std::out_of_range);

  EXPECT_EQ(config.fossils().size(), 1);

  auto fossil = config.fossils().front();
  EXPECT_EQ(fossil.mrca_name, "baz.bar");
  EXPECT_EQ(fossil.type, FossilType::INCLUDE);
  EXPECT_EQ(fossil.area, convert_dist_binary_string_to_dist("100"));

  EXPECT_EQ(config.prefix(), "foo");
}

TEST_F(ConfigFileTest, lines) {
  auto failure_lines = {
      "treefile = ",
      "treefile",
      "treefile = foo bar",
      "treefile = ../test tree.nkw",
      "datafile = ",
      "datafile",
      "prefix",
      "areanames",
      "areanames = ",
      "workers = asdbs",
      "workers = ",
      "workers ",
      "foo",
      "period exclude = 101",
      "period foo\nperiod foo",
      "period foo bar",
      "period foo\nperiod foo start = 1.0 2.0",
      "period foo\nperiod foo start = bar",
      "period foo\nperiod foo matrix = foo bar",
  };

  for (auto line : failure_lines) {
    std::istringstream iss{line};
    EXPECT_THROW(ConfigFile(iss, true), ConfigFileParsingError);
  }

  auto success_lines = {
      "treefile = test.nwk",
      "datafile = test.phy",
      "datafile = test.fasta",
      "datafile = 'foo bar'",
      "areanames = a b c",
      "areanames = 1 2 3",
      "states",
      "splits",
      "workers = 10",
      "threads-per-worker = 10",
      "dispersion = 1.2",
      "extinction = 1.2",
      "lh-epsilon = 1e-4",
      "maxareas = 3",
      "expm-mode = krylov",
      "mode = optimize",
      "mrca foo = a b c",
      "mrca foo = 'a d' b c",
      "mrca foo = a b c \nfossil include foo = 011",
      "mrca foo = a b c \nfossil exclude foo = 011",
      "mrca foo = a b c \nfossil fixed foo = 011",
      "mrca foo = a b c \nfossil node foo = 011",
      "mrca foo = a b c \nfossil branch foo = 011 0.1",
      "logfile = test.log",
      "logfile = 'test a'",
      "logfile = 'test\" a'",
      "period foo\nperiod foo start = 1.0",
      "period exclude\nperiod exclude exclude = 101",
      "period foo\nperiod foo matrix = 'foo bar'",
  };

  for (auto line : success_lines) {
    std::istringstream iss{line};
    EXPECT_NO_THROW(ConfigFile(iss, true));
  }
}
