#include "CSV.hpp"

#include <sstream>

#include "gtest/gtest.h"

TEST(CSVTest, simple) {
  std::string csv_string{" a , b \n 1 , 2\n"};
  CSVReader reader(std::istringstream{csv_string});

  EXPECT_EQ(reader.header()[0], "a");
  EXPECT_EQ(reader.header()[1], "b");

  auto row = *reader.begin();

  EXPECT_EQ(row.get<std::string>("a"), "1");
  EXPECT_EQ(row.get<std::string>("b"), "2");

  EXPECT_EQ(row.get<uint64_t>("a"), 1);
  EXPECT_EQ(row.get<uint64_t>("b"), 2);
}

TEST(CSVTest, errors) {
  std::string csv_string{" a , b \n c , d\n"};
  CSVReader reader(std::istringstream{csv_string});

  EXPECT_EQ(reader.header()[0], "a");
  EXPECT_EQ(reader.header()[1], "b");

  auto row = *reader.begin();

  EXPECT_EQ(row.get<std::string>("a"), "c");
  EXPECT_EQ(row.get<std::string>("b"), "d");

  EXPECT_ANY_THROW(row.get<uint64_t>("a"));
  EXPECT_ANY_THROW(row.get<uint64_t>("b"));
}
