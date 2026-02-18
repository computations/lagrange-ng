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

TEST(CSVTest, multipleRows) {
  std::string csv_string{"a,b,c\n1,2,3\n4,5,6\n7,8,9\n"};
  CSVReader reader(std::istringstream{csv_string});

  EXPECT_EQ(reader.header().size(), 3);

  auto it = reader.begin();
  auto row1 = *it;
  EXPECT_EQ(row1.get<uint64_t>("a"), 1);
  EXPECT_EQ(row1.get<uint64_t>("b"), 2);
  EXPECT_EQ(row1.get<uint64_t>("c"), 3);

  ++it;
  auto row2 = *it;
  EXPECT_EQ(row2.get<uint64_t>("a"), 4);
  EXPECT_EQ(row2.get<uint64_t>("b"), 5);
  EXPECT_EQ(row2.get<uint64_t>("c"), 6);

  ++it;
  auto row3 = *it;
  EXPECT_EQ(row3.get<uint64_t>("a"), 7);
  EXPECT_EQ(row3.get<uint64_t>("b"), 8);
  EXPECT_EQ(row3.get<uint64_t>("c"), 9);
}

TEST(CSVTest, doubleValues) {
  std::string csv_string{"x,y\n1.5,2.7\n3.14,2.718\n"};
  CSVReader reader(std::istringstream{csv_string});

  auto row = *reader.begin();
  EXPECT_DOUBLE_EQ(row.get<double>("x"), 1.5);
  EXPECT_DOUBLE_EQ(row.get<double>("y"), 2.7);

  auto row2 = *(++reader.begin());
  EXPECT_DOUBLE_EQ(row2.get<double>("x"), 3.14);
  EXPECT_DOUBLE_EQ(row2.get<double>("y"), 2.718);
}

TEST(CSVTest, longValues) {
  std::string csv_string{"val\n1234567890123\n"};
  CSVReader reader(std::istringstream{csv_string});

  auto row = *reader.begin();
  EXPECT_EQ(row.get<long>("val"), 1234567890123L);
}

TEST(CSVTest, booleanValues) {
  std::string csv_string{"flag\n1\n0\n5\n"};
  CSVReader reader(std::istringstream{csv_string});

  auto it = reader.begin();
  auto row1 = *it;
  EXPECT_TRUE(row1.get<bool>("flag"));

  ++it;
  auto row2 = *it;
  EXPECT_FALSE(row2.get<bool>("flag"));

  ++it;
  auto row3 = *it;
  EXPECT_TRUE(row3.get<bool>("flag"));
}

TEST(CSVTest, whitespaceTabs) {
  std::string csv_string{"a\t,\tb\n\t1\t,\t2\t\n"};
  CSVReader reader(std::istringstream{csv_string});

  EXPECT_EQ(reader.header()[0], "a");
  EXPECT_EQ(reader.header()[1], "b");

  auto row = *reader.begin();
  EXPECT_EQ(row.get<std::string>("a"), "1");
  EXPECT_EQ(row.get<std::string>("b"), "2");
}

TEST(CSVTest, multipleSpaces) {
  std::string csv_string{"a   ,   b\n  1  ,  2  \n"};
  CSVReader reader(std::istringstream{csv_string});

  EXPECT_EQ(reader.header()[0], "a");
  EXPECT_EQ(reader.header()[1], "b");

  auto row = *reader.begin();
  EXPECT_EQ(row.get<std::string>("a"), "1");
  EXPECT_EQ(row.get<std::string>("b"), "2");
}

TEST(CSVTest, missingKey) {
  std::string csv_string{"a,b\n1,2\n"};
  CSVReader reader(std::istringstream{csv_string});

  auto row = *reader.begin();
  EXPECT_EQ(row.get<std::string>("c"), "");
}

TEST(CSVTest, fileNotFound) {
  EXPECT_THROW(CSVReader("/nonexistent/path/to/file.csv"), CSVOpenError);
}

TEST(CSVTest, numericIndexAccess) {
  std::string csv_string{"a,b,c\n1,2,3\n"};
  CSVReader reader(std::istringstream{csv_string});

  auto row = *reader.begin();

  auto [h0, v0] = row.get<std::string>(0);
  EXPECT_EQ(h0, "a");
  EXPECT_EQ(v0, "1");

  auto [h1, v1] = row.get<std::string>(1);
  EXPECT_EQ(h1, "b");
  EXPECT_EQ(v1, "2");

  auto [h2, v2] = row.get<std::string>(2);
  EXPECT_EQ(h2, "c");
  EXPECT_EQ(v2, "3");
}
