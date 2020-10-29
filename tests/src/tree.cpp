#include <tree_reader.h>

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>

#include "Common.h"
#include "Workspace.h"
#include "environment.hpp"
#include "gtest/gtest.h"
#include "node.h"
#include "tree.h"

class TreeTest : public ::testing::Test {
 protected:
  void SetUp() override {
    _basic_tree_newick = "((a:1.0,b:1.0)1.0, c:2.0);";
    _basic_tree_dist_data = {{"a", 0b01}, {"b", 0b11}, {"c", 0b01}};
    _basic_ws = std::make_shared<Workspace>(/*taxa=*/3,
                                            /*regions=*/2);
  }

  std::shared_ptr<Tree> parse_tree(const std::string& newick) {
    return TreeReader().readTree(newick);
  }

  std::string _basic_tree_newick;
  std::unordered_map<std::string, lagrange_dist_t> _basic_tree_dist_data;
  std::shared_ptr<Workspace> _basic_ws;
};

TEST_F(TreeTest, simple0) { auto t = parse_tree(_basic_tree_newick); }

TEST_F(TreeTest, generateOperationsSimple0) {
  auto t = parse_tree(_basic_tree_newick);
  auto ops =
      t->generateOperations(*_basic_ws, _basic_tree_dist_data, true, false, false);

  EXPECT_EQ(ops._forward_ops.size(), 2);
  EXPECT_EQ(ops._backwards_ops.size(), 0);
  EXPECT_EQ(ops._lh.size(), 1);
  EXPECT_EQ(ops._state.size(), 0);
  EXPECT_EQ(ops._split.size(), 0);
}

TEST_F(TreeTest, generateOperationsSimple1) {
  auto t = parse_tree(_basic_tree_newick);
  auto ops =
      t->generateOperations(*_basic_ws, _basic_tree_dist_data, true, true, false);

  EXPECT_EQ(ops._forward_ops.size(), 2);
  EXPECT_EQ(ops._backwards_ops.size(), 1);
  EXPECT_EQ(ops._lh.size(), 1);
  EXPECT_EQ(ops._state.size(), 2);
  EXPECT_EQ(ops._split.size(), 0);
}
