#include <stdexcept>
#include <string>

#include "TreeReader.hpp"
#include "environment.hpp"
#include "gtest/gtest.h"

using namespace lagrange;

TEST(TreeReader, simple0) {
  auto tr = TreeReader();
  auto t = tr.readTree("(a:1.0,b:1.0);");
  EXPECT_EQ(t->getTipCount(), 2);
  EXPECT_EQ(t->getNode(1)->getName(), "a");
  EXPECT_EQ(t->getNode(2)->getName(), "b");
  EXPECT_EQ(t->getNode(1)->getBL(), 1.0);
  EXPECT_EQ(t->getNode(2)->getBL(), 1.0);

  EXPECT_EQ(t->getInternalCount(), 1);
  EXPECT_EQ(t->getNode(0)->getName(), "");
}

TEST(TreeReader, simple1) {
  auto tr = TreeReader();
  auto t = tr.readTree("((a,b)ab,(c,d)cd)root;");
  EXPECT_EQ(t->getNodeCount(), 7);
  EXPECT_EQ(t->getTipCount(), 4);
  EXPECT_EQ(t->getInternalCount(), 3);

  EXPECT_EQ(t->getNode(0)->getName(), "root");
  EXPECT_EQ(t->getNode(1)->getName(), "ab");
  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_EQ(t->getNode(4)->getName(), "cd");
  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_EQ(t->getNode(6)->getName(), "d");
}

TEST(TreeReader, simple2) {
  auto tr = TreeReader();
  auto t = tr.readTree("((a,b)13,(c,d)4cd)root;");
  EXPECT_EQ(t->getNodeCount(), 7);
  EXPECT_EQ(t->getTipCount(), 4);
  EXPECT_EQ(t->getInternalCount(), 3);

  EXPECT_EQ(t->getNode(0)->getName(), "root");
  EXPECT_EQ(t->getNode(1)->getName(), "13");
  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_EQ(t->getNode(4)->getName(), "4cd");
  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_EQ(t->getNode(6)->getName(), "d");
}

TEST(TreeReader, simple3) {
  auto tr = TreeReader();
  auto t = tr.readTree("((a:30.5,b:0.03):48.0,(c:0,d:3)cd)root;");
  EXPECT_EQ(t->getNodeCount(), 7);
  EXPECT_EQ(t->getTipCount(), 4);

  EXPECT_EQ(t->getNode(0)->getName(), "root");
  EXPECT_DOUBLE_EQ(t->getNode(0)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(1)->getName(), "");
  EXPECT_DOUBLE_EQ(t->getNode(1)->getBL(), 48.0);

  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_DOUBLE_EQ(t->getNode(2)->getBL(), 30.5);

  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0.03);

  EXPECT_EQ(t->getInternalCount(), 3);

  EXPECT_EQ(t->getNode(4)->getName(), "cd");
  EXPECT_DOUBLE_EQ(t->getNode(4)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_DOUBLE_EQ(t->getNode(5)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(6)->getName(), "d");
  EXPECT_DOUBLE_EQ(t->getNode(6)->getBL(), 3.0);
}

TEST(TreeReader, simple4) {
  auto tr = TreeReader();
  auto t = tr.readTree("((a:1e-10,b:0.03)ab,(c:0,d:3E-5)cd)root;");
  EXPECT_EQ(t->getNodeCount(), 7);
  EXPECT_EQ(t->getTipCount(), 4);

  EXPECT_EQ(t->getNode(0)->getName(), "root");
  EXPECT_DOUBLE_EQ(t->getNode(0)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(1)->getName(), "ab");
  EXPECT_DOUBLE_EQ(t->getNode(1)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_DOUBLE_EQ(t->getNode(2)->getBL(), 1e-10);

  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0.03);

  EXPECT_EQ(t->getInternalCount(), 3);

  EXPECT_EQ(t->getNode(4)->getName(), "cd");
  EXPECT_DOUBLE_EQ(t->getNode(4)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_DOUBLE_EQ(t->getNode(5)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(6)->getName(), "d");
  EXPECT_DOUBLE_EQ(t->getNode(6)->getBL(), 3e-5);
}

TEST(TreeReader, simple5) {
  auto tr = TreeReader();
  auto t = tr.readTree("( (a : 30.5 , b : 0.03 ) ab , (c :0,d : 3 ) cd )root;");
  EXPECT_EQ(t->getNodeCount(), 7);
  EXPECT_EQ(t->getTipCount(), 4);

  EXPECT_EQ(t->getNode(0)->getName(), "root");
  EXPECT_DOUBLE_EQ(t->getNode(0)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(1)->getName(), "ab");
  EXPECT_DOUBLE_EQ(t->getNode(1)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_DOUBLE_EQ(t->getNode(2)->getBL(), 30.5);

  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0.03);

  EXPECT_EQ(t->getInternalCount(), 3);

  EXPECT_EQ(t->getNode(4)->getName(), "cd");
  EXPECT_DOUBLE_EQ(t->getNode(4)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_DOUBLE_EQ(t->getNode(5)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(6)->getName(), "d");
  EXPECT_DOUBLE_EQ(t->getNode(6)->getBL(), 3.0);
}

TEST(TreeReader, simple6) {
  auto tr = TreeReader();
  auto t = tr.readTree(
      "(\t(a\t:\t30.5\t,\tb\t:\t0.03\t)\tab\t,\t(c\t:0,d\t:"
      "\t3\t)\tcd\t)root\t;");
  EXPECT_EQ(t->getNodeCount(), 7);
  EXPECT_EQ(t->getTipCount(), 4);

  EXPECT_EQ(t->getNode(0)->getName(), "root");
  EXPECT_DOUBLE_EQ(t->getNode(0)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(1)->getName(), "ab");
  EXPECT_DOUBLE_EQ(t->getNode(1)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_DOUBLE_EQ(t->getNode(2)->getBL(), 30.5);

  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0.03);

  EXPECT_EQ(t->getInternalCount(), 3);

  EXPECT_EQ(t->getNode(4)->getName(), "cd");
  EXPECT_DOUBLE_EQ(t->getNode(4)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_DOUBLE_EQ(t->getNode(5)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(6)->getName(), "d");
  EXPECT_DOUBLE_EQ(t->getNode(6)->getBL(), 3.0);
}

TEST(TreeReader, simple7) {
  auto tr = TreeReader();
  auto t = tr.readTree("((a,b)ab\n,(c,d\n)cd)root;");
  EXPECT_EQ(t->getNodeCount(), 7);
  EXPECT_EQ(t->getTipCount(), 4);

  EXPECT_EQ(t->getNode(0)->getName(), "root");
  EXPECT_DOUBLE_EQ(t->getNode(0)->getBL(), 0);

  EXPECT_EQ(t->getNode(1)->getName(), "ab");
  EXPECT_DOUBLE_EQ(t->getNode(1)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_DOUBLE_EQ(t->getNode(2)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0);

  EXPECT_EQ(t->getInternalCount(), 3);

  EXPECT_EQ(t->getNode(4)->getName(), "cd");
  EXPECT_DOUBLE_EQ(t->getNode(4)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_DOUBLE_EQ(t->getNode(5)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(6)->getName(), "d");
  EXPECT_DOUBLE_EQ(t->getNode(6)->getBL(), 0.0);
}

TEST(TreeReader, simple8) {
  auto tr = TreeReader();
  auto t = tr.readTree("((a,b)ab\r\n,(c,d\r\n)cd)root;");
  EXPECT_EQ(t->getNodeCount(), 7);
  EXPECT_EQ(t->getTipCount(), 4);

  EXPECT_EQ(t->getNode(0)->getName(), "root");
  EXPECT_DOUBLE_EQ(t->getNode(0)->getBL(), 0);

  EXPECT_EQ(t->getNode(1)->getName(), "ab");
  EXPECT_DOUBLE_EQ(t->getNode(1)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_DOUBLE_EQ(t->getNode(2)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0);

  EXPECT_EQ(t->getInternalCount(), 3);

  EXPECT_EQ(t->getNode(4)->getName(), "cd");
  EXPECT_DOUBLE_EQ(t->getNode(4)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_DOUBLE_EQ(t->getNode(5)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(6)->getName(), "d");
  EXPECT_DOUBLE_EQ(t->getNode(6)->getBL(), 0.0);
}

TEST(TreeReader, simple9) {
  auto tr = TreeReader();
  auto t = tr.readTree("((!a+7=5,b^o&)ab,($$£*c,d/\\?!_-|)cd)ro#~ot;");
  EXPECT_EQ(t->getTipCount(), 4);
  EXPECT_EQ(t->getNodeCount(), 7);

  EXPECT_EQ(t->getNode(0)->getName(), "ro#~ot");
  EXPECT_DOUBLE_EQ(t->getNode(0)->getBL(), 0);

  EXPECT_EQ(t->getNode(1)->getName(), "ab");
  EXPECT_DOUBLE_EQ(t->getNode(1)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(2)->getName(), "!a+7=5");
  EXPECT_DOUBLE_EQ(t->getNode(2)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(3)->getName(), "b^o&");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0);

  EXPECT_EQ(t->getInternalCount(), 3);

  EXPECT_EQ(t->getNode(4)->getName(), "cd");
  EXPECT_DOUBLE_EQ(t->getNode(4)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(5)->getName(), "$$£*c");
  EXPECT_DOUBLE_EQ(t->getNode(5)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(6)->getName(), "d/\\?!_-|");
  EXPECT_DOUBLE_EQ(t->getNode(6)->getBL(), 0.0);
}

TEST(TreeReader, simple10) {
  auto tr = TreeReader();
  auto t = tr.readTree("((a[comment],b),(c,(d, e):0.5));");
  EXPECT_EQ(t->getTipCount(), 5);
  EXPECT_EQ(t->getNodeCount(), 9);

  EXPECT_EQ(t->getNode(0)->getName(), "");
  EXPECT_DOUBLE_EQ(t->getNode(0)->getBL(), 0);

  EXPECT_EQ(t->getNode(1)->getName(), "");
  EXPECT_DOUBLE_EQ(t->getNode(1)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_DOUBLE_EQ(t->getNode(2)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0);

  EXPECT_EQ(t->getNode(4)->getName(), "");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0);

  EXPECT_EQ(t->getInternalCount(), 4);

  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_DOUBLE_EQ(t->getNode(5)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(6)->getName(), "");
  EXPECT_DOUBLE_EQ(t->getNode(6)->getBL(), 0.5);

  EXPECT_EQ(t->getNode(7)->getName(), "d");
  EXPECT_DOUBLE_EQ(t->getNode(7)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(8)->getName(), "e");
  EXPECT_DOUBLE_EQ(t->getNode(8)->getBL(), 0.0);
}

TEST(TreeReader, simple11) {
  auto tr = TreeReader();
  auto t = tr.readTree("((a,b),(c,(d, e)hello world));");
  EXPECT_EQ(t->getTipCount(), 5);
  EXPECT_EQ(t->getNodeCount(), 9);

  EXPECT_EQ(t->getNode(0)->getName(), "");
  EXPECT_DOUBLE_EQ(t->getNode(0)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(1)->getName(), "");
  EXPECT_DOUBLE_EQ(t->getNode(1)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(2)->getName(), "a");
  EXPECT_DOUBLE_EQ(t->getNode(2)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(3)->getName(), "b");
  EXPECT_DOUBLE_EQ(t->getNode(3)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(4)->getName(), "");
  EXPECT_DOUBLE_EQ(t->getNode(4)->getBL(), 0.0);

  EXPECT_EQ(t->getInternalCount(), 4);

  EXPECT_EQ(t->getNode(5)->getName(), "c");
  EXPECT_DOUBLE_EQ(t->getNode(5)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(6)->getName(), "hello world");
  EXPECT_DOUBLE_EQ(t->getNode(6)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(7)->getName(), "d");
  EXPECT_DOUBLE_EQ(t->getNode(7)->getBL(), 0.0);

  EXPECT_EQ(t->getNode(8)->getName(), "e");
  EXPECT_DOUBLE_EQ(t->getNode(8)->getBL(), 0.0);
}

TEST(TreeReader, badtrees1) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b),(c,(d, e)))"), std::runtime_error);
}

TEST(TreeReader, badtrees2) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b);,(c,(d, e)))"), std::runtime_error);
}

TEST(TreeReader, badtrees3) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b)(c,(d, e):0.5));"), std::runtime_error);
}

TEST(TreeReader, badtrees4) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b),(c,(d, e)));wtf"), std::runtime_error);
}

TEST(TreeReader, badtrees5) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b),(c,(d, e:0.1));"), std::runtime_error);
}

TEST(TreeReader, badtrees6) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("(a,b),(c,(d, e:0.1)));"), std::runtime_error);
}

TEST(TreeReader, badtrees7) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b),(c,(d, ())));"), std::runtime_error);
}

TEST(TreeReader, badtrees8) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b),(c,(d, (e))));"), std::runtime_error);
}

TEST(TreeReader, badtrees9) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b),(c,(d, e):0.5 label));"),
               std::runtime_error);
}

TEST(TreeReader, badtrees10) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b),(c,(d, e):0.a5));"), std::runtime_error);
}

TEST(TreeReader, badtrees11) {
  auto tr = TreeReader();
  EXPECT_THROW(tr.readTree("((a,b),(c,(d, e:0.0:0.1)));"), std::runtime_error);
}

TEST(TreeReader, manytrees) {
  auto treefile = env->get_datafile();

  if (!treefile.is_open()) {
    throw std::runtime_error{"Could not find data file"};
  }

  size_t line_number = 1;
  std::string line;
  while (std::getline(treefile, line)) {
    if (!treefile) { break; }
    auto tr = TreeReader();
    try {
      auto t = tr.readTree(line);
      EXPECT_GT(t->getTipCount(), 0);
    } catch (const std::exception &e) {
      throw std::runtime_error{std::string("Got error on line ")
                               + std::to_string(line_number)
                               + " with the error: \n" + e.what()};
    }
    line_number++;
    line.clear();
  }
}

TEST(TreeReader, DISABLED_pathologic0) {
  auto treefile = env->get_pathological_data();

  if (!treefile.is_open()) {
    throw std::runtime_error{"Could not find data file"};
  }

  size_t line_number = 1;
  std::string line;
  while (std::getline(treefile, line)) {
    if (!treefile) { break; }
    auto tr = TreeReader();
    try {
      auto t = tr.readTree(line);
      EXPECT_GT(t->getTipCount(), 0);
    } catch (const std::exception &e) {
      throw std::runtime_error{std::string("Got error on line ")
                               + std::to_string(line_number)
                               + " with the error: \n" + e.what()};
    }
    line_number++;
    line.clear();
  }
}

TEST(TreeReader, SlothsTree) {
  auto treefile = env->get_sloth_tree();

  if (!treefile.is_open()) {
    throw std::runtime_error{"Could not find data file"};
  }

  size_t line_number = 1;
  std::string line;
  while (std::getline(treefile, line)) {
    if (!treefile) { break; }
    auto tr = TreeReader();
    try {
      auto t = tr.readTree(line);
      EXPECT_GT(t->getTipCount(), 0);
    } catch (const std::exception &e) {
      throw std::runtime_error{std::string("Got error on line ")
                               + std::to_string(line_number)
                               + " with the error: \n" + e.what()};
    }
    line_number++;
    line.clear();
  }
}
