#include "ThreadState.h"
#include "gtest/gtest.h"

class ThreadStateTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(ThreadStateTest, simple0) {
  ThreadState ts;
  EXPECT_EQ(ts.thread_id(), 0);
}
