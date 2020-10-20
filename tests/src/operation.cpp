#include <memory>
#include <stdexcept>

#include "Operation.h"
#include "RateModel.h"
#include "Workspace.h"
#include "environment.hpp"
#include "gtest/gtest.h"

class OperationTest : public ::testing::Test {
protected:
  void SetUp() override {
    _arbitrary_rate_matrix = {
        {-2.426898, 1.094290, 0.849836, 0.482772},
        {0.088512, -1.800944, 0.083282, 1.629150},
        {0.759663, 0.547961, -1.367337, 0.059712},
        {0.360901, 0.374056, 1.267820, -2.002777},
    };

    _correct_prob_matrix = {
        {1.949195922460272e-01, 2.632718625819546e-01, 3.005026886371733e-01,
         2.413056510928967e-01},
        {1.107111041626160e-01, 3.018364882518785e-01, 2.646053198768748e-01,
         3.228469513771057e-01},
        {1.710689022436381e-01, 2.362891506596847e-01, 4.115737326038578e-01,
         1.810676035510108e-01},
        {1.494294989412788e-01, 2.111069803010049e-01, 3.574843271565785e-01,
         2.819789228671793e-01},
    };

    _correct_root_clv = {
        0,
        0.0483959,
        0.0600372,
        0.322861,
    };

    _arbitrary_clv1 = {
        0.491885,
        0.180252,
        0.136036,
        0.191827,
    };

    _arbitrary_clv2 = {
        0.2921,
        0.297896,
        0.353511,
        0.0564923,
    };

    _ws = make_shared<Workspace>(_taxa, _regions);
    _ws->rate_matrix(_rate_matrix) = _arbitrary_rate_matrix;
    _ws->clv(_rbot_clv) = _arbitrary_clv1;
    _ws->clv(_lbot_clv) = _arbitrary_clv2;
    _ws->clv(_root_clv) = {0, 0, 0, 0};
  }

  size_t _lbot_clv = 0;
  size_t _ltop_clv = 1;
  size_t _rbot_clv = 2;
  size_t _rtop_clv = 3;

  size_t _root_clv = 4;

  size_t _prob_matrix = 0;
  size_t _rate_matrix = 0;
  size_t _regions = 2;
  size_t _taxa = 2;

  double _t = 1.0;

  lagrange_matrix_t _arbitrary_rate_matrix;
  lagrange_matrix_t _correct_prob_matrix;
  lagrange_col_vector_t _arbitrary_clv1;
  lagrange_col_vector_t _arbitrary_clv2;
  lagrange_col_vector_t _correct_root_clv;
  std::shared_ptr<Workspace> _ws;
};

TEST_F(OperationTest, ExpmSimple0) {
  auto expm_op = ExpmOperation(_rate_matrix, _prob_matrix, _t);
  expm_op.eval(_ws);
}

TEST_F(OperationTest, ExpmSimple1) {
  auto expm_op = ExpmOperation(_rate_matrix, _prob_matrix, _t);

  expm_op.eval(_ws);

  auto diff = _ws->prob_matrix(_prob_matrix) - _correct_prob_matrix;
  double error = blaze::norm(diff);
  EXPECT_NEAR(error, 0.0, 1e-7);
}

TEST_F(OperationTest, DispersionSimple0) {
  auto ws = std::make_shared<Workspace>(_taxa, _regions);

  DispersionOperation disp_op(_rtop_clv, _rbot_clv, _prob_matrix);
  disp_op.eval(ws);
}

TEST_F(OperationTest, DispersionSimple1) {
  auto expm_op =
      std::make_shared<ExpmOperation>(_rate_matrix, _prob_matrix, _t);
  DispersionOperation disp_op(_rtop_clv, _rbot_clv, expm_op);
  disp_op.eval(_ws);

  double error =
      blaze::norm(_ws->prob_matrix(_prob_matrix) - _correct_prob_matrix);
  EXPECT_NEAR(error, 0.0, 1e-7);
}

TEST_F(OperationTest, DispersionSimple2) {
  auto expm_op =
      std::make_shared<ExpmOperation>(_rate_matrix, _prob_matrix, _t);
  DispersionOperation disp_op(_rtop_clv, _rbot_clv, expm_op);
  disp_op.eval(_ws);

  double error =
      blaze::norm(_ws->prob_matrix(_prob_matrix) - _correct_prob_matrix);
  EXPECT_NEAR(error, 0.0, 1e-7);

  lagrange_col_vector_t correct_clv = {
      2.305014262897032e-01,
      2.067903735879926e-01,
      2.174603184396995e-01,
      2.142764932658623e-01,
  };

  error = blaze::norm(_ws->clv(_rtop_clv) - correct_clv);

  EXPECT_NEAR(error, 0.0, 1e-7);
}

TEST_F(OperationTest, SplitSimple0) {
  SplitOperation split_op(_ltop_clv, _lbot_clv, _rtop_clv, _rbot_clv, _t, _t,
                          _prob_matrix, _rate_matrix, _root_clv);

  split_op.eval(_ws);

  auto &root_clv = _ws->clv(_root_clv);
  double error = blaze::norm(root_clv - _correct_root_clv);
  EXPECT_NEAR(error, 0.0, 1e-4);
}
