#include <iomanip>
#include <memory>
#include <type_traits>

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
        {6.959163838915934e-01, 1.382996456185051e-01, 1.382996456185052e-01,
         2.748432487139660e-02},
        {4.973561703302694e-01, 3.368598591798291e-01, 9.883972226405867e-02,
         6.694424822584306e-02},
        {4.973561703302696e-01, 9.883972226405867e-02, 3.368598591798295e-01,
         6.694424822584309e-02},
        {3.554495423463473e-01, 2.407463502479806e-01, 2.407463502479807e-01,
         1.630577571576918e-01},

    };

    /*
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
      */

    _correct_root_clv = {
        0,
        0.09430434611,
        0.09556265233,
        0.08840423149,
    };

    _correct_reverse_bot_clv = {
        0,
        0.05978042766,
        0.04368454415,
        0.05291602473,
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

    _lbot_clv = _ws->register_generic_clv();
    _ltop_clv = _ws->register_generic_clv();
    _rbot_clv = _ws->register_generic_clv();
    _rtop_clv = _ws->register_generic_clv();

    _root_clv = _ws->register_generic_clv();

    _reverse_bot_clv = _ws->register_generic_clv();
    _reverse_ltop_clv = _ws->register_generic_clv();
    _reverse_lbot_clv = _ws->register_generic_clv();

    _ws->reserve();

    _ws->rate_matrix(_rate_matrix) = _arbitrary_rate_matrix;

    _ws->clv(_rbot_clv) = _arbitrary_clv1;
    _ws->clv(_lbot_clv) = _arbitrary_clv2;
    _ws->clv(_root_clv) = {0, 0, 0, 0};

    _ws->clv(_reverse_bot_clv) = {0, 0, 0, 0};
    _ws->clv(_reverse_lbot_clv) = _arbitrary_clv1;
    _ws->set_period_params(0, .3123, 1.1231);

    _rate_matrix_op = std::make_shared<MakeRateMatrixOperation>(_rate_matrix);
    //_rate_matrix_op->update_rates(_ws, .3123, 1.1231);
  }

  size_t _lbot_clv;
  size_t _ltop_clv;
  size_t _rbot_clv;
  size_t _rtop_clv;

  size_t _root_clv;

  size_t _reverse_bot_clv;
  size_t _reverse_ltop_clv;
  size_t _reverse_lbot_clv;

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
  lagrange_col_vector_t _correct_reverse_bot_clv;
  std::shared_ptr<Workspace> _ws;
  std::shared_ptr<MakeRateMatrixOperation> _rate_matrix_op;
};

TEST_F(OperationTest, ExpmSimple0) {
  auto expm_op = ExpmOperation(_prob_matrix, _t, _rate_matrix_op);
  expm_op.eval(_ws);
}

TEST_F(OperationTest, ExpmSimple1) {
  auto expm_op = ExpmOperation(_prob_matrix, _t, _rate_matrix_op);

  expm_op.eval(_ws);

  auto diff = _ws->prob_matrix(_prob_matrix) - _correct_prob_matrix;
  double error = blaze::norm(diff);
  EXPECT_NEAR(error, 0.0, 1e-7);
}

TEST_F(OperationTest, DispersionSimple0) {
  DispersionOperation disp_op(_rtop_clv, _rbot_clv, _prob_matrix);
  disp_op.eval(_ws);
}

TEST_F(OperationTest, DispersionSimple1) {
  auto expm_op =
      std::make_shared<ExpmOperation>(_prob_matrix, _t, _rate_matrix_op);

  DispersionOperation disp_op(_rtop_clv, _rbot_clv, expm_op);
  disp_op.eval(_ws);

  double error =
      blaze::norm(_ws->prob_matrix(_prob_matrix) - _correct_prob_matrix);
  EXPECT_NEAR(error, 0.0, 1e-7);
}

TEST_F(OperationTest, DispersionSimple2) {
  auto expm_op =
      std::make_shared<ExpmOperation>(_prob_matrix, _t, _rate_matrix_op);
  DispersionOperation disp_op(_rtop_clv, _rbot_clv, expm_op);
  disp_op.eval(_ws);

  double error =
      blaze::norm(_ws->prob_matrix(_prob_matrix) - _correct_prob_matrix);
  EXPECT_NEAR(error, 0.0, 1e-7);

  lagrange_col_vector_t correct_clv = {
      3.913255843910076e-01,
      3.316491779421193e-01,
      3.211248795682519e-01,
      2.822643601465549e-01,
  };

  error = blaze::norm(_ws->clv(_rtop_clv) - correct_clv);

  EXPECT_NEAR(error, 0.0, 1e-7);
}

TEST_F(OperationTest, SplitSimple0) {
  SplitOperation split_op(_ltop_clv, _lbot_clv, _rtop_clv, _rbot_clv, _t, _t,
                          _prob_matrix, _rate_matrix_op, _root_clv);

  split_op.eval(_ws);

  auto &root_clv = _ws->clv(_root_clv);
  double error = blaze::norm(root_clv - _correct_root_clv);
  EXPECT_NEAR(error, 0.0, 1e-4);
}

TEST_F(OperationTest, ReverseSplitSimple0) {
  ReverseSplitOperation rsplit_op(_reverse_bot_clv, _reverse_ltop_clv,
                                  _rbot_clv, _rate_matrix_op, _prob_matrix,
                                  _reverse_lbot_clv, _t);

  _ws->clv(_reverse_ltop_clv) = 0.0;

  rsplit_op.eval(_ws);

  auto &root_clv = _ws->clv(_reverse_bot_clv);
  double error = blaze::norm(root_clv - _correct_reverse_bot_clv);
  EXPECT_NEAR(error, 0.0, 1e-4);
}
