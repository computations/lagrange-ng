#include "Operation.hpp"

#include <cmath>
#include <memory>

#include "Common.hpp"
#include "Workspace.hpp"
#include "gtest/gtest.h"

using namespace lagrange;

#define lagrange_compute_error_from_matrix(calc_buffer, ref_buffer, workspace) \
  {                                                                            \
    double error = 0.0;                                                        \
    for (size_t i = 0; i < workspace->matrixSize(); i++) {                     \
      double cur_error = calc_buffer[i] - ref_buffer.get()[i];                 \
      error += cur_error * cur_error;                                          \
    }                                                                          \
    error = std::sqrt(error);                                                  \
    EXPECT_NEAR(error, 0.0, 1e-7);                                             \
  }

#define lagrange_compute_error_from_vector(calc_buffer, ref_buffer, workspace) \
  {                                                                            \
    double error = 0.0;                                                        \
    for (size_t i = 0; i < workspace->states(); i++) {                         \
      double cur_error = calc_buffer[i] - ref_buffer.get()[i];                 \
      error += cur_error * cur_error;                                          \
    }                                                                          \
    error = std::sqrt(error);                                                  \
    EXPECT_NEAR(error, 0.0, 1e-7);                                             \
  }

class OperationTest : public ::testing::Test {
 protected:
  void SetUp() override {
    /*
    _arbitrary_rate_matrix = {
        {-2.426898, 1.094290, 0.849836, 0.482772},
        {0.088512, -1.800944, 0.083282, 1.629150},
        {0.759663, 0.547961, -1.367337, 0.059712},
        {0.360901, 0.374056, 1.267820, -2.002777},
    };
    */

    _arbitrary_rate_matrix.reset(new LagrangeMatrixBase[4 * 4]{
        -2.426898,
        1.094290,
        0.849836,
        0.482772,
        0.088512,
        -1.800944,
        0.083282,
        1.629150,
        0.759663,
        0.547961,
        -1.367337,
        0.059712,
        0.360901,
        0.374056,
        1.267820,
        -2.002777,
    });

    /*
    _correct_prob_matrix = {
        {0.1949195896, 0.2632718648, 0.3005026902, 0.24130565},
        {0.110711104, 0.3018364877, 0.2646053177, 0.3228469542},
        {0.1710689037, 0.2362891503, 0.4115737324, 0.1810676026},
        {0.1494294987, 0.21110698, 0.3574843286, 0.2819789219},
    };
    */

    _correct_prob_matrix.reset(new LagrangeMatrixBase[4 * 4]{
        0.1949195896,
        0.2632718648,
        0.3005026902,
        0.24130565,
        0.110711104,
        0.3018364877,
        0.2646053177,
        0.3228469542,
        0.1710689037,
        0.2362891503,
        0.4115737324,
        0.1810676026,
        0.1494294987,
        0.21110698,
        0.3574843286,
        0.2819789219,
    });

    /* _correct_root_clv = {
     *  0, 0.04839594797,        0.06003722477,        0.05381024353};
        0, 0.048395948120807594, 0.060037224600835798, 0.053855745689636927});
   */
    /*_correct_root_clv.reset(new LagrangeMatrixBase[4]{*/
    /*    0, 0.04839594797, 0.06003722477, 0.05381024353});*/

    _correct_root_clv.reset(new LagrangeMatrixBase[4]{
        0, 0.048395948120807594, 0.060037224600835798, 0.053855745689636927});
    /* _correct_reverse_bot_clv = {0, 0.04898326741, 0.04287039582,
     * 0.01129551391};
     */

    /*_correct_reverse_bot_clv.reset(new LagrangeMatrixBase[4]{*/
    /*    0, 0.04898326741, 0.04287039582, 0.01129551391});*/
    
    _correct_reverse_bot_clv.reset(new LagrangeMatrixBase[4]{
        0, 0.049699709788676508, 0.044376375604528666, 0.0084716354377575281});
    /*
      _arbitrary_clv1 = {
          0.491885,
          0.180252,
          0.136036,
          0.191827,
      };
     */

    _arbitrary_clv1.reset(new LagrangeMatrixBase[4]{
        0.491885,
        0.180252,
        0.136036,
        0.191827,
    });

    /*
    _arbitrary_clv2 = {
        0.2921,
        0.297896,
        0.353511,
        0.0564923,
    };
    */

    _arbitrary_clv2.reset(new LagrangeMatrixBase[4]{
        0.2921,
        0.297896,
        0.353511,
        0.0564923,
    });

    _ws = std::make_shared<Workspace>(_taxa, _regions, _regions);

    _lbot_clv = _ws->registerGenericCLV();
    _ltop_clv = _ws->registerGenericCLV();
    _rbot_clv = _ws->registerGenericCLV();
    _rtop_clv = _ws->registerGenericCLV();

    _root_clv = _ws->registerGenericCLV();

    _reverse_bot_clv = _ws->registerGenericCLV();
    _reverse_ltop_clv = _ws->registerGenericCLV();
    _reverse_lbot_clv = _ws->registerGenericCLV();

    _ws->reserve();

    _ws->updateRateMatrix(_rate_matrix, _arbitrary_rate_matrix.get());

    _ws->updateCLV(_arbitrary_clv1.get(), _rbot_clv);
    _ws->updateCLV(_arbitrary_clv2.get(), _lbot_clv);

    _ws->updateCLV(_arbitrary_clv1.get(), _reverse_lbot_clv);
    _ws->setPeriodParams(0, .3123, 1.1231);

    _rate_matrix_op =
        std::make_shared<MakeRateMatrixOperation>(_rate_matrix, _period);
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
  size_t _period = 0;

  std::unique_ptr<LagrangeMatrixBase[]> _arbitrary_rate_matrix;
  std::unique_ptr<LagrangeMatrixBase[]> _correct_prob_matrix;
  std::unique_ptr<LagrangeMatrixBase[]> _arbitrary_clv1;
  std::unique_ptr<LagrangeMatrixBase[]> _arbitrary_clv2;
  std::unique_ptr<LagrangeMatrixBase[]> _correct_root_clv;
  std::unique_ptr<LagrangeMatrixBase[]> _correct_reverse_bot_clv;

  std::shared_ptr<Workspace> _ws;
  std::shared_ptr<MakeRateMatrixOperation> _rate_matrix_op;
};

TEST_F(OperationTest, ExpmSimple0) {
  ExpmOperation expm_op(_prob_matrix, _t, _rate_matrix_op);
  expm_op.eval(_ws);
}

TEST_F(OperationTest, ExpmSimple1) {
  ExpmOperation expm_op(_prob_matrix, _t, _rate_matrix_op);

  expm_op.eval(_ws);

  lagrange_compute_error_from_matrix(
      _ws->probMatrix(_prob_matrix), _correct_prob_matrix, _ws);
}

TEST_F(OperationTest, DispersionSimple0) {
  auto expm_op =
      std::make_shared<ExpmOperation>(_prob_matrix, _t, _rate_matrix_op);
  DispersionOperation disp_op(_rtop_clv, _rbot_clv, expm_op);
  disp_op.fallback();
  disp_op.eval(_ws);
}

TEST_F(OperationTest, DispersionSimple1) {
  auto expm_op =
      std::make_shared<ExpmOperation>(_prob_matrix, _t, _rate_matrix_op);

  DispersionOperation disp_op(_rtop_clv, _rbot_clv, expm_op);
  disp_op.eval(_ws);

  /*
  double error =
      blaze::norm(_ws->prob_matrix(_prob_matrix) - _correct_prob_matrix);
  */

  lagrange_compute_error_from_vector(
      _ws->probMatrix(_prob_matrix), _correct_prob_matrix, _ws);
}

/* Regression */
TEST_F(OperationTest, DispersionSimple2) {
  auto expm_op =
      std::make_shared<ExpmOperation>(_prob_matrix, _t, _rate_matrix_op);
  DispersionOperation disp_op(_rtop_clv, _rbot_clv, expm_op);
  disp_op.eval(_ws);

  lagrange_compute_error_from_vector(
      _ws->probMatrix(_prob_matrix), _correct_prob_matrix, _ws);

  std::unique_ptr<LagrangeMatrixBase[]> correct_clv(new LagrangeMatrixBase[4]{
      0.2305014254, 0.2067903737, 0.2174603189, 0.2142764931});

  lagrange_compute_error_from_vector(_ws->CLV(_rtop_clv), correct_clv, _ws);
}

TEST_F(OperationTest, SplitSimple0) {
  SplitOperation split_op(_ltop_clv,
                          _lbot_clv,
                          _rtop_clv,
                          _rbot_clv,
                          _t,
                          _t,
                          _prob_matrix,
                          _rate_matrix_op,
                          _root_clv);

  split_op.eval(_ws);

  lagrange_compute_error_from_vector(
      _ws->CLV(_root_clv), _correct_root_clv, _ws);
}

TEST_F(OperationTest, ReverseSplitSimple0) {
  ReverseSplitOperation rsplit_op(_reverse_bot_clv,
                                  _reverse_ltop_clv,
                                  _rbot_clv,
                                  _rate_matrix_op,
                                  _prob_matrix,
                                  _reverse_lbot_clv,
                                  _t);

  rsplit_op.eval(_ws);

  lagrange_compute_error_from_vector(
      _ws->CLV(_reverse_bot_clv), _correct_reverse_bot_clv, _ws);
}

TEST_F(OperationTest, MakeRateMatrixOperationSimple0) {
  MakeRateMatrixOperation make_op(_rate_matrix, _period);

  make_op.eval(_ws);
}

TEST_F(OperationTest, MakeRateMatrixOperationSimple1) {
  auto local_ws = std::make_shared<Workspace>(_taxa, 3, 3);
  local_ws->reserve();
  local_ws->setPeriodParams(_period, .3123, 1.1231);
  MakeRateMatrixOperation make_op(_rate_matrix, _period);

  make_op.eval(local_ws);
}
