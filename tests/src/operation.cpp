#include <iomanip>
#include <memory>
#include <type_traits>
#include <unordered_map>

#include "Operation.h"
#include "Utils.h"
#include "Workspace.h"
#include "environment.hpp"
#include "gtest/gtest.h"

#define lagrange_test_check_blis_obj(calc_matrix, ref_matrix)            \
  {                                                                      \
    obj_t diff;                                                          \
    bli_obj_create_conf_to(&ref_matrix, &diff);                          \
    bli_copym(&ref_matrix, &diff);                                       \
    bli_subm(calc_matrix, &diff);                                        \
                                                                         \
    obj_t norm;                                                          \
    double error;                                                        \
    bli_obj_create_1x1_with_attached_buffer(BLIS_DOUBLE, &error, &norm); \
    bli_normfm(&diff, &norm);                                            \
    bli_obj_free(&diff);                                                 \
    EXPECT_NEAR(error, 0.0, 1e-7);                                       \
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

    lagrange_set_bli_matrix(&_arbitrary_rate_matrix, 4, 4,
                            {
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

    lagrange_set_bli_matrix(&_correct_prob_matrix, 4, 4,
                            {
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

    /* _correct_root_clv = {0, 0.04839594797, 0.06003722477, 0.05381024353}; */
    lagrange_set_bli_vector(&_correct_root_clv,
                            {0, 0.04839594797, 0.06003722477, 0.05381024353});

    /* _correct_reverse_bot_clv = {0, 0.04898326741, 0.04287039582,
     * 0.01129551391};
     */
    lagrange_set_bli_vector(&_correct_reverse_bot_clv,
                            {0, 0.04898326741, 0.04287039582, 0.01129551391});

    /*
      _arbitrary_clv1 = {
          0.491885,
          0.180252,
          0.136036,
          0.191827,
      };
     */

    lagrange_set_bli_vector(&_arbitrary_clv1, {
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

    lagrange_set_bli_vector(&_arbitrary_clv2, {
                                                  0.2921,
                                                  0.297896,
                                                  0.353511,
                                                  0.0564923,
                                              });

    _ws = std::make_shared<Workspace>(_taxa, _regions);

    _lbot_clv = _ws->register_generic_clv();
    _ltop_clv = _ws->register_generic_clv();
    _rbot_clv = _ws->register_generic_clv();
    _rtop_clv = _ws->register_generic_clv();

    _root_clv = _ws->register_generic_clv();

    _reverse_bot_clv = _ws->register_generic_clv();
    _reverse_ltop_clv = _ws->register_generic_clv();
    _reverse_lbot_clv = _ws->register_generic_clv();

    _ws->reserve();

    bli_copym(&_arbitrary_rate_matrix, _ws->rate_matrix(_rate_matrix));

    _ws->update_clv(&_arbitrary_clv1, _rbot_clv);
    _ws->update_clv(&_arbitrary_clv2, _lbot_clv);
    //_ws->update_clv({0, 0, 0, 0}, &_root_clv);

    //_ws->update_clv({0, 0, 0, 0}, _reverse_bot_clv);
    _ws->update_clv(&_arbitrary_clv1, _reverse_lbot_clv);
    _ws->set_period_params(0, .3123, 1.1231);

    _rate_matrix_op = std::make_shared<MakeRateMatrixOperation>(_rate_matrix);

    _blis_runtime = BLIS_RNTM_INITIALIZER;
    _blis_context = bli_gks_query_cntx();
  }

  void TearDown() override {
    bli_obj_free(&_arbitrary_rate_matrix);
    bli_obj_free(&_correct_prob_matrix);
    bli_obj_free(&_arbitrary_clv1);
    bli_obj_free(&_arbitrary_clv2);
    bli_obj_free(&_correct_root_clv);
    bli_obj_free(&_correct_reverse_bot_clv);
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

  lagrange_matrix_base_t _arbitrary_rate_matrix;
  lagrange_matrix_base_t _correct_prob_matrix;
  lagrange_matrix_base_t _arbitrary_clv1;
  lagrange_matrix_base_t _arbitrary_clv2;
  lagrange_matrix_base_t _correct_root_clv;
  lagrange_matrix_base_t _correct_reverse_bot_clv;

  std::shared_ptr<Workspace> _ws;
  std::shared_ptr<MakeRateMatrixOperation> _rate_matrix_op;

  cntx_t *_blis_context;
  rntm_t _blis_runtime;
};

TEST_F(OperationTest, ExpmSimple0) {
  ExpmOperation expm_op(_prob_matrix, _t, _rate_matrix_op);
  expm_op.eval(_ws, _blis_context, &_blis_runtime);
}

TEST_F(OperationTest, ExpmSimple1) {
  ExpmOperation expm_op(_prob_matrix, _t, _rate_matrix_op);

  expm_op.eval(_ws, _blis_context, &_blis_runtime);

  /*
  auto diff = _ws->prob_matrix(_prob_matrix) - _correct_prob_matrix;
  double error = blaze::norm(diff);
  */

  lagrange_test_check_blis_obj(_ws->prob_matrix(_prob_matrix),
                               _correct_prob_matrix);
}

TEST_F(OperationTest, DispersionSimple0) {
  DispersionOperation disp_op(_rtop_clv, _rbot_clv, _prob_matrix);
  disp_op.eval(_ws, _blis_context, &_blis_runtime);
}

TEST_F(OperationTest, DispersionSimple1) {
  auto expm_op =
      std::make_shared<ExpmOperation>(_prob_matrix, _t, _rate_matrix_op);

  DispersionOperation disp_op(_rtop_clv, _rbot_clv, expm_op);
  disp_op.eval(_ws, _blis_context, &_blis_runtime);

  /*
  double error =
      blaze::norm(_ws->prob_matrix(_prob_matrix) - _correct_prob_matrix);
  */

  lagrange_test_check_blis_obj(_ws->prob_matrix(_prob_matrix),
                               _correct_prob_matrix);
}

/* Regression */
TEST_F(OperationTest, DispersionSimple2) {
  auto expm_op =
      std::make_shared<ExpmOperation>(_prob_matrix, _t, _rate_matrix_op);
  DispersionOperation disp_op(_rtop_clv, _rbot_clv, expm_op);
  disp_op.eval(_ws, _blis_context, &_blis_runtime);

  lagrange_test_check_blis_obj(_ws->prob_matrix(_prob_matrix),
                               _correct_prob_matrix);

  lagrange_matrix_base_t correct_clv;
  lagrange_set_bli_vector(
      &correct_clv, {0.2305014254, 0.2067903737, 0.2174603189, 0.2142764931});

  lagrange_test_check_blis_obj(_ws->clv(_rtop_clv), correct_clv);

  bli_obj_free(&correct_clv);
}

TEST_F(OperationTest, SplitSimple0) {
  SplitOperation split_op(_ltop_clv, _lbot_clv, _rtop_clv, _rbot_clv, _t, _t,
                          _prob_matrix, _rate_matrix_op, _root_clv);

  split_op.eval(_ws, _blis_context, &_blis_runtime);

  lagrange_test_check_blis_obj(_ws->clv(_root_clv), _correct_root_clv);
}

TEST_F(OperationTest, ReverseSplitSimple0) {
  ReverseSplitOperation rsplit_op(_reverse_bot_clv, _reverse_ltop_clv,
                                  _rbot_clv, _rate_matrix_op, _prob_matrix,
                                  _reverse_lbot_clv, _t);

  rsplit_op.eval(_ws, _blis_context, &_blis_runtime);

  lagrange_test_check_blis_obj(_ws->clv(_reverse_bot_clv),
                               _correct_reverse_bot_clv);
}

TEST_F(OperationTest, MakeRateMatrixOperationSimple0) {
  MakeRateMatrixOperation make_op(_rate_matrix);

  make_op.eval(_ws);
}

TEST_F(OperationTest, MakeRateMatrixOperationSimple1) {
  auto local_ws = std::make_shared<Workspace>(_taxa, 3);
  local_ws->reserve();
  local_ws->set_period_params(0, .3123, 1.1231);
  MakeRateMatrixOperation make_op(_rate_matrix);

  make_op.eval(local_ws);
}
