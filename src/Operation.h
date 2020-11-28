/* Operation.h
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#ifndef _LAGRANGE_OPERATION_H
#define _LAGRANGE_OPERATION_H

#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "Common.h"
#include "RateModel.h"
#include "Workspace.h"

class MakeRateMatrixOperation {
 public:
  MakeRateMatrixOperation(size_t index)
      : _rate_matrix_index{index},
        _period_index{0},
        _last_execution{0},
        _last_update{0} {}

  void eval(std::shared_ptr<Workspace> ws);

  lagrange_clock_t last_update() const {
    return std::max(_last_execution, _last_update);
  }

  inline void update_rates(const std::shared_ptr<Workspace>& ws, double disp,
                           double ext) {
    ws->set_period_params(_period_index, disp, ext);
    _last_update = ws->advance_clock();
  }

  inline void update_rates(const std::shared_ptr<Workspace>& ws, period_t p) {
    update_rates(ws, p.dispersion_rate, p.extinction_rate);
  }

  inline size_t rate_matrix_index() const { return _rate_matrix_index; }

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  std::string printStatus(const std::shared_ptr<Workspace>& ws,
                          size_t tabLevel = 0) const;

 private:
  size_t _rate_matrix_index;
  size_t _period_index;
  lagrange_clock_t _last_execution;
  lagrange_clock_t _last_update;
};

class ExpmOperation {
 public:
  ExpmOperation(size_t prob_matrix, double t,
                const std::shared_ptr<MakeRateMatrixOperation>& rm_op,
                bool transpose)
      : _prob_matrix_index{prob_matrix},
        _rate_matrix_index{rm_op->rate_matrix_index()},
        _t{t},
        _last_execution{0},
        _rate_matrix_op{rm_op},
        _transposed{transpose} {}

  ExpmOperation(size_t prob_matrix, double t,
                const std::shared_ptr<MakeRateMatrixOperation>& rm_op)
      : ExpmOperation(prob_matrix, t, rm_op, false) {}

  void eval(std::shared_ptr<Workspace> ws);

  size_t prob_matrix() const { return _prob_matrix_index; }

  lagrange_clock_t last_execution() const { return _last_execution; }

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;
  std::string printStatus(const std::shared_ptr<Workspace>& ws,
                          size_t tabLevel = 0) const;

 private:
  size_t _prob_matrix_index;
  size_t _rate_matrix_index;
  double _t;
  lagrange_clock_t _last_execution;
  std::shared_ptr<MakeRateMatrixOperation> _rate_matrix_op;
  bool _transposed;
};

class DispersionOperation {
 public:
  DispersionOperation(size_t top, size_t bot, double brlen, size_t prob_matrix,
                      const std::shared_ptr<MakeRateMatrixOperation>& rm_op)
      : _top_clv{top},
        _bot_clv{bot},
        _prob_matrix_index{prob_matrix},
        _expm_op{new ExpmOperation{prob_matrix, brlen, rm_op}} {}

  DispersionOperation(size_t bot, double brlen, size_t prob_matrix,
                      const std::shared_ptr<MakeRateMatrixOperation>& rm_op)
      : _top_clv{std::numeric_limits<size_t>::max()},
        _bot_clv{bot},
        _prob_matrix_index{prob_matrix},
        _expm_op{new ExpmOperation{prob_matrix, brlen, rm_op}} {}

  DispersionOperation(size_t bot, double brlen, size_t prob_matrix,
                      const std::shared_ptr<MakeRateMatrixOperation>& rm_op,
                      bool transpose)
      : _top_clv{std::numeric_limits<size_t>::max()},
        _bot_clv{bot},
        _prob_matrix_index{prob_matrix},
        _expm_op{new ExpmOperation{prob_matrix, brlen, rm_op, transpose}} {}

  DispersionOperation(size_t top, size_t bot,
                      const std::shared_ptr<ExpmOperation>& op)
      : _top_clv{top},
        _bot_clv{bot},
        _prob_matrix_index{op->prob_matrix()},
        _expm_op{op} {}

  DispersionOperation(size_t top, size_t bot, size_t prob_matrix)
      : _top_clv{top},
        _bot_clv{bot},
        _prob_matrix_index{prob_matrix},
        _expm_op{nullptr} {}

  void eval(std::shared_ptr<Workspace> ws);

  void terminate(size_t clv_index) {
    if (_top_clv != std::numeric_limits<size_t>::max()) {
      throw std::runtime_error{"Dispersion Operation already terminated"};
    }
    _top_clv = clv_index;
  }

  size_t top_clv_index() const { return _top_clv; };
  size_t bot_clv_index() const { return _bot_clv; };

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  std::string printStatus(const std::shared_ptr<Workspace>& ws,
                          size_t tabLevel = 0) const;

 private:
  /* Remember, the top and bottom clv indexes could be the same. This is to save
   * on storage when computing different periods along a single branch. The idea
   * is that we just apply the matrix multiplication over and over again to the
   * same CLV. Matvecs can't be done in place, so there will need to be
   * temporary storage if we do that.
   */
  size_t _top_clv;
  size_t _bot_clv;
  size_t _prob_matrix_index;

  /* _expm_op is an owning pointer, and can be null */
  std::shared_ptr<ExpmOperation> _expm_op;
};

class SplitOperation {
 public:
  SplitOperation(size_t lchild_clv_top, size_t rchild_clv_top,
                 size_t lchild_clv_bot, size_t rchild_clv_bot, double lbrlen,
                 double rbrlen, size_t lprob_mat, size_t rprob_mat,
                 std::shared_ptr<MakeRateMatrixOperation> lrate_matrix,
                 std::shared_ptr<MakeRateMatrixOperation> rrate_matrix,
                 size_t parent_clv)
      : _lbranch_clv_index{lchild_clv_top},
        _rbranch_clv_index{rchild_clv_top},
        _parent_clv_index{parent_clv},
        _lbranch_ops{{std::make_shared<DispersionOperation>(
            lchild_clv_top, lchild_clv_bot, lbrlen, lprob_mat, lrate_matrix)}},
        _rbranch_ops{{std::make_shared<DispersionOperation>(
            rchild_clv_top, rchild_clv_bot, rbrlen, rprob_mat, rrate_matrix)}} {
  }

  SplitOperation(size_t lchild_clv_top, size_t lchild_clv_bot,
                 size_t rchild_clv_top, size_t rchild_clv_bot, double lbrlen,
                 double rbrlen, size_t prob_mat,
                 std::shared_ptr<MakeRateMatrixOperation> rate_matrix,
                 size_t parent_clv)
      : SplitOperation(lchild_clv_top, rchild_clv_top, lchild_clv_bot,
                       rchild_clv_bot, lbrlen, rbrlen, prob_mat, prob_mat,
                       rate_matrix, rate_matrix, parent_clv) {}

  // NOTE: Fix this for periods
  SplitOperation(size_t parent_clv, std::shared_ptr<DispersionOperation> l_ops,
                 std::shared_ptr<DispersionOperation> r_ops)
      : _lbranch_clv_index{l_ops->top_clv_index()},
        _rbranch_clv_index{r_ops->top_clv_index()},
        _parent_clv_index{parent_clv},
        _lbranch_ops{{l_ops}},
        _rbranch_ops{{r_ops}} {}

  void eval(std::shared_ptr<Workspace>) const;

  size_t get_parent_clv() const { return _parent_clv_index; }

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  std::string printStatus(const std::shared_ptr<Workspace>& ws,
                          size_t tabLevel = 0) const;

 private:
  size_t _lbranch_clv_index;
  size_t _rbranch_clv_index;
  size_t _parent_clv_index;

  std::vector<std::shared_ptr<DispersionOperation>> _lbranch_ops;
  std::vector<std::shared_ptr<DispersionOperation>> _rbranch_ops;

  std::vector<lagrange_dist_t> _excl_dists;
};

class ReverseSplitOperation {
 public:
  ReverseSplitOperation(
      size_t bot_clv,  /* Where the result is stored*/
      size_t ltop_clv, /* Where the result of the dispop is stored */
      size_t rtop_clv, /* Should be a _non_ reverse CLV */
      std::shared_ptr<MakeRateMatrixOperation> rate_matrix_op,
      size_t prob_matrix_index, size_t disp_clv_index, double brlen)
      : _bot_clv_index{bot_clv},
        _ltop_clv_index{ltop_clv},
        _rtop_clv_index{rtop_clv},
        _eval_clvs{true},
        _branch_ops{{std::make_shared<DispersionOperation>(
            ltop_clv, disp_clv_index, brlen, prob_matrix_index,
            rate_matrix_op)}} {}

  ReverseSplitOperation(size_t bot_clv, size_t rtop_clv,
                        std::shared_ptr<DispersionOperation> branch_op)
      : _bot_clv_index{bot_clv},
        _ltop_clv_index{branch_op->bot_clv_index()},
        _rtop_clv_index{rtop_clv},
        _eval_clvs{true},
        _branch_ops{{branch_op}} {}

  ReverseSplitOperation(std::shared_ptr<DispersionOperation> branch_op)
      : _bot_clv_index{std::numeric_limits<size_t>::max()},
        _ltop_clv_index{std::numeric_limits<size_t>::max()},
        _rtop_clv_index{std::numeric_limits<size_t>::max()},
        _eval_clvs{false},
        _branch_ops{{branch_op}} {}

  void eval(std::shared_ptr<Workspace>) const;

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  std::string printStatus(const std::shared_ptr<Workspace>& ws,
                          size_t tabLevel = 0) const;

  void makeRootOperation(size_t clv_index) {
    _branch_ops.clear();
    if (_ltop_clv_index == std::numeric_limits<size_t>::max()) {
      _ltop_clv_index = clv_index;
    }
    if (_rtop_clv_index == std::numeric_limits<size_t>::max()) {
      _rtop_clv_index = clv_index;
    }
  }

  size_t getStableCLV() const{
    return _ltop_clv_index;
  }

 private:
  size_t _bot_clv_index;
  size_t _ltop_clv_index;
  size_t _rtop_clv_index;
  bool _eval_clvs;

  std::vector<std::shared_ptr<DispersionOperation>> _branch_ops;
  std::vector<lagrange_dist_t> _excl_dists;
};

class LHGoal {
 public:
  LHGoal(size_t root_clv, size_t prior_index)
      : _root_clv_index{root_clv}, _prior_index{prior_index} {}
  double eval(std::shared_ptr<Workspace>) const;

  size_t _root_clv_index;
  size_t _prior_index;

 private:
};

class StateLHGoal {
 public:
  StateLHGoal(size_t parent_clv, size_t lchild_clv, size_t rchild_clv)
      : _parent_clv_index{parent_clv},
        _lchild_clv_index{lchild_clv},
        _rchild_clv_index{rchild_clv} {}

  lagrange_col_vector_t eval(std::shared_ptr<Workspace>) const;

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;
};

class SplitLHGoal {
 public:
  std::unordered_map<lagrange_dist_t, std::vector<AncSplit>> eval(
      std::shared_ptr<Workspace> ws) const;

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;
};

#endif
