/* Operation.h
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#ifndef _LAGRANGE_OPERATION_H
#define _LAGRANGE_OPERATION_H

#include <limits>
#include <memory>
#include <stdexcept>

#include "Workspace.h"

class MakeRateMatrixOperation {
 public:
  void eval(Workspace &ws) const;

 private:
  size_t _rate_matrix_index;
  double _dispersion_rate;
  double _extinction_rate;
};

class ExpmOperation {
 public:
  ExpmOperation(size_t prob_matrix, size_t rate_matrix, double t,
                bool transpose)
      : _prob_matrix_index{prob_matrix},
        _rate_matrix_index{rate_matrix},
        _t{t},
        _transposed{transpose} {}

  ExpmOperation(size_t prob_matrix, size_t rate_matrix, double t)
      : ExpmOperation(prob_matrix, rate_matrix, t, false) {}

  void eval(std::shared_ptr<Workspace> ws);

  size_t prob_matrix() const { return _prob_matrix_index; }

 private:
  size_t _prob_matrix_index;
  size_t _rate_matrix_index;
  double _t;
  bool _transposed;
};

class DispersionOperation {
 public:
  DispersionOperation(size_t top, size_t bot, double brlen, size_t prob_matrix,
                      size_t rate_matrix)
      : _top_clv{top},
        _bot_clv{bot},
        _prob_matrix_index{prob_matrix},
        _expm_op{new ExpmOperation{prob_matrix, rate_matrix, brlen}} {}

  DispersionOperation(size_t bot, double brlen, size_t prob_matrix,
                      size_t rate_matrix)
      : _top_clv{std::numeric_limits<size_t>::max()},
        _bot_clv{bot},
        _prob_matrix_index{prob_matrix},
        _expm_op{new ExpmOperation{prob_matrix, rate_matrix, brlen}} {}

  DispersionOperation(size_t top, size_t bot, std::shared_ptr<ExpmOperation> op)
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
  SplitOperation(size_t lchild_clv, size_t rchild_clv, size_t parent_clv,
                 std::shared_ptr<DispersionOperation> lops,
                 std::shared_ptr<DispersionOperation> rops)
      : _lbranch_clv_index{lchild_clv},
        _rbranch_clv_index{rchild_clv},
        _parent_clv_index{parent_clv},
        _lbranch_ops{{lops}},
        _rbranch_ops{{rops}} {}

  SplitOperation(size_t lchild_clv_top, size_t rchild_clv_top,
                 size_t lchild_clv_bot, size_t rchild_clv_bot, double lbrlen,
                 double rbrlen, size_t lprob_mat, size_t rprob_mat,
                 size_t lrate_matrix, size_t rrate_matrix, size_t parent_clv)
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
                 double rbrlen, size_t prob_mat, size_t rate_matrix,
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

 private:
  size_t _lbranch_clv_index;
  size_t _rbranch_clv_index;
  size_t _parent_clv_index;

  std::vector<std::shared_ptr<DispersionOperation>> _lbranch_ops;
  std::vector<std::shared_ptr<DispersionOperation>> _rbranch_ops;
};

class ReverseSplitOperation {
 public:
  ReverseSplitOperation(
      size_t bot_clv,  /* Where the result is stored*/
      size_t ltop_clv, /* Where the result of the dispop is stored */
      size_t rtop_clv, /* Should be a _non_ reverse CLV */
      size_t rate_matrix_index, size_t prob_matrix_index, size_t disp_clv_index,
      double brlen)
      : _bot_clv_index{bot_clv},
        _ltop_clv_index{ltop_clv},
        _rtop_clv_index{rtop_clv},
        _eval_clvs{true},
        _branch_ops{{std::make_shared<DispersionOperation>(
            ltop_clv, disp_clv_index, brlen, prob_matrix_index,
            rate_matrix_index)}} {}

  ReverseSplitOperation(size_t bot_clv, size_t rtop_clv,
                        std::shared_ptr<DispersionOperation> branch_op)
      : _bot_clv_index{bot_clv},
        _ltop_clv_index{branch_op->top_clv_index()},
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

 private:
  size_t _bot_clv_index;
  size_t _ltop_clv_index;
  size_t _rtop_clv_index;
  bool _eval_clvs;

  std::vector<std::shared_ptr<DispersionOperation>> _branch_ops;
};

class LHGoal {
 public:
  LHGoal(size_t root_clv, size_t prior_index)
      : _root_clv_index{root_clv}, _prior_index{prior_index} {}
  double eval(std::shared_ptr<Workspace>) const;

 private:
  size_t _root_clv_index;
  size_t _prior_index;
};

class StateLHGoal {
 public:
  StateLHGoal(size_t parent_clv, size_t lchild_clv, size_t rchild_clv)
      : _parent_clv_index{parent_clv},
        _lchild_clv_index{lchild_clv},
        _rchild_clv_index{rchild_clv} {}

  void eval(std::shared_ptr<Workspace>) const;

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;
};

class SplitLHGoal {
 public:
  void eval(std::shared_ptr<Workspace>) const;

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;
};

struct OperationWrapper {
  std::vector<LHGoal> _lh;
  std::vector<StateLHGoal> _state;
  std::vector<SplitLHGoal> _split;
  std::vector<SplitOperation> _forward_ops;
  std::vector<ReverseSplitOperation> _backwards_ops;

  void registerLHGoal() {
    auto root_clv = (_forward_ops.end() - 1)->get_parent_clv();
    size_t frequency_index = 0;
    _lh.emplace_back(root_clv, frequency_index);
  }

  void registerStateLHGoal(size_t parent_clv, size_t c1_clv, size_t c2_clv) {
    _state.emplace_back(parent_clv, c1_clv, c2_clv);
  }
};

#endif
