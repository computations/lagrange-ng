/* Operation.h
 *
 * Created: 27 Oct 2020
 * Author: Ben Bettisworth
 */

#ifndef LAGRANGE_OPERATION_H
#define LAGRANGE_OPERATION_H

#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "AncSplit.h"
#include "Common.h"
#include "Workspace.h"

class MakeRateMatrixOperation {
 public:
  explicit MakeRateMatrixOperation(size_t index)
      : _rate_matrix_index{index},
        _period_index{0},
        _last_execution{0},
        _last_update{0} {}

  void eval(const std::shared_ptr<Workspace>& ws);

  inline auto last_update() const -> lagrange_clock_tick_t {
    return _last_execution;
  }

  inline void update_rates(const std::shared_ptr<Workspace>& ws, double disp,
                           double ext) {
    ws->set_period_params(_period_index, disp, ext);
    _last_update = ws->advance_clock();
  }

  inline void update_rates(const std::shared_ptr<Workspace>& ws,
                           const period_t& p) {
    update_rates(ws, p.dispersion_rate, p.extinction_rate);
  }

  inline auto rate_matrix_index() const -> size_t { return _rate_matrix_index; }

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  auto printStatus(const std::shared_ptr<Workspace>& ws,
                   size_t tabLevel = 0) const -> std::string;

 private:
  size_t _rate_matrix_index;
  size_t _period_index;

  std::unique_ptr<std::mutex> _lock{new std::mutex};

  lagrange_clock_tick_t _last_execution;
  lagrange_clock_tick_t _last_update;
};

class ExpmOperation {
 public:
  ExpmOperation(size_t prob_matrix, double t,
                const std::shared_ptr<MakeRateMatrixOperation>& rm_op,
                bool transpose)
      : _prob_matrix_index{prob_matrix},
        _rate_matrix_index{rm_op->rate_matrix_index()},
        _t{t},
        _rate_matrix_op{rm_op},
        _transposed{transpose} {}

  ExpmOperation(size_t prob_matrix, double t,
                const std::shared_ptr<MakeRateMatrixOperation>& rm_op)
      : ExpmOperation(prob_matrix, t, rm_op, false) {}

  ExpmOperation(const ExpmOperation&) = delete;
  auto operator=(const ExpmOperation&) -> ExpmOperation& = delete;

  void eval(const std::shared_ptr<Workspace>& ws);

  auto prob_matrix() const -> size_t { return _prob_matrix_index; }

  auto last_execution() const -> lagrange_clock_tick_t {
    return _last_execution;
  }

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  auto printStatus(const std::shared_ptr<Workspace>& ws,
                   size_t tabLevel = 0) const -> std::string;

  auto getLock() -> std::mutex& { return *_lock; }

 private:
  size_t _prob_matrix_index;
  size_t _rate_matrix_index;
  double _t;

  std::shared_ptr<MakeRateMatrixOperation> _rate_matrix_op;

  /* Temp matrices for the computation of the exponential */
  std::unique_ptr<lagrange_matrix_base_t[]> _A = nullptr;
  std::unique_ptr<lagrange_matrix_base_t[]> _X_1 = nullptr;
  std::unique_ptr<lagrange_matrix_base_t[]> _X_2 = nullptr;
  std::unique_ptr<lagrange_matrix_base_t[]> _N = nullptr;
  std::unique_ptr<lagrange_matrix_base_t[]> _D = nullptr;
  std::unique_ptr<lagrange_matrix_base_t[]> _lapack_work_buffer = nullptr;

  lagrange_clock_tick_t _last_execution = 0;
  std::unique_ptr<std::mutex> _lock{new std::mutex};

  bool _transposed;
  size_t _execution_count = 0;
};

class DispersionOperation {
 public:
  DispersionOperation(size_t top, size_t bot, double brlen, size_t prob_matrix,
                      const std::shared_ptr<MakeRateMatrixOperation>& rm_op)
      : _top_clv{top},
        _bot_clv{bot},
        _prob_matrix_index{prob_matrix},
        _expm_op{new ExpmOperation{prob_matrix, brlen, rm_op}} {}

  DispersionOperation(size_t top, size_t bot, double brlen, size_t prob_matrix,
                      const std::shared_ptr<MakeRateMatrixOperation>& rm_op,
                      bool transpose)
      : _top_clv{top},
        _bot_clv{bot},
        _prob_matrix_index{prob_matrix},
        _expm_op{new ExpmOperation{prob_matrix, brlen, rm_op, transpose}} {}

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

  DispersionOperation(size_t bot, const std::shared_ptr<ExpmOperation>& op)
      : _top_clv{std::numeric_limits<size_t>::max()},
        _bot_clv{bot},
        _prob_matrix_index{op->prob_matrix()},
        _expm_op{op} {}

  DispersionOperation(size_t top, size_t bot, size_t prob_matrix)
      : _top_clv{top},
        _bot_clv{bot},
        _prob_matrix_index{prob_matrix},
        _expm_op{nullptr} {}

  void eval(const std::shared_ptr<Workspace>& ws);

  void terminate_top(size_t clv_index) {
    if (_top_clv != std::numeric_limits<size_t>::max()) {
      throw std::runtime_error{
          "Dispersion Operation already terminated at the top"};
    }
    _top_clv = clv_index;
  }

  void terminate_bot(size_t clv_index) {
    if (_bot_clv != std::numeric_limits<size_t>::max()) {
      throw std::runtime_error{
          "Dispersion Operation already terminated at the bot"};
    }
    _bot_clv = clv_index;
  }

  auto top_clv_index() const -> size_t { return _top_clv; }
  auto bot_clv_index() const -> size_t { return _bot_clv; }

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  auto printStatus(const std::shared_ptr<Workspace>& ws,
                   size_t tabLevel = 0) const -> std::string;

  auto ready(const std::shared_ptr<Workspace>& ws,
             lagrange_clock_tick_t deadline) const -> bool {
    return ws->last_update_clv(_bot_clv) > deadline;
  }

  auto getLock() -> std::mutex& { return *_lock; }

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

  std::unique_ptr<std::mutex> _lock{new std::mutex};
  lagrange_clock_tick_t _last_execution = 0;
};

class SplitOperation {
 public:
  SplitOperation(size_t lchild_clv_top, size_t rchild_clv_top,
                 size_t lchild_clv_bot, size_t rchild_clv_bot, double lbrlen,
                 double rbrlen, size_t lprob_mat, size_t rprob_mat,
                 const std::shared_ptr<MakeRateMatrixOperation>& lrate_matrix,
                 const std::shared_ptr<MakeRateMatrixOperation>& rrate_matrix,
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
                 const std::shared_ptr<MakeRateMatrixOperation>& rate_matrix,
                 size_t parent_clv)
      : SplitOperation(lchild_clv_top, rchild_clv_top, lchild_clv_bot,
                       rchild_clv_bot, lbrlen, rbrlen, prob_mat, prob_mat,
                       rate_matrix, rate_matrix, parent_clv) {}

  // TODO: Fix this for periods
  SplitOperation(size_t parent_clv, std::shared_ptr<DispersionOperation> l_ops,
                 std::shared_ptr<DispersionOperation> r_ops)
      : _lbranch_clv_index{l_ops->top_clv_index()},
        _rbranch_clv_index{r_ops->top_clv_index()},
        _parent_clv_index{parent_clv},
        _lbranch_ops{{l_ops}},
        _rbranch_ops{{r_ops}} {}

  void eval(const std::shared_ptr<Workspace>&);

  auto get_parent_clv() const -> size_t { return _parent_clv_index; }

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  auto printStatus(const std::shared_ptr<Workspace>& ws,
                   size_t tabLevel = 0) const -> std::string;

  auto ready(const std::shared_ptr<Workspace>& ws) const -> bool {
    return _lbranch_ops[_lbranch_ops.size() - 1]->ready(ws, _last_execution) &&
           _rbranch_ops[_rbranch_ops.size() - 1]->ready(ws, _last_execution);
  }

  auto getLock() -> std::mutex& { return *_lock; }

 private:
  size_t _lbranch_clv_index;
  size_t _rbranch_clv_index;
  size_t _parent_clv_index;

  std::vector<std::shared_ptr<DispersionOperation>> _lbranch_ops;
  std::vector<std::shared_ptr<DispersionOperation>> _rbranch_ops;

  std::unique_ptr<std::mutex> _lock{new std::mutex};
  lagrange_clock_tick_t _last_execution = 0;
};

class ReverseSplitOperation {
 public:
  ReverseSplitOperation(
      size_t bot_clv,  /* Where the result is stored*/
      size_t ltop_clv, /* Where the result of the dispop is stored */
      size_t rtop_clv, /* Should be a _non_ reverse CLV */
      const std::shared_ptr<MakeRateMatrixOperation>& rate_matrix_op,
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
        _ltop_clv_index{branch_op->top_clv_index()},
        _rtop_clv_index{rtop_clv},
        _eval_clvs{true},
        _branch_ops{{branch_op}} {}

  explicit ReverseSplitOperation(std::shared_ptr<DispersionOperation> branch_op)
      : _bot_clv_index{std::numeric_limits<size_t>::max()},
        _ltop_clv_index{std::numeric_limits<size_t>::max()},
        _rtop_clv_index{std::numeric_limits<size_t>::max()},
        _eval_clvs{false},
        _branch_ops{{branch_op}} {}

  void eval(const std::shared_ptr<Workspace>&);

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  auto printStatus(const std::shared_ptr<Workspace>& ws,
                   size_t tabLevel = 0) const -> std::string;

  void makeRootOperation(size_t clv_index) {
    _branch_ops.clear();
    if (_ltop_clv_index == std::numeric_limits<size_t>::max()) {
      _ltop_clv_index = clv_index;
    }
    if (_rtop_clv_index == std::numeric_limits<size_t>::max()) {
      _rtop_clv_index = clv_index;
    }
  }

  auto getStableCLV() const -> size_t { return _ltop_clv_index; }

  auto ready(const std::shared_ptr<Workspace>& ws) const -> bool {
    bool branch_ops_ready = true;
    for (const auto& op : _branch_ops) {
      branch_ops_ready = branch_ops_ready && op->ready(ws, _last_execution);
    }

    return branch_ops_ready &&
           (ws->last_update_clv(_ltop_clv_index) >= _last_execution);
  }

  auto getLock() -> std::mutex& { return *_lock; }

 private:
  size_t _bot_clv_index;
  size_t _ltop_clv_index;
  size_t _rtop_clv_index;
  bool _eval_clvs;

  std::vector<std::shared_ptr<DispersionOperation>> _branch_ops;
  std::vector<lagrange_dist_t> _excl_dists;

  std::unique_ptr<std::mutex> _lock{new std::mutex};
  lagrange_clock_tick_t _last_execution = 0;
};

class LLHGoal {
 public:
  LLHGoal(size_t root_clv, size_t prior_index)
      : _root_clv_index{root_clv}, _prior_index{prior_index} {}
  void eval(const std::shared_ptr<Workspace>&);

  inline auto result() const -> double { return _result; }

  size_t _root_clv_index;
  size_t _prior_index;
  double _result = -std::numeric_limits<double>::infinity();

  auto ready(const std::shared_ptr<Workspace>&) const -> bool;

 private:
  lagrange_clock_tick_t _last_execution = 0;
};

class StateLHGoal {
 public:
  StateLHGoal(size_t parent_clv, size_t lchild_clv, size_t rchild_clv)
      : _parent_clv_index{parent_clv},
        _lchild_clv_index{lchild_clv},
        _rchild_clv_index{rchild_clv},
        _result{nullptr},
        _states{0} {}
  StateLHGoal(const StateLHGoal&) = delete;
  auto operator=(const StateLHGoal&) -> StateLHGoal& = delete;
  StateLHGoal(StateLHGoal&& other) noexcept
      : _parent_clv_index{other._parent_clv_index},
        _lchild_clv_index{other._lchild_clv_index},
        _rchild_clv_index{other._rchild_clv_index},
        _result{std::move(other._result)},
        _states{other._states} {}
  auto operator=(StateLHGoal&&) -> StateLHGoal& = delete;

  void eval(const std::shared_ptr<Workspace>&);

  inline auto result() {
    decltype(_result) ret{new decltype(_result)::element_type[_states]};
    for (size_t i = 0; i < _states; i++) { ret.get()[i] = _result.get()[i]; }
    return ret;
  }

  auto ready(const std::shared_ptr<Workspace>&) const -> bool;

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;

  lagrange_clock_tick_t _last_execution = 0;

  std::unique_ptr<lagrange_matrix_base_t[]> _result;
  size_t _states;
};

class SplitLHGoal {
 public:
  SplitLHGoal(size_t parent_clv, size_t lchild_clv, size_t rchild_clv)
      : _parent_clv_index{parent_clv},
        _lchild_clv_index{lchild_clv},
        _rchild_clv_index{rchild_clv} {}

  void eval(const std::shared_ptr<Workspace>&);

  inline auto result() const -> lagrange_split_return_t { return _result; }

  auto ready(const std::shared_ptr<Workspace>&) const -> bool;

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;

  lagrange_clock_tick_t _last_execution = 0;

  lagrange_split_return_t _result;
};

using PeriodRateMatrixMap =
    std::unordered_map<size_t, std::shared_ptr<MakeRateMatrixOperation>>;

using BranchProbMatrixMap = std::unordered_map<std::pair<size_t, double>,
                                               std::shared_ptr<ExpmOperation>>;

#endif
