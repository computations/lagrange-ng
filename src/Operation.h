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
  MakeRateMatrixOperation(size_t index)
      : _rate_matrix_index{index}, _period_index{0} {}

  void eval(const std::shared_ptr<Workspace>& ws);

  inline lagrange_clock_tick_t last_update() const { return _last_execution; }

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

  lagrange_clock_tick_t _last_execution = 0;
  lagrange_clock_tick_t _last_update = 0;
};

class ExpmOperation {
 public:
  ExpmOperation(size_t prob_matrix, size_t rate_matrix_index, double t,
                bool transpose)
      : _prob_matrix_index{prob_matrix},
        _rate_matrix_index{rate_matrix_index},
        _t{t},
        _transposed{transpose} {}

  ExpmOperation(const ExpmOperation& op) {
    _prob_matrix_index = op._prob_matrix_index;
    _rate_matrix_index = op._rate_matrix_index;
    _t = op._t;
    _transposed = op._transposed;
  }

  ExpmOperation& operator=(const ExpmOperation&) = delete;

  void eval(const std::shared_ptr<Workspace>& ws);

  size_t prob_matrix() const { return _prob_matrix_index; }

  lagrange_clock_tick_t last_execution() const { return _last_execution; }

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  size_t get_index() const { return _prob_matrix_index; }

  std::string printStatus(const std::shared_ptr<Workspace>& ws,
                          size_t tabLevel = 0) const;

  std::mutex& getLock() { return *_lock; }

 private:
  size_t _prob_matrix_index;
  size_t _rate_matrix_index;
  double _t;

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
};

class DispersionOperation {
 public:
  DispersionOperation(size_t top, size_t bot, size_t prob_matrix)
      : _top_clv{top}, _bot_clv{bot}, _prob_matrix_index{prob_matrix} {}

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

  bool is_bot_terminated() const {
    return _bot_clv == std::numeric_limits<size_t>::max();
  }

  size_t top_clv_index() const { return _top_clv; };
  size_t bot_clv_index() const { return _bot_clv; };

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  std::string printStatus(const std::shared_ptr<Workspace>& ws,
                          size_t tabLevel = 0) const;

  bool ready(const std::shared_ptr<Workspace>& ws,
             lagrange_clock_tick_t deadline) const {
    return ws->last_update_clv(_bot_clv) > deadline;
  }

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

  lagrange_clock_tick_t _last_execution = 0;
};

class SplitOperation {
 public:
  SplitOperation(size_t lchild_clv_top, size_t rchild_clv_top,
                 size_t parent_clv)
      : _lbranch_clv_index{lchild_clv_top},
        _rbranch_clv_index{rchild_clv_top},
        _parent_clv_index{parent_clv} {}

  void eval(const std::shared_ptr<Workspace>&);

  size_t get_parent_clv() const { return _parent_clv_index; }

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  std::string printStatus(const std::shared_ptr<Workspace>& ws,
                          size_t tabLevel = 0) const;

  bool ready(const std::shared_ptr<Workspace>& ws) const {
    throw std::runtime_error{"Not implemented"};
    return false;
  }

 private:
  size_t _lbranch_clv_index;
  size_t _rbranch_clv_index;
  size_t _parent_clv_index;

  lagrange_clock_tick_t _last_execution = 0;
};

class Operation {
 public:
  Operation(const SplitOperation& sp_op,
            const std::vector<DispersionOperation>& left_disp_ops,
            const std::vector<DispersionOperation>& right_disp_ops,
            const std::vector<ExpmOperation>& expm_ops)
      : _split_op{sp_op},
        _lbranch_ops{left_disp_ops},
        _rbranch_ops{right_disp_ops},
        _expm_ops{expm_ops} {}

  void eval(const std::shared_ptr<Workspace>& workspace);
  bool ready(const std::shared_ptr<Workspace>& workspace) { return true; }

  void set_split_operation(const SplitOperation& split_op) {
    _split_op = split_op;
  }

  void add_lbranch_op(const DispersionOperation& disp_op) {
    _lbranch_ops.push_back(disp_op);
  }

  void add_rbranch_op(const DispersionOperation& disp_op) {
    _lbranch_ops.push_back(disp_op);
  }

  void add_expm_op(const ExpmOperation& expm_op) {
    _expm_ops.push_back(expm_op);
  }

  size_t get_parent_clv() const { return _split_op.get_parent_clv(); }

 private:
  SplitOperation _split_op;

  std::vector<DispersionOperation> _lbranch_ops;
  std::vector<DispersionOperation> _rbranch_ops;

  std::vector<ExpmOperation> _expm_ops;

  lagrange_clock_tick_t _last_execution;
};

class ReverseSplitOperation {
 public:
  ReverseSplitOperation(
      size_t bot_clv,  /* Where the result is stored*/
      size_t ltop_clv, /* Where the result of the dispop is stored */
      size_t rtop_clv  /* Should be a _non_ reverse CLV */
      )
      : _bot_clv_index{bot_clv},
        _ltop_clv_index{ltop_clv},
        _rtop_clv_index{rtop_clv},
        _eval_clvs{true} {}

  void eval(const std::shared_ptr<Workspace>&);

  void printStatus(const std::shared_ptr<Workspace>& ws, std::ostream& os,
                   size_t tabLevel = 0) const;

  std::string printStatus(const std::shared_ptr<Workspace>& ws,
                          size_t tabLevel = 0) const;

  size_t getStableCLV() const { return _ltop_clv_index; }

  bool ready(const std::shared_ptr<Workspace>& ws) const {
    throw std::runtime_error{
        "Ready is Not Implementd for ReverseSplitOperation"};
  }

  void makeRootOperation(size_t clv_index) {
    if (_ltop_clv_index == std::numeric_limits<size_t>::max()) {
      _ltop_clv_index = clv_index;
    } else if (_rtop_clv_index == std::numeric_limits<size_t>::max()) {
      _rtop_clv_index = clv_index;
    }
  }

  size_t getLTopCLVIndex() const { return _ltop_clv_index; }
  size_t getRTopCLVIndex() const { return _rtop_clv_index; }
  size_t getBotCLVIndex() const { return _bot_clv_index; }

  void setLTopCLVIndex(size_t index) { _ltop_clv_index = index; }

  void setRTopCLVIndex(size_t index) { _rtop_clv_index = index; }

 private:
  size_t _bot_clv_index;
  size_t _ltop_clv_index;
  size_t _rtop_clv_index;
  bool _eval_clvs;

  lagrange_clock_tick_t _last_execution = 0;
};

class ReverseOperation {
 public:
  ReverseOperation(ReverseSplitOperation rs_op,
                   std::vector<DispersionOperation> disp_ops,
                   std::vector<ExpmOperation> expm_ops)
      : _rsplit_op{rs_op}, _branch_ops{disp_ops}, _expm_ops{expm_ops} {}

  void eval(const std::shared_ptr<Workspace>& ws);

  void makeRootOperation(size_t clv_index) {
    _branch_ops.clear();
    _rsplit_op.makeRootOperation(clv_index);
  }

  void terminate_bot(size_t clv_index) {
    _branch_ops.back().terminate_bot(clv_index);
  }

  bool is_bot_terminated() const {
    return _branch_ops.back().is_bot_terminated();
  }

  size_t getStableCLV() const { return _rsplit_op.getStableCLV(); }

  bool ready(const std::shared_ptr<Workspace>& ws) const { return true; }

 private:
  ReverseSplitOperation _rsplit_op;
  std::vector<DispersionOperation> _branch_ops;

  std::vector<ExpmOperation> _expm_ops;

  lagrange_clock_tick_t _last_execution;
};

class LLHGoal {
 public:
  LLHGoal(size_t root_clv, size_t prior_index)
      : _root_clv_index{root_clv}, _prior_index{prior_index} {}
  void eval(const std::shared_ptr<Workspace>&);

  inline double result() const { return _result; }

  size_t _root_clv_index;
  size_t _prior_index;
  double _result = -std::numeric_limits<double>::infinity();

  bool ready(const std::shared_ptr<Workspace>&) const;

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
  StateLHGoal& operator=(const StateLHGoal&) = delete;
  StateLHGoal(StateLHGoal&& other)
      : _parent_clv_index{other._parent_clv_index},
        _lchild_clv_index{other._lchild_clv_index},
        _rchild_clv_index{other._rchild_clv_index},
        _result{std::move(other._result)} {}
  StateLHGoal& operator=(StateLHGoal&&) = delete;

  void eval(const std::shared_ptr<Workspace>&);

  inline auto result() {
    decltype(_result) ret{new decltype(_result)::element_type[_states]};
    for (size_t i = 0; i < _states; i++) { ret.get()[i] = _result.get()[i]; }
    return ret;
  }

  bool ready(const std::shared_ptr<Workspace>&) const;

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
  void eval(const std::shared_ptr<Workspace>&);

  inline lagrange_split_return_t result() const { return _result; }

  bool ready(const std::shared_ptr<Workspace>&) const;

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;

  lagrange_clock_tick_t _last_execution = 0;

  lagrange_split_return_t _result;
};

typedef std::unordered_map<size_t, MakeRateMatrixOperation> PeriodRateMatrixMap;

typedef std::unordered_map<std::pair<size_t, double>, ExpmOperation>
    BranchProbMatrixMap;

#endif
