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
#include <stdexcept>
#include <unordered_map>

#include "AncSplit.hpp"
#include "Common.hpp"
#include "Utils.hpp"
#include "Workspace.hpp"
#include "logger.hpp"

namespace lagrange {

class MakeRateMatrixOperation {
 public:
  explicit MakeRateMatrixOperation(size_t index, size_t period_index) :
      _rate_matrix_index{index},
      _period_index{period_index} {}

  void eval(const std::shared_ptr<Workspace>& ws);

  [[nodiscard]] auto lastUpdate() const -> ClockTick { return _last_execution; }

  void updateRates(const std::shared_ptr<Workspace>& ws,
                   double disp,
                   double ext) {
    ws->setPeriodParams(_period_index, disp, ext);
    _last_update = ws->advanceClock();
  }

  void updateRates(const std::shared_ptr<Workspace>& ws,
                   const PeriodParams& p) {
    updateRates(ws, p.dispersion_rate, p.extinction_rate);
  }

  [[nodiscard]] auto rateMatrixIndex() const -> size_t {
    return _rate_matrix_index;
  }

  void printStatus(const std::shared_ptr<Workspace>& ws,
                   std::ostream& os,
                   size_t tabLevel = 0) const;

  [[nodiscard]] auto printStatus(const std::shared_ptr<Workspace>& ws,
                                 size_t tabLevel = 0) const -> std::string;

  [[nodiscard]] auto getPeriodIndex() const -> size_t { return _period_index; }

  auto evaluated() const -> bool { return _last_execution > _last_update; }

 private:
  size_t _rate_matrix_index;
  size_t _period_index;

  std::unique_ptr<std::mutex> _lock{new std::mutex};

  ClockTick _last_execution{0};
  ClockTick _last_update{0};
};

class ExpmOperation {
 public:
  ExpmOperation(size_t prob_matrix,
                double t,
                const std::shared_ptr<MakeRateMatrixOperation>& rm_op,
                bool transpose) :
      _prob_matrix_index{prob_matrix},
      _rate_matrix_index{rm_op->rateMatrixIndex()},
      _t{t},
      _rate_matrix_op{rm_op},
      _transposed{transpose} {}

  ExpmOperation(size_t prob_matrix,
                double t,
                const std::shared_ptr<MakeRateMatrixOperation>& rm_op) :
      ExpmOperation(prob_matrix, t, rm_op, false) {}

  ExpmOperation(const ExpmOperation&) = delete;
  auto operator=(const ExpmOperation&) -> ExpmOperation& = delete;

  void eval(const std::shared_ptr<Workspace>& ws);

  [[nodiscard]] auto probMatrix() const -> size_t { return _prob_matrix_index; }

  [[nodiscard]] auto rateMatrix() const -> size_t { return _rate_matrix_index; }

  [[nodiscard]] auto transposed() const -> bool { return _transposed; }

  [[nodiscard]] auto getT() const { return _t; }

  [[nodiscard]] auto lastExecution() const -> ClockTick {
    return _last_execution;
  }

  void printStatus(const std::shared_ptr<Workspace>& ws,
                   std::ostream& os,
                   size_t tabLevel = 0) const;

  [[nodiscard]] auto printStatus(const std::shared_ptr<Workspace>& ws,
                                 size_t tabLevel = 0) const -> std::string;

  auto getLock() -> std::mutex& { return *_lock; }

  void setArnoldiMode(bool d) { _arnoldi_mode = d; }

  [[nodiscard]] auto isArnoldiMode() const -> bool { return _arnoldi_mode; }

  void setAdaptive(bool m) { _adaptive_mode = m; }

  [[nodiscard]] auto isAdaptive() const -> bool { return _adaptive_mode; }

  void emplaceRateMatrixOperations(
      std::unordered_map<size_t, std::shared_ptr<MakeRateMatrixOperation>>&
          vec) {
    if (!_rate_matrix_op) { return; }

    size_t period = _rate_matrix_op->getPeriodIndex();
    if (!vec.contains(period)) { vec[period] = _rate_matrix_op; }
  }

  auto evaluated(const std::shared_ptr<Workspace>& ws) const -> bool {
    auto rm_eval = _rate_matrix_op ? _rate_matrix_op->evaluated() : true;
    auto pm_eval =
        _last_execution > ws->lastUpdateProbMatrix(_rate_matrix_index);

    return rm_eval && pm_eval;
  }

 private:
  size_t _prob_matrix_index;
  size_t _rate_matrix_index;
  double _t;

  std::shared_ptr<MakeRateMatrixOperation> _rate_matrix_op;

  /* Temp matrices for the computation of the exponential */
  std::unique_ptr<LagrangeMatrixBase[]> _A = nullptr;
  std::unique_ptr<LagrangeMatrixBase[]> _X_1 = nullptr;
  std::unique_ptr<LagrangeMatrixBase[]> _X_2 = nullptr;
  std::unique_ptr<LagrangeMatrixBase[]> _N = nullptr;
  std::unique_ptr<LagrangeMatrixBase[]> _D = nullptr;
  std::unique_ptr<LagrangeMatrixBase[]> _lapack_work_buffer = nullptr;

  ClockTick _last_execution = 0;
  std::unique_ptr<std::mutex> _lock{new std::mutex};

  bool _transposed;
  bool _arnoldi_mode = false;
  bool _adaptive_mode = false;
  size_t _execution_count = 0;
};

class DispersionOperation {
 public:
  DispersionOperation(size_t top,
                      size_t bot,
                      double brlen,
                      size_t prob_matrix,
                      const std::shared_ptr<MakeRateMatrixOperation>& rm_op) :
      _top_clv{top},
      _bot_clv{bot},
      _prob_matrix_index{prob_matrix},
      _expm_op{new ExpmOperation{prob_matrix, brlen, rm_op}} {}

  DispersionOperation(size_t top,
                      size_t bot,
                      double brlen,
                      size_t prob_matrix,
                      const std::shared_ptr<MakeRateMatrixOperation>& rm_op,
                      bool transpose) :
      _top_clv{top},
      _bot_clv{bot},
      _prob_matrix_index{prob_matrix},
      _expm_op{new ExpmOperation{prob_matrix, brlen, rm_op, transpose}} {}

  DispersionOperation(size_t bot,
                      double brlen,
                      size_t prob_matrix,
                      const std::shared_ptr<MakeRateMatrixOperation>& rm_op) :
      _top_clv{std::numeric_limits<size_t>::max()},
      _bot_clv{bot},
      _prob_matrix_index{prob_matrix},
      _expm_op{new ExpmOperation{prob_matrix, brlen, rm_op}} {}

  DispersionOperation(size_t bot,
                      double brlen,
                      size_t prob_matrix,
                      const std::shared_ptr<MakeRateMatrixOperation>& rm_op,
                      bool transpose) :
      _top_clv{std::numeric_limits<size_t>::max()},
      _bot_clv{bot},
      _prob_matrix_index{prob_matrix},
      _expm_op{new ExpmOperation{prob_matrix, brlen, rm_op, transpose}} {}

  DispersionOperation(size_t top,
                      size_t bot,
                      const std::shared_ptr<ExpmOperation>& op) :
      _top_clv{top},
      _bot_clv{bot},
      _prob_matrix_index{op->probMatrix()},
      _expm_op{op} {}

  DispersionOperation(size_t bot, const std::shared_ptr<ExpmOperation>& op) :
      _top_clv{std::numeric_limits<size_t>::max()},
      _bot_clv{bot},
      _prob_matrix_index{op->probMatrix()},
      _expm_op{op} {}

  DispersionOperation(size_t top, size_t bot, size_t prob_matrix) :
      _top_clv{top},
      _bot_clv{bot},
      _prob_matrix_index{prob_matrix},
      _expm_op{nullptr} {}

  void eval(const std::shared_ptr<Workspace>& ws);

  void terminateTop(size_t clv_index) {
    if (_top_clv != std::numeric_limits<size_t>::max()) {
      throw std::runtime_error{
          "Dispersion Operation already terminated at the top"};
    }
    _top_clv = clv_index;
  }

  void terminateBotChildren(size_t clv_index) {
    if (_child_op == nullptr) {
      terminateBot(clv_index);
    } else {
      _child_op->terminateBotChildren(clv_index);
    }
  }

  void unterminateTop() {
    if (_top_clv == std::numeric_limits<size_t>::max()) {
      throw std::runtime_error{"Tried to unterminate an unterminated top"};
    }
    _top_clv = std::numeric_limits<size_t>::max();
  }

  void terminateBot(size_t clv_index) {
    if (_bot_clv != std::numeric_limits<size_t>::max()) {
      throw std::runtime_error{
          "Dispersion Operation already terminated at the bot"};
    }
    _bot_clv = clv_index;
  }

  void unterminateBot() {
    if (_bot_clv == std::numeric_limits<size_t>::max()) {
      throw std::runtime_error{"Tried to unterminate an unterminated bot"};
    }
    _bot_clv = std::numeric_limits<size_t>::max();
  }

  [[nodiscard]] auto topCLVIndex() const -> size_t { return _top_clv; }

  [[nodiscard]] auto botCLVIndex() const -> size_t { return _bot_clv; }

  void printStatus(const std::shared_ptr<Workspace>& ws,
                   std::ostream& os,
                   size_t tabLevel = 0) const;

  [[nodiscard]] auto printStatus(const std::shared_ptr<Workspace>& ws,
                                 size_t tabLevel = 0) const -> std::string;

  [[nodiscard]] auto ready(const std::shared_ptr<Workspace>& ws,
                           ClockTick deadline) const -> bool {
    if (ws->lastUpdateCLV(_bot_clv) > deadline) { return true; }
    return (_child_op != nullptr && _child_op->ready(ws, deadline));
  }

  auto getLock() -> std::mutex& { return _lock; }

  [[nodiscard]] auto getExpmOperation() const
      -> std::shared_ptr<ExpmOperation> {
    return _expm_op;
  }

  void fallback() {
    LOG_WARNING("Falling back to Pade exponentiation for probability matrix {}",
                _prob_matrix_index);
    _expm_op->setArnoldiMode(false);
  }

  void setChildOp(const std::shared_ptr<DispersionOperation>& op) {
    if (_child_op != nullptr) {
      throw std::runtime_error{"Child is already set"};
    }
    _child_op = op;
  }

  void setChildOpAndTerminate(std::shared_ptr<DispersionOperation>& op) {
    setChildOp(op);
    terminateBot(op->topCLVIndex());
  }

  auto getExpmOperations() -> std::vector<std::shared_ptr<ExpmOperation>> {
    size_t op_count = countChildOps();
    std::vector<std::shared_ptr<ExpmOperation>> expm_ops;
    expm_ops.reserve(op_count);
    insertExpmOpRecursive(expm_ops);
    return expm_ops;
  }

  auto countChildOps() -> size_t {
    if (!_child_op) { return 1; }
    return _child_op->countChildOps() + 1;
  }

  void emplaceRateMatrixOperations(
      std::unordered_map<size_t, std::shared_ptr<MakeRateMatrixOperation>>&
          vec) {
    if (_child_op) { _child_op->emplaceRateMatrixOperations(vec); }
    _expm_op->emplaceRateMatrixOperations(vec);
  }

  void printGraph(std::ostream& os, size_t& index) const;

  auto evaluated(const std::shared_ptr<Workspace>& ws) const -> bool {
    auto eval = ws->lastUpdateCLV(_bot_clv) < _last_execution
                && _last_execution < ws->lastUpdateCLV(_top_clv)
                && _expm_op->evaluated(ws);
    eval &= _child_op ? _last_execution > _child_op->evaluated(ws) : true;

    return eval;
  }

  auto getFinalDest() const {
    if (_child_op == nullptr) { return _top_clv; }
    return _child_op->getFinalDest();
  }

 private:
  void insertExpmOpRecursive(
      std::vector<std::shared_ptr<ExpmOperation>>& expm_ops) const {
    if (_child_op != nullptr) { _child_op->insertExpmOpRecursive(expm_ops); }
    expm_ops.push_back(_expm_op);
  }

  /* Remember, the top and bottom clv indexes could be the same. This is to save
   * on storage when computing different periods along a single branch. The idea
   * is that we just apply the matrix multiplication over and over again to the
   * same CLV. Matvecs can't be done in place, so there will need to be
   * temporary storage if we do that.
   */
  size_t _top_clv;
  size_t _bot_clv;
  size_t _prob_matrix_index;

  /* _expm_op is an non-owning pointer, and can be null */
  std::shared_ptr<ExpmOperation> _expm_op;
  std::shared_ptr<DispersionOperation> _child_op;

  std::mutex _lock;
  ClockTick _last_execution = 0;
};

class SplitOperation {
 public:
  SplitOperation(size_t lchild_clv_top,
                 size_t rchild_clv_top,
                 size_t lchild_clv_bot,
                 size_t rchild_clv_bot,
                 double lbrlen,
                 double rbrlen,
                 size_t lprob_mat,
                 size_t rprob_mat,
                 const std::shared_ptr<MakeRateMatrixOperation>& lrate_matrix,
                 const std::shared_ptr<MakeRateMatrixOperation>& rrate_matrix,
                 size_t parent_clv) :
      _lbranch_clv_index{lchild_clv_top},
      _rbranch_clv_index{rchild_clv_top},
      _parent_clv_index{parent_clv},
      _lbranch_op{{std::make_shared<DispersionOperation>(
          lchild_clv_top, lchild_clv_bot, lbrlen, lprob_mat, lrate_matrix)}},
      _rbranch_op{{std::make_shared<DispersionOperation>(
          rchild_clv_top, rchild_clv_bot, rbrlen, rprob_mat, rrate_matrix)}} {}

  SplitOperation(size_t lchild_clv_top,
                 size_t lchild_clv_bot,
                 size_t rchild_clv_top,
                 size_t rchild_clv_bot,
                 double lbrlen,
                 double rbrlen,
                 size_t prob_mat,
                 const std::shared_ptr<MakeRateMatrixOperation>& rate_matrix,
                 size_t parent_clv) :
      SplitOperation(lchild_clv_top,
                     rchild_clv_top,
                     lchild_clv_bot,
                     rchild_clv_bot,
                     lbrlen,
                     rbrlen,
                     prob_mat,
                     prob_mat,
                     rate_matrix,
                     rate_matrix,
                     parent_clv) {}

  SplitOperation(size_t parent_clv,
                 const std::shared_ptr<DispersionOperation>& l_ops,
                 const std::shared_ptr<DispersionOperation>& r_ops) :
      _lbranch_clv_index{l_ops->topCLVIndex()},
      _rbranch_clv_index{r_ops->topCLVIndex()},
      _parent_clv_index{parent_clv},
      _lbranch_op{l_ops},
      _rbranch_op{r_ops} {}

  void eval(const std::shared_ptr<Workspace>&);

  [[nodiscard]] auto getParentCLV() const -> size_t {
    return _parent_clv_index;
  }

  void printStatus(const std::shared_ptr<Workspace>& ws,
                   std::ostream& os,
                   size_t tabLevel = 0) const;

  [[nodiscard]] auto printStatus(const std::shared_ptr<Workspace>& ws,
                                 size_t tabLevel = 0) const -> std::string;

  [[nodiscard]] auto ready(const std::shared_ptr<Workspace>& ws) const -> bool {
    return _lbranch_op->ready(ws, _last_execution)
           && _rbranch_op->ready(ws, _last_execution);
  }

  auto getLock() -> std::mutex& { return _lock; }

  [[nodiscard]] auto getExpmOperations() const
      -> std::vector<std::shared_ptr<ExpmOperation>> {
    auto ret = _lbranch_op->getExpmOperations();
    auto tmp = _rbranch_op->getExpmOperations();

    ret.reserve(ret.size() + tmp.size());
    for (auto& op : tmp) { ret.push_back(op); }

    return ret;
  }

  void getRateMatrixOperations(
      std::unordered_map<size_t, std::shared_ptr<MakeRateMatrixOperation>>&
          ret) {
    _lbranch_op->emplaceRateMatrixOperations(ret);
    _rbranch_op->emplaceRateMatrixOperations(ret);
  }

  void fixDist(Range fix_dist) { _fixed_dist = fix_dist; }

  auto getFixedDist() -> Option<Range> const { return _fixed_dist; }

  void setExclAreas(Range e) { _excl_area_mask = e; }

  void setInclAreas(Range i) { _incl_area_mask = i; }

  void printGraph(std::ostream& os, size_t& index) const;

  auto evaluated(const std::shared_ptr<Workspace>& ws) const {
    auto eval = ws->lastUpdateCLV(_rbranch_clv_index) < _last_execution
                && ws->lastUpdateCLV(_lbranch_clv_index) < _last_execution
                && _last_execution < ws->lastUpdateCLV(_parent_clv_index);
    eval &= _lbranch_op ? _lbranch_op->evaluated(ws) : true;
    eval &= _rbranch_op ? _lbranch_op->evaluated(ws) : true;
    return eval;
  }

 private:
  size_t _lbranch_clv_index;
  size_t _rbranch_clv_index;
  size_t _parent_clv_index;

  std::shared_ptr<DispersionOperation> _lbranch_op;
  std::shared_ptr<DispersionOperation> _rbranch_op;

  Option<Range> _fixed_dist;
  Option<Range> _excl_area_mask;
  Option<Range> _incl_area_mask;

  std::mutex _lock;
  ClockTick _last_execution = 0;
};

class ReverseSplitOperation {
 public:
  ReverseSplitOperation(
      size_t bot_clv,  /* Where the result is stored */
      size_t ltop_clv, /* Where the result of the dispop is stored */
      size_t rtop_clv, /* Should be a _non_ reverse CLV */
      const std::shared_ptr<MakeRateMatrixOperation>& rate_matrix_op,
      size_t prob_matrix_index,
      size_t disp_clv_index,
      double brlen) :
      _bot_clv_index{bot_clv},
      _ltop_clv_index{ltop_clv},
      _rtop_clv_index{rtop_clv},
      _eval_clvs{true},
      _branch_op{{std::make_shared<DispersionOperation>(ltop_clv,
                                                        disp_clv_index,
                                                        brlen,
                                                        prob_matrix_index,
                                                        rate_matrix_op)}} {}

  ReverseSplitOperation(size_t bot_clv,
                        size_t rtop_clv,
                        const std::shared_ptr<DispersionOperation>& branch_op) :
      _bot_clv_index{bot_clv},
      _ltop_clv_index{branch_op->topCLVIndex()},
      _rtop_clv_index{rtop_clv},
      _eval_clvs{true},
      _branch_op{branch_op} {}

  explicit ReverseSplitOperation(
      const std::shared_ptr<DispersionOperation>& branch_op) :
      _bot_clv_index{std::numeric_limits<size_t>::max()},
      _ltop_clv_index{std::numeric_limits<size_t>::max()},
      _rtop_clv_index{std::numeric_limits<size_t>::max()},
      _eval_clvs{false},
      _branch_op{branch_op} {}

  explicit ReverseSplitOperation(size_t bot_clv,
                                 size_t ltop_clv,
                                 size_t rtop_clv) :
      _bot_clv_index{bot_clv},
      _ltop_clv_index{ltop_clv},
      _rtop_clv_index{rtop_clv},
      _eval_clvs{true},
      _branch_op{nullptr} {}

  void eval(const std::shared_ptr<Workspace>&);

  void printStatus(const std::shared_ptr<Workspace>& ws,
                   std::ostream& os,
                   size_t tabLevel = 0) const;

  [[nodiscard]] auto printStatus(const std::shared_ptr<Workspace>& ws,
                                 size_t tabLevel = 0) const -> std::string;

  void printGraph(std::ostream& os, size_t& index) const;

  void makeRootOperation(size_t clv_index) {
    _branch_op = nullptr;
    if (_ltop_clv_index == std::numeric_limits<size_t>::max()) {
      _ltop_clv_index = clv_index;
    }
    if (_rtop_clv_index == std::numeric_limits<size_t>::max()) {
      _rtop_clv_index = clv_index;
    }
  }

  [[nodiscard]] auto getStableCLV() const -> size_t { return _ltop_clv_index; }

  [[nodiscard]] auto ready(const std::shared_ptr<Workspace>& ws) const -> bool {
    bool ltop_clv_ready =
        (_branch_op != nullptr ? _branch_op->ready(ws, _last_execution) : true)
        || ws->lastUpdateCLV(_ltop_clv_index) > _last_execution;
    return ltop_clv_ready
           && (ws->lastUpdateCLV(_rtop_clv_index) > _last_execution);
  }

  auto getLock() -> std::mutex& { return _lock; }

  [[nodiscard]] auto getExpmOperations() const
      -> std::vector<std::shared_ptr<ExpmOperation>> {
    if (_branch_op) { return _branch_op->getExpmOperations(); }
    return {};
  }

  void fixDist(Range fix_dist) { _fixed_dist = fix_dist; }

  auto getFixedDist() -> Option<Range> const { return _fixed_dist; }

  void setInclAreas(Range i) { _incl_area_mask = i; }

  void setExclAreas(Range i) { _excl_area_mask = i; }

  auto evaluated(const std::shared_ptr<Workspace>& ws) const -> bool {
    return _last_execution > 0;
    bool eval = _last_execution > ws->lastUpdateCLV(_ltop_clv_index)
                && _last_execution > ws->lastUpdateCLV(_rtop_clv_index)
                && ws->lastUpdateCLV(_bot_clv_index) > _last_execution;
    eval &= _branch_op ? _branch_op->evaluated(ws) : true;
    return eval;
  }

 private:
  size_t _bot_clv_index;
  size_t _ltop_clv_index;
  size_t _rtop_clv_index;
  bool _eval_clvs;

  std::shared_ptr<DispersionOperation> _branch_op;

  Option<Range> _fixed_dist;
  Option<Range> _incl_area_mask;
  Option<Range> _excl_area_mask;

  std::mutex _lock;
  ClockTick _last_execution = 0;
};

class LLHGoal {
 public:
  LLHGoal(size_t root_clv, size_t prior_index) :
      _root_clv_index{root_clv},
      _prior_index{prior_index} {}

  void eval(const std::shared_ptr<Workspace>&);

  [[nodiscard]] auto result() const -> double { return _result; }

  size_t _root_clv_index;
  size_t _prior_index;
  double _result = -std::numeric_limits<double>::infinity();

  [[nodiscard]] auto ready(const std::shared_ptr<Workspace>&) const -> bool;
  void printGraph(std::ostream& os, size_t& index) const;

 private:
  ClockTick _last_execution = 0;
};

class StateLHGoal {
 public:
  StateLHGoal(size_t node_id,
              size_t parent_clv,
              size_t lchild_clv,
              size_t rchild_clv) :
      _parent_clv_index{parent_clv},
      _lchild_clv_index{lchild_clv},
      _rchild_clv_index{rchild_clv},
      _node_id{node_id},
      _result{nullptr} {}

  StateLHGoal(const StateLHGoal&) = delete;
  auto operator=(const StateLHGoal&) -> StateLHGoal& = delete;

  StateLHGoal(StateLHGoal&& other) noexcept :
      _parent_clv_index{other._parent_clv_index},
      _lchild_clv_index{other._lchild_clv_index},
      _rchild_clv_index{other._rchild_clv_index},
      _node_id{other._node_id},
      _fixed_dist{other._fixed_dist},
      _excl_area_mask{other._excl_area_mask},
      _incl_area_mask{other._incl_area_mask},
      _result{std::move(other._result)},
      _states{other._states} {}

  auto operator=(StateLHGoal&&) -> StateLHGoal& = delete;

  void eval(const std::shared_ptr<Workspace>&);

  [[nodiscard]] auto result() const -> std::unique_ptr<LagrangeMatrixBase[]> {
    std::unique_ptr<LagrangeMatrixBase[]> ret{new LagrangeMatrixBase[_states]};
    for (size_t i = 0; i < _states; i++) { ret[i] = _result[i]; }
    return ret;
  }

  [[nodiscard]] auto ready(const std::shared_ptr<Workspace>&) const -> bool;

  void fixDist(Range dist) { _fixed_dist = dist; }

  auto getFixedDist() -> Option<Range> { return _fixed_dist; }

  void setInclAreas(Range dist) { _incl_area_mask = dist; }

  void setExclAreas(Range dist) { _excl_area_mask = dist; }

  [[nodiscard]] auto node_id() const -> size_t { return _node_id; }

  void printGraph(std::ostream& os, size_t& index) const;

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;

  size_t _node_id;

  Option<Range> _fixed_dist;
  Range _excl_area_mask = 0;
  Range _incl_area_mask = 0;

  ClockTick _last_execution = 0;

  std::unique_ptr<LagrangeMatrixBase[]> _result;
  size_t _states{0};
};

class SplitLHGoal {
 public:
  SplitLHGoal(size_t node_id,
              size_t parent_clv,
              size_t lchild_clv,
              size_t rchild_clv) :
      _parent_clv_index{parent_clv},
      _lchild_clv_index{lchild_clv},
      _rchild_clv_index{rchild_clv},
      _node_id{node_id} {}

  void eval(const std::shared_ptr<Workspace>&);

  auto result() const -> SplitReturn { return _result; }

  auto ready(const std::shared_ptr<Workspace>&) const -> bool;

  void fixDist(Range dist) { _fixed_dist = dist; }

  void setInclAreas(Range dist) { _incl_area_mask = dist; }

  void setExclAreas(Range dist) { _excl_area_mask = dist; }

  auto node_id() const -> size_t { return _node_id; }

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;

  size_t _node_id;

  Option<Range> _fixed_dist;
  Range _excl_area_mask = 0;
  Range _incl_area_mask = 0;

  ClockTick _last_execution = 0;

  SplitReturn _result;
};

using PeriodRateMatrixMap =
    std::unordered_map<size_t, std::shared_ptr<MakeRateMatrixOperation>>;

using BranchProbMatrixMap = std::unordered_map<std::pair<size_t, double>,
                                               std::shared_ptr<ExpmOperation>>;

}  // namespace lagrange
#endif
