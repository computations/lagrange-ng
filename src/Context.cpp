#include "Context.hpp"

#include <nlopt.hpp>

#include "AncSplit.hpp"
#include "Common.hpp"
#include "Operation.hpp"
#include "Utils.hpp"
#include "WorkerState.hpp"

void Context::registerForwardOperations() {
  _forward_operations = _tree->generateForwardOperations(*_workspace);
  extractRateMatrixOperations();
}

void Context::extractRateMatrixOperations() {
  std::unordered_map<size_t, std::shared_ptr<MakeRateMatrixOperation>> tmp;
  for (const auto& op : _forward_operations) {
    op->getRateMatrixOperations(tmp);
  }

  _rate_matrix_ops.resize(tmp.size());
  for (const auto& kv : tmp) { _rate_matrix_ops[kv.first] = kv.second; }
}

void Context::registerBackwardOperations() {
  _reverse_operations = _tree->generateBackwardOperations(*_workspace);
}

void Context::registerLHGoal() {
  if (_forward_operations.empty()) { registerForwardOperations(); }

  auto root_clv = _forward_operations.back()->getParentCLV();
  size_t frequency_index = _workspace->suggestFreqVectorIndex();
  _llh_goal.emplace_back(root_clv, frequency_index);
}

void Context::registerStateLHGoal() {
  if (_reverse_operations.empty()) { registerBackwardOperations(); }

  auto cb = [&](Node& n) {
    _state_lh_goal.emplace_back(_workspace->getTopCLVReverse(n.getId()),
                                _workspace->getLeftChildCLV(n.getId()),
                                _workspace->getRightChildCLV(n.getId()));
    if (n.getIncludedAreas().hasValue()) {
      _state_lh_goal.back().setInclAreas(n.getIncludedAreas().get());
    }
    if (n.getFixedDist().hasValue()) {
      _state_lh_goal.back().fixDist(n.getFixedDist().get());
    }
  };

  _tree->applyPreorderInternalOnly(cb);
}

void Context::registerSplitLHGoal() {
  if (_reverse_operations.empty()) { registerBackwardOperations(); }

  auto cb = [&](Node& n) {
    _split_lh_goal.emplace_back(_workspace->getTopCLVReverse(n.getId()),
                                _workspace->getLeftChildCLV(n.getId()),
                                _workspace->getRightChildCLV(n.getId()));
    if (n.getIncludedAreas().hasValue()) {
      _split_lh_goal.back().setInclAreas(n.getIncludedAreas().get());
    }
    if (n.getFixedDist().hasValue()) {
      _split_lh_goal.back().fixDist(n.getFixedDist().get());
    }
  };

  _tree->applyPreorderInternalOnly(cb);
}

void Context::updateRates(const std::vector<PeriodParams>& params) {
  for (size_t i = 0; i < _rate_matrix_ops.size(); ++i) {
    _rate_matrix_ops[i]->updateRates(_workspace, params[i]);
    _rate_matrix_ops[i]->eval(_workspace);
  }
}

void Context::init() {
  _workspace->reserve();
  std::vector<PeriodParams> initial_rates(_workspace->rateMatrixCount(),
                                             {0.01, 0.01});
  _workspace->setPeriodParamsCount(_workspace->rateMatrixCount());
  updateRates(initial_rates);

  if (!_reverse_operations.empty()) {
    size_t prior_index = _reverse_operations.front()->getStableCLV();
    _workspace->setReversePrior(prior_index);
  }

  auto fixed_dist = _forward_operations.back()->getFixedDist();

  if (fixed_dist.hasValue()) {
    _workspace->setBaseFrequenciesByDist(0, fixed_dist.get());
  }
}

void Context::registerTipClvs(
    const std::unordered_map<std::string, lagrange_dist_t>& dist_data) {
  if (_forward_operations.empty()) {
    throw std::runtime_error{
        "The forward operations need to be generated first"};
  }
  if (!_workspace->reserved()) { _workspace->reserve(); }
  _tree->assignTipData(*_workspace, dist_data);
}

void Context::computeBackwardOperations(WorkerState& ts, WorkerContext& tc) {
  ts.work(WorkerMode::ComputeReverse, tc, _workspace);
}

auto Context::computeLLH(WorkerState& ts) -> double {
  auto tc = makeThreadContext();
  return computeLLH(ts, tc);
}

auto Context::computeLLH(WorkerState& ts, WorkerContext& tc) -> double {
  ts.work(WorkerMode::ComputeForward, tc, _workspace);
  ts.work(WorkerMode::ComputeLH, tc, _workspace);
  return _llh_goal.begin()->result();
}

void Context::optimizeAndComputeValues(WorkerState& ts, bool states,
                                       bool splits, bool output,
                                       const LagrangeOperationMode& mode) {
  ts.assign_threads();
  WorkerContext tc = makeThreadContext();
  /* This blocks all but the main thread from proceeding until the halt mode
   * is set, which means that all further code is only executed by one thread
   */
  if (!ts.masterThread()) {
    ts.work(tc, _workspace);
    return;
  }

  double initial_lh = computeLLH(ts, tc);

  if (mode == LagrangeOperationMode::EVALUATE) {
    if (output) {
      std::cout << "LLH: " << initial_lh << std::endl;
      auto params = currentParams();
      for (const auto& p : params) {
        std::cout << "Dispersion: " << p.dispersion_rate
                  << " Extinction: " << p.extinction_rate << std::endl;
      }
    }
  }

  if (mode == LagrangeOperationMode::OPTIMIZE) {
    if (output) { std::cout << "Initial LLH: " << initial_lh << std::endl; }

    double final_lh = optimize(ts, tc);

    if (output) { std::cout << "Final LLH: " << final_lh << std::endl; }
  }

  if (states || splits) {
    std::cout << "Computing reverse operations" << std::endl;
    computeBackwardOperations(ts, tc);
  }

  if (states) {
    std::cout << "Computing state goals" << std::endl;
    computeStateGoal(ts, tc);
  }
  if (splits) {
    std::cout << "Computing split goals" << std::endl;
    computeSplitGoal(ts, tc);
  }
  ts.haltThreads();
}

auto Context::optimize(WorkerState& ts, WorkerContext& tc) -> double {
  struct OptContext {
    Context& context;
    WorkerContext& tc;
    WorkerState& ts;
    size_t iter = 0;
  } oc{*this, tc, ts};

  const size_t dims = _workspace->rateMatrixCount() * 2;
  nlopt::opt opt(nlopt::LN_NELDERMEAD, dims);
  auto objective = [](const std::vector<double>& x, std::vector<double>& grad,
                      void* f_data) -> double {
    (void)(grad);
    auto* obj = static_cast<OptContext*>(f_data);
    std::vector<PeriodParams> period_paramters(
        obj->context._workspace->rateMatrixCount());
    for (size_t i = 0; i < period_paramters.size(); ++i) {
      period_paramters[i] = {x[2 * i], x[2 * i + 1]};
    }

    obj->context.updateRates(period_paramters);
    double llh = obj->context.computeLLH(obj->ts, obj->tc);
    if (obj->iter % 10 == 0) {
      std::cout << "Current LLH: " << llh << std::endl;
    }
    if (std::isnan(llh)) {
      throw std::runtime_error{"Log likelihood is not not a number"};
    }
    obj->iter += 1;
    return llh;
  };

  opt.set_max_objective(objective, &oc);

  std::vector<double> lower_bounds(dims, 1e-7);
  opt.set_lower_bounds(lower_bounds);

  opt.set_ftol_rel(_lh_epsilon);

  std::vector<double> results(dims, 0.01);
  double obj_val = 0;

  opt.optimize(results, obj_val);

  return obj_val;
}

/* Computes the registered state goals
 * If there are no reverse state goals, this operation does nothing. It also
 * requires that the forward operations be computed first.
 *
 * Returns:
 *   a map (really a vector) from node id to result. The result is a col
 * vector indexed by lagrange_dist_t, where each entry contains the likelihood
 * of that particular distribution at that node.
 */
auto Context::getStateResults() const
    -> std::vector<std::unique_ptr<LagrangeMatrixBase[]>> {
  std::vector<std::unique_ptr<LagrangeMatrixBase[]>> states;
  states.reserve(_state_lh_goal.size());

  for (auto& op : _state_lh_goal) { states.push_back(op.result()); }

  return states;
}

auto Context::getSplitResults() const -> lagrange_split_list_t {
  lagrange_split_list_t splits;
  splits.reserve(_split_lh_goal.size());

  for (auto& op : _split_lh_goal) { splits.push_back(op.result()); }

  return splits;
}

void Context::computeStateGoal(WorkerState& ts, WorkerContext& tc) {
  ts.work(WorkerMode::ComputeStateGoal, tc, _workspace);
}
void Context::computeSplitGoal(WorkerState& ts, WorkerContext& tc) {
  ts.work(WorkerMode::ComputeSplitGoal, tc, _workspace);
}

auto Context::currentParams() const -> std::vector<PeriodParams> {
  return _workspace->getPeriodParams();
}

#if 0
std::string Context::treeCLVStatus() const {
  std::stringstream oss;
  auto node_ids = _tree->traversePreorderInternalNodesOnly();
  for (auto nid : node_ids) {
    oss << "node id: " << nid << "\n"
        << _workspace->report_node_vecs(nid) << std::endl;
  }
  return oss.str();
}
#endif

auto Context::computeStateGoal(WorkerState& ts)
    -> std::vector<std::unique_ptr<LagrangeMatrixBase[]>> {
  auto tc = makeThreadContext();
  computeBackwardOperations(ts, tc);
  computeStateGoal(ts, tc);
  return getStateResults();
}

auto Context::computeSplitGoal(WorkerState& ts) -> lagrange_split_list_t {
  auto tc = makeThreadContext();
  computeBackwardOperations(ts, tc);
  computeSplitGoal(ts, tc);
  return getSplitResults();
}

void Context::useArnoldi(bool mode_set, bool adaptive) const {
  for (auto& op : _forward_operations) {
    auto expm_ops = op->getExpmOperations();
    for (auto& eop : expm_ops) {
      eop->setArnoldiMode(mode_set);
      eop->setAdaptive(adaptive);
    }
  }
  for (auto& op : _reverse_operations) {
    auto expm_ops = op->getExpmOperations();
    for (auto& eop : expm_ops) {
      eop->setArnoldiMode(mode_set);
      eop->setAdaptive(adaptive);
    }
  }
}

size_t Context::getPeriodCount() const { return _rate_matrix_ops.size(); }
