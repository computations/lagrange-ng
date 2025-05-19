#include "Context.hpp"

#include <nlopt.hpp>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "AncSplit.hpp"
#include "Common.hpp"
#include "ConfigFile.hpp"
#include "MRCA.hpp"
#include "Operation.hpp"
#include "Utils.hpp"
#include "WorkerState.hpp"
#include "Workspace.hpp"
#include "logger.hpp"

namespace lagrange {

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

auto Context::makeStateGoalCB() -> std::function<void(Node&)> {
  auto cb = [&](Node& n) {
    _state_lh_goal.emplace_back(n.getId(),
                                _workspace->getTopCLVReverse(n.getId()),
                                _workspace->getLeftChildCLV(n.getId()),
                                _workspace->getRightChildCLV(n.getId()));
    if (n.getIncludedAreas().hasValue()) {
      _state_lh_goal.back().setInclAreas(n.getIncludedAreas().get());
    }
    if (n.getExcludedAreas().hasValue()) {
      _state_lh_goal.back().setExclAreas(n.getExcludedAreas().get());
    }
    if (n.getFixedDist().hasValue()) {
      _state_lh_goal.back().fixDist(n.getFixedDist().get());
    }
  };
  return cb;
}

auto Context::makeSplitGoalCB() -> std::function<void(Node&)> {
  auto cb = [&](Node& n) {
    _split_lh_goal.emplace_back(n.getId(),
                                _workspace->getTopCLVReverse(n.getId()),
                                _workspace->getLeftChildCLV(n.getId()),
                                _workspace->getRightChildCLV(n.getId()));
    if (n.getIncludedAreas().hasValue()) {
      _split_lh_goal.back().setInclAreas(n.getIncludedAreas().get());
    }
    if (n.getExcludedAreas().hasValue()) {
      _split_lh_goal.back().setExclAreas(n.getExcludedAreas().get());
    }
    if (n.getFixedDist().hasValue()) {
      _split_lh_goal.back().fixDist(n.getFixedDist().get());
    }
  };
  return cb;
}

void Context::registerStateLHGoal() {
  if (!_state_lh_goal.empty()) { return; }
  if (_reverse_operations.empty()) { registerBackwardOperations(); }

  auto cb = makeStateGoalCB();

  _tree->applyPreorderInternalOnly(cb);
}

void Context::registerStateLHGoal(
    const std::unordered_set<std::string>& mrca_keys,
    const std::unordered_map<std::string, std::shared_ptr<MRCAEntry>>&
        mrca_map) {
  if (!_state_lh_goal.empty()) { return; }
  auto cb = makeStateGoalCB();
  registerGoals(mrca_keys, mrca_map, cb);
}

void Context::registerSplitLHGoal() {
  auto cb = makeSplitGoalCB();
  registerGoals(cb);
}

void Context::registerSplitLHGoal(
    const std::unordered_set<std::string>& mrca_keys,
    const std::unordered_map<std::string, std::shared_ptr<MRCAEntry>>&
        mrca_map) {
  auto cb = makeSplitGoalCB();

  registerGoals(mrca_keys, mrca_map, cb);
}

void Context::registerGoals(const std::function<void(Node&)>& cb) {
  if (_reverse_operations.empty()) { registerBackwardOperations(); }
  _tree->applyPreorderInternalOnly(cb);
}

void Context::registerGoals(
    const std::unordered_set<std::string>& mrca_keys,
    const std::unordered_map<std::string, std::shared_ptr<MRCAEntry>>& mrca_map,
    const std::function<void(Node&)>& cb) {
  if (_reverse_operations.empty()) { registerBackwardOperations(); }

  for (const auto& mrca_key : mrca_keys) {
    auto n = _tree->getMRCA(mrca_map.at(mrca_key));
    n->applyCB(cb);
  }
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

SetCLVStatus Context::registerTipClvs(
    const std::unordered_map<std::string, Range>& dist_data) {
  if (_forward_operations.empty()) {
    throw std::runtime_error{
        "The forward operations need to be generated first"};
  }
  if (!_workspace->reserved()) { _workspace->reserve(); }
  auto res = _tree->assignTipData(*_workspace, dist_data);
  if (res == SetCLVStatus::failed) { LOG(ERROR, "Failed to assign data"); }
  return res;
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

void Context::haltThreads(WorkerState& ts, WorkerContext& tc) {
  ts.work(WorkerMode::Halt, tc, _workspace);
}

void Context::optimizeAndComputeValues(WorkerState& ts,
                                       WorkerContext& tc,
                                       bool states,
                                       bool splits,
                                       const LagrangeOperationMode& mode) {
  ts.assign_threads();
  /* This blocks all but the main thread from proceeding until the halt mode
   * is set, which means that all further code is only executed by one thread
   */
  if (!ts.masterThread()) {
    ts.work(tc, _workspace);
    return;
  }

  double initial_lh = computeLLH(ts, tc);

  if (mode == LagrangeOperationMode::EVALUATE) {
    LOG(INFO, "LLH: {:.7}", initial_lh);
    auto params = currentParams();
    for (const auto& p : params) {
      LOG(INFO,
          "Dispersion: {:.7}, Extinction: {:.7}",
          p.dispersion_rate,
          p.extinction_rate);
    }
  }

  if (mode == LagrangeOperationMode::OPTIMIZE) {
    LOG(INFO, "Initial LLH: {:.7}", initial_lh);

    double final_lh = optimize(ts, tc);

    LOG(INFO, "Final LLH: {:.7}", final_lh);
  }

  if (states || splits) {
    LOG(INFO, "Computing reverse operations");
    computeBackwardOperations(ts, tc);
  }

  if (states) {
    LOG(INFO, "Computing ancestral states");
    computeStateGoal(ts, tc);
  }

  if (splits) {
    LOG(INFO, "Computing ancestral splits");
    computeSplitGoal(ts, tc);
  }

  LOG_INFO("Halting workers");
  haltThreads(ts, tc);
}

auto Context::optimize(WorkerState& ts, WorkerContext& tc) -> double {
  struct OptContext {
    Context& context;
    WorkerContext& tc;
    WorkerState& ts;
    size_t iter = 0;
    size_t f_evals = 1;
  } oc{*this, tc, ts};

  const size_t dims = _workspace->rateMatrixCount() * 2;
  nlopt::opt opt(_opt_method, dims);
  auto objective = [](const std::vector<double>& x,
                      std::vector<double>& grad,
                      void* f_data) -> double {
    auto* obj = static_cast<OptContext*>(f_data);
    std::vector<PeriodParams> period_paramters(
        obj->context._workspace->rateMatrixCount());
    for (size_t i = 0; i < period_paramters.size(); ++i) {
      period_paramters[i] = {x[2 * i], x[2 * i + 1]};
    }

    obj->context.updateRates(period_paramters);
    double llh = obj->context.computeLLH(obj->ts, obj->tc);
    obj->f_evals += 1;

    if (!grad.empty()) {
      constexpr double step = 1e-7;

      for (size_t i = 0; i < grad.size(); ++i) {
        auto tmp_params = period_paramters;
        tmp_params[i / 2][i % 2] += step;
        obj->context.updateRates(tmp_params);
        grad[i] = (obj->context.computeLLH(obj->ts, obj->tc) - llh) / step;
        obj->f_evals += 1;
      }
    }

    if (obj->iter % 10 == 0) { LOG(PROGRESS, "Current LLH: {:.7}", llh); }
    if (std::isnan(llh)) {
      LOG(ERROR, "Log liklihood is not a number");
      throw std::runtime_error{"Log likelihood is not a number"};
    }
    obj->iter += 1;
    return llh;
  };

  opt.set_max_objective(objective, &oc);

  std::vector<double> lower_bounds(dims, 1e-7);
  opt.set_lower_bounds(lower_bounds);

  std::vector<double> upper_bounds(dims, 1e2);
  opt.set_upper_bounds(upper_bounds);

  opt.set_ftol_rel(_lh_epsilon);

  std::vector<double> results(dims, 0.01);
  double obj_val = 0;

  try {
    opt.optimize(results, obj_val);
  } catch (nlopt::roundoff_limited& err) {
    LOG(WARNING,
        "NLopt finished with limited roundoff, results might be incorrect");
  }

  LOG(INFO, "Finished optimization with {} likelihood evaluations", oc.f_evals);

  return obj_val;
}

/* Computes the registered state goals
 * If there are no reverse state goals, this operation does nothing. It also
 * requires that the forward operations be computed first.
 *
 * Returns:
 *   a map from node id to result. The result is a col
 * vector indexed by lagrange_dist_t, where each entry contains the likelihood
 * of that particular distribution at that node.
 */
auto Context::getStateResults() const
    -> std::unordered_map<size_t, std::unique_ptr<LagrangeMatrixBase[]>> {
  std::unordered_map<size_t, std::unique_ptr<LagrangeMatrixBase[]>> states;
  states.reserve(_state_lh_goal.size());

  for (const auto& op : _state_lh_goal) { states[op.node_id()] = op.result(); }

  return states;
}

auto Context::getSplitResults() const -> SplitReturnList {
  SplitReturnList splits;
  splits.reserve(_split_lh_goal.size());

  for (const auto& op : _split_lh_goal) { splits[op.node_id()] = op.result(); }

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

auto Context::computeAndGetStateGoals(WorkerState& ts, WorkerContext& tc)
    -> std::unordered_map<size_t, std::unique_ptr<LagrangeMatrixBase[]>> {
  computeBackwardOperations(ts, tc);
  computeStateGoal(ts, tc);
  return getStateResults();
}

auto Context::computeSplitGoal(WorkerState& ts) -> SplitReturnList {
  auto tc = makeThreadContext();
  computeBackwardOperations(ts, tc);
  computeSplitGoal(ts, tc);
  return getSplitResults();
}

void Context::useArnoldi(bool mode_set, bool adaptive) const {
  for (const auto& op : _forward_operations) {
    auto expm_ops = op->getExpmOperations();
    for (auto& eop : expm_ops) {
      eop->setArnoldiMode(mode_set);
      eop->setAdaptive(adaptive);
    }
  }
  for (const auto& op : _reverse_operations) {
    auto expm_ops = op->getExpmOperations();
    for (auto& eop : expm_ops) {
      eop->setArnoldiMode(mode_set);
      eop->setAdaptive(adaptive);
    }
  }
}

auto Context::getPeriodCount() const -> size_t {
  return _rate_matrix_ops.size();
}

void Context::set_opt_method(const OptimizationMethod& m) {
  switch (m) {
    case OptimizationMethod::NELDER_MEAD:
      _opt_method = nlopt::LN_NELDERMEAD;
      break;
    case OptimizationMethod::COBYLA:
      _opt_method = nlopt::LN_COBYLA;
      break;
    case OptimizationMethod::BOBYQA:
      _opt_method = nlopt::LN_BOBYQA;
      break;
    case OptimizationMethod::BFGS:
      _opt_method = nlopt::LD_LBFGS;
      break;
    case OptimizationMethod::DIRECT:
      _opt_method = nlopt::GN_DIRECT;
      break;
    case OptimizationMethod::STOGO:
      _opt_method = nlopt::GD_STOGO;
      break;
    case OptimizationMethod::UNKNOWN:
    default:
      LOG(ERROR, "Unknown optimization method");
      throw std::runtime_error{"Unknown optimization method"};
  }
  LOG_INFO("Using {} for optimization", nlopt::algorithm_name(_opt_method));
}

auto Context::computeDerivative() const -> bool {
  return _opt_method == nlopt::LD_LBFGS;
}

void Context::dumpReverseGraph(std::ostream& os) const {
  size_t index = 0;
  os << "digraph g {\n";
  os << R"(node [shape = box];)" << "\n";
  os << R"(graph [rankdir = "BT"];)" << "\n";
  for (const auto& slhg : _state_lh_goal) { slhg.printGraph(os, index); }
  for (const auto& rop : _reverse_operations) { rop->printGraph(os, index); }
  os << "}";
}

void Context::dumpForwardGraph(std::ostream& os) const {
  size_t index = 0;
  os << "digraph g {\n";
  os << R"(node [shape = box];)" << "\n";
  os << R"(graph [rankdir = "TB"];)" << "\n";
  for (const auto& lhg : _llh_goal) { lhg.printGraph(os, index); }
  for (const auto& fop : _forward_operations) { fop->printGraph(os, index); }
  os << "}";
}

}  // namespace lagrange
