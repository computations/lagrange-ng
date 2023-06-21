#include "Context.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <nlopt.hpp>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "AncSplit.h"
#include "Common.h"
#include "Operation.h"
#include "Utils.h"
#include "WorkerState.h"
#include "Workspace.h"

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

  auto root_clv = _forward_operations.back()->get_parent_clv();
  size_t frequency_index = _workspace->suggest_freq_vector_index();
  _llh_goal.emplace_back(root_clv, frequency_index);
}

void Context::registerStateLHGoal() {
  if (_reverse_operations.empty()) { registerBackwardOperations(); }

  auto cb = [&](Node& n) {
    _state_lh_goal.emplace_back(_workspace->get_top_clv_reverse(n.getId()),
                                _workspace->get_lchild_clv(n.getId()),
                                _workspace->get_rchild_clv(n.getId()));
    if (n.getInclAreas().has_value()) {
      _state_lh_goal.back().setInclAreas(n.getInclAreas().get());
    }
    if (n.getFixedDist().has_value()) {
      _state_lh_goal.back().fixDist(n.getFixedDist().get());
    }
  };

  _tree->applyPreorderInternalOnly(cb);
}

void Context::registerSplitLHGoal() {
  if (_reverse_operations.empty()) { registerBackwardOperations(); }

  auto cb = [&](Node& n) {
    _split_lh_goal.emplace_back(_workspace->get_top_clv_reverse(n.getId()),
                                _workspace->get_lchild_clv(n.getId()),
                                _workspace->get_rchild_clv(n.getId()));
    if (n.getInclAreas().has_value()) {
      _split_lh_goal.back().setInclAreas(n.getInclAreas().get());
    }
    if (n.getFixedDist().has_value()) {
      _split_lh_goal.back().fixDist(n.getFixedDist().get());
    }
  };

  _tree->applyPreorderInternalOnly(cb);
}

void Context::updateRates(const std::vector<period_params_t>& params) {
  for (size_t i = 0; i < _rate_matrix_ops.size(); ++i) {
    _rate_matrix_ops[i]->update_rates(_workspace, params[i]);
    _rate_matrix_ops[i]->eval(_workspace);
  }
}

void Context::init() {
  _workspace->reserve();
  std::vector<period_params_t> initial_rates(_workspace->rate_matrix_count(),
                                             {0.01, 0.01});
  _workspace->set_period_params_count(_workspace->rate_matrix_count());
  updateRates(initial_rates);

  if (!_reverse_operations.empty()) {
    size_t prior_index = _reverse_operations.front()->getStableCLV();
    _workspace->set_reverse_prior(prior_index);
  }

  auto fixed_dist = _forward_operations.back()->getFixedDist();

  if (fixed_dist.has_value()) {
    _workspace->set_base_frequencies_by_dist(0, fixed_dist.get());
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
                                       const lagrange_operation_mode& mode) {
  ts.assign_threads();
  WorkerContext tc = makeThreadContext();
  /* This blocks all but the main thread from proceeding until the halt mode
   * is set, which means that all further code is only executed by one thread
   */
  if (!ts.master_thread()) {
    ts.work(tc, _workspace);
    return;
  }

  double initial_lh = computeLLH(ts, tc);

  if (mode == lagrange_operation_mode::EVALUATE) {
    if (output) {
      std::cout << "LLH: " << initial_lh << std::endl;
      auto params = currentParams();
      for (const auto& p : params) {
        std::cout << "Dispersion: " << p.dispersion_rate
                  << " Extinction: " << p.extinction_rate << std::endl;
      }
    }
  }

  if (mode == lagrange_operation_mode::OPTIMIZE) {
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
  ts.halt_threads();
}

auto Context::optimize(WorkerState& ts, WorkerContext& tc) -> double {
  struct OptContext {
    Context& context;
    WorkerContext& tc;
    WorkerState& ts;
    size_t iter = 0;
  } oc{*this, tc, ts};

  const size_t dims = _workspace->rate_matrix_count() * 2;
  nlopt::opt opt(nlopt::LN_NELDERMEAD, dims);
  auto objective = [](const std::vector<double>& x, std::vector<double>& grad,
                      void* f_data) -> double {
    (void)(grad);
    auto* obj = static_cast<OptContext*>(f_data);
    std::vector<period_params_t> period_paramters(
        obj->context._workspace->rate_matrix_count());
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
auto Context::getStateResults()
    -> std::vector<std::unique_ptr<lagrange_matrix_base_t[]>> {
  std::vector<std::unique_ptr<lagrange_matrix_base_t[]>> states;
  states.reserve(_state_lh_goal.size());

  for (auto& op : _state_lh_goal) { states.push_back(op.result()); }

  return states;
}

auto Context::getSplitResults() -> lagrange_split_list_t {
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

auto Context::currentParams() const -> std::vector<period_params_t> {
  return _workspace->get_period_params();
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
    -> std::vector<std::unique_ptr<lagrange_matrix_base_t[]>> {
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
