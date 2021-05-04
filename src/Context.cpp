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
#include "Context.h"
#include "Operation.h"
#include "ThreadState.h"
#include "Workspace.h"

void Context::registerForwardOperations() {
  _forward_operations =
      _tree->generateForwardOperations(*_workspace, _rate_matrix_op);
}

void Context::registerBackwardOperations() {
  _reverse_operations =
      _tree->generateBackwardOperations(*_workspace, _rate_matrix_op);
}

void Context::registerLHGoal() {
  if (_forward_operations.size() == 0) { registerForwardOperations(); }

  auto root_clv = (*(_forward_operations.end() - 1))->get_parent_clv();
  size_t frequency_index = _workspace->suggest_freq_vector_index();
  _llh_goal.emplace_back(root_clv, frequency_index);
}

void Context::registerStateLHGoal() {
  if (_reverse_operations.size() == 0) { registerBackwardOperations(); }

  auto node_ids = _tree->traversePreorderInternalNodesOnly();
  for (auto nid : node_ids) {
    _state_lh_goal.emplace_back(_workspace->get_top_clv_reverse(nid),
                                _workspace->get_lchild_clv(nid),
                                _workspace->get_rchild_clv(nid));
  }
}

void Context::updateRates(const period_t& params) {
  _rate_matrix_op->update_rates(_workspace, params);
  _rate_matrix_op->eval(_workspace);
}

void Context::init() {
  _workspace->reserve();
  updateRates({0.01, 0.01});

  if (_reverse_operations.size() != 0) {
    size_t prior_index = _reverse_operations.front()->getStableCLV();
    _workspace->set_reverse_prior(prior_index);
  }
}

void Context::registerTipClvs(
    const std::unordered_map<std::string, lagrange_dist_t>& dist_data) {
  if (_forward_operations.size() == 0) {
    throw std::runtime_error{
        "The forward operations need to be generated first"};
  }
  if (!_workspace->reserved()) { _workspace->reserve(); }
  _tree->assignTipData(*_workspace, dist_data);
}

void Context::computeForwardOperations(ThreadState& ts, ThreadContext& tc) {
  ts.work(ThreadMode::ComputeForward, tc, _workspace);
}

void Context::computeBackwardOperations(ThreadState& ts, ThreadContext& tc) {
  ts.work(ThreadMode::ComputeReverse, tc, _workspace);
}

double Context::computeLLH(ThreadState& ts) {
  auto tc = makeThreadContext();
  return computeLLH(ts, tc);
}

double Context::computeLLH(ThreadState& ts, ThreadContext& tc) {
  ts.work(ThreadMode::ComputeForward, tc, _workspace);
  ts.work(ThreadMode::ComputeLH, tc, _workspace);
  return _llh_goal.begin()->result();
}

void Context::optimizeAndComputeValues(ThreadState& ts, bool states,
                                       bool splits, bool output) {
  ThreadContext tc = makeThreadContext();
  /* This blocks all but the main thread from proceeding until the halt mode
   * is set, which means that
   */
  if (!ts.master_thread()) {
    ts.work(tc, _workspace);
    return;
  }

  double initial_lh = computeLLH(ts, tc);

  if (output) { std::cout << "Initial LH: " << initial_lh << std::endl; }

  double final_lh = optimize(ts, tc);

  if (output) { std::cout << "Final LH: " << final_lh << std::endl; }

  if (states || splits) { computeBackwardOperations(ts, tc); }

  if (states) { computeStateGoal(ts, tc); }
  if (splits) { computeSplitGoal(ts, tc); }
  ts.halt_threads();
}

double Context::optimize(ThreadState& ts, ThreadContext& tc) {
  struct OptContext {
    Context& context;
    ThreadContext& tc;
    ThreadState& ts;
  } oc{*this, tc, ts};

  nlopt::opt opt(nlopt::LN_SBPLX, 2);
  auto objective = [](const std::vector<double>& x, std::vector<double>& grad,
                      void* f_data) -> double {
    (void)(grad);
    auto obj = static_cast<OptContext*>(f_data);
    period_t p{x[0], x[1]};
    obj->context.updateRates(p);
    double llh = obj->context.computeLLH(obj->ts, obj->tc);
    // std::cout << p.toString() << ": " << llh << std::endl;
    if (std::isnan(llh)) {
      throw std::runtime_error{"Log likelihood is not not a number"};
    }
    return llh;
  };

  opt.set_max_objective(objective, &oc);
  opt.set_lower_bounds({1e-7, 1e-7});

  std::vector<double> results(2, 0.01);
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
std::vector<std::unique_ptr<lagrange_matrix_base_t>>
Context::getStateResults() {
  std::vector<std::unique_ptr<lagrange_matrix_base_t>> states;
  states.reserve(_state_lh_goal.size());

  for (auto& op : _state_lh_goal) { states.push_back(op.result()); }

  return states;
}

lagrange_split_list_t Context::getSplitResults() {
  lagrange_split_list_t splits;
  splits.reserve(_split_lh_goal.size());

  for (auto& op : _split_lh_goal) { splits.push_back(op.result()); }

  return splits;
}

void Context::computeStateGoal(ThreadState& ts, ThreadContext& tc) {
  ts.work(ThreadMode::ComputeStateGoal, tc, _workspace);
}
void Context::computeSplitGoal(ThreadState& ts, ThreadContext& tc) {
  ts.work(ThreadMode::ComputeSplitGoal, tc, _workspace);
}

period_t Context::currentParams() const {
  return _workspace->get_period_params(0);
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

std::vector<std::unique_ptr<lagrange_matrix_base_t>> Context::computeStateGoal(
    ThreadState& ts) {
  auto tc = makeThreadContext();
  computeStateGoal(ts, tc);
  return getStateResults();
}

lagrange_split_list_t Context::computeSplitGoal(ThreadState& ts) {
  auto tc = makeThreadContext();
  computeSplitGoal(ts, tc);
  return getSplitResults();
}
