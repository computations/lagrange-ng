#include <cmath>
#include <iostream>
#include <string>

#include "Common.h"
#include "Context.h"
#include "RateModel.h"
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
  if (_forward_operations.size() == 0) {
    registerForwardOperations();
  }

  auto root_clv = (_forward_operations.end() - 1)->get_parent_clv();
  size_t frequency_index = _workspace->suggest_freq_vector_index();
  _lh_goal.emplace_back(root_clv, frequency_index);
}

void Context::registerStateLHGoal() {
  if (_reverse_operations.size() == 0) {
    registerBackwardOperations();
  }

  auto node_ids = _tree->traversePreorderInternalNodesOnly();
  for (auto nid : node_ids) {
    _state_lh_goal.emplace_back(_workspace->get_top_clv_reverse(nid),
                                _workspace->get_lchild_clv(nid),
                                _workspace->get_rchild_clv(nid));
  }
}

void Context::updateRates(const period_t& params) {
  _rate_matrix_op->update_rates(_workspace, params);
}

void Context::init() {
  _workspace->reserve();
  updateRates({1.0, 1.0});
}

void Context::registerTipClvs(
    const std::unordered_map<std::string, lagrange_dist_t>& dist_data) {
  if (_forward_operations.size() == 0) {
    throw std::runtime_error{
        "The forward operations need to be generated first"};
  }
  _workspace->reserve();
  _tree->assignTipData(*_workspace, dist_data);
}

void Context::computeForwardOperations() {
  for (auto& op : _forward_operations) {
    op.eval(_workspace);
  }
}

void Context::computeBackwardOperations() {
  size_t prior_index = _reverse_operations.begin()->getStableCLV();
  _workspace->clv(prior_index) = 1.0 / _workspace->states();
  for (auto& op : _reverse_operations) {
    op.eval(_workspace);
    std::cout << op.printStatus(_workspace) << std::endl;
  }
}

double Context::computeLH() {
  computeForwardOperations();
  return _lh_goal.begin()->eval(_workspace);
}

double Context::computeLLH() { return std::log(computeLH()); }

period_derivative_t Context::computeDLLH(double initial_lh) {
  period_derivative_t derivative;

  auto cur = _workspace->get_period_params(0);
  auto tmp = cur;
  tmp.dispersion_rate += 1e-4;
  updateRates(tmp);
  derivative.d_dispersion = computeLLH() - initial_lh;

  tmp = cur;
  tmp.extinction_rate += 1e-4;
  updateRates(tmp);
  derivative.d_extinction = computeLLH() - initial_lh;

  return derivative;
}

void Context::optimize() {
  while (true) {
    period_t current_params = _workspace->get_period_params(0);
    double curlh = computeLLH();
    auto gradient = computeDLLH(curlh);

    if (gradient.norm() < 1e-10) {
      break;
    }
    current_params.applyDerivative(gradient);
    updateRates(current_params);
  }
}

double Context::computeLHGoal() { return computeLH(); }

std::vector<lagrange_col_vector_t> Context::computeStateGoal() {
  computeBackwardOperations();
  std::vector<lagrange_col_vector_t> states;
  states.reserve(_state_lh_goal.size());
  for (auto& op : _state_lh_goal) {
    states.push_back(op.eval(_workspace));
  }
  return states;
}
