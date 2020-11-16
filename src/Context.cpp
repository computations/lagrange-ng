#include <cmath>
#include <iostream>
#include <string>

#include "Common.h"
#include "Context.h"
#include "Workspace.h"

void Context::registerForwardOperations() {
  _forward_operations =
      _tree->generateForwardOperations(*_workspace, _rate_matrix_op);
}

void Context::registerBackwardOperations() {
  _forward_operations =
      _tree->generateForwardOperations(*_workspace, _rate_matrix_op);
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
  if (_reverse_operations.size() != 0) {
    registerBackwardOperations();
  }

  std::vector<size_t> node_ids;
  _tree->traversePreorderInternalNodesOnly();
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

double Context::computeLH() {
  for (auto& op : _forward_operations) {
    op.eval(_workspace);
  }
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
