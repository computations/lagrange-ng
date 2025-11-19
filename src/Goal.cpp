#include "Goal.hpp"

#include "Operation.hpp"

namespace lagrange {

void LLHGoal::eval(const std::shared_ptr<Workspace> &ws) {
  double rho = cblas_ddot(ws->CLVSize(),
                          ws->CLV(_root_clv_index),
                          1,
                          ws->getBaseFrequencies(_prior_index),
                          1);
  lagrange_assert(rho > 0.0,
                  "Likelihood computations failed with an invalid likelihood");
  _result = std::log(rho)
            - lagrange_scaling_factor_log * ws->CLVScalar(_root_clv_index);

  _last_execution = ws->advanceClock();
}

auto LLHGoal::ready(const std::shared_ptr<Workspace> &ws) const -> bool {
  return ws->lastUpdateCLV(_root_clv_index) > _last_execution;
}

void LLHGoal::printGraph(std::ostream &os, size_t &index) const {
  os << std::format(R"({} [label = "LLHGoal {}"];)", index, index) << "\n";
  os << std::format("clv{} -> {};\n", _root_clv_index, index);
  os << std::format("bf{} -> {};\n", _prior_index, index);

  index++;
}

void StateLHGoal::eval(const std::shared_ptr<const Workspace> &ws) {
  if (_last_execution == 0) {
    _result.reset(new LagrangeMatrixBase[ws->restrictedStateCount()]);
    _states = ws->restrictedStateCount();
    for (Range i = 0; i < _states; ++i) { _result[i] = 0.0; }
  }

  size_t tmp_scalar = 0;

  weighted_combine(ws->CLV(_lchild_clv_index),
                   ws->CLV(_rchild_clv_index),
                   _states,
                   ws->regions(),
                   ws->maxAreas(),
                   _result.get(),
                   ws->CLVScalar(_lchild_clv_index),
                   ws->CLVScalar(_rchild_clv_index),
                   tmp_scalar,
                   _fixed_dist,
                   _excl_area_mask,
                   _incl_area_mask);

  tmp_scalar += ws->CLVScalar(_parent_clv_index);
  auto *result = _result.get();

  size_t dist_index = 0;
  size_t dist_val = 0;
  for (dist_index = 0, dist_val = 0; dist_index < ws->restrictedStateCount();
       dist_val = next_dist(dist_val, ws->maxAreas()), ++dist_index) {
    double tmp_val = result[dist_index];
    double parent_val = ws->CLV(_parent_clv_index)[dist_index];
    result[dist_index] = std::log(tmp_val * parent_val)
                         - tmp_scalar * lagrange_scaling_factor_log;
    assert(dist_index == 0 || (dist_val & _incl_area_mask)
           || !(dist_val & _excl_area_mask)
           || std::isfinite(result[dist_index]));
  }
  _last_execution = ws->readClock();
}

auto StateLHGoal::ready(const std::shared_ptr<const Workspace> &ws) const
    -> bool {
  return ws->lastUpdateCLV(_lchild_clv_index) > _last_execution
         && ws->lastUpdateCLV(_rchild_clv_index) > _last_execution
         && ws->lastUpdateCLV(_parent_clv_index) > _last_execution;
}

void StateLHGoal::printGraph(std::ostream &os, size_t &index) const {
  os << std::format(R"({} [label = "State Goal {}"];)", index, index) << "\n";

  os << std::format("clv{} -> {};\n", _parent_clv_index, index);
  os << std::format("clv{} -> {};\n", _lchild_clv_index, index);
  os << std::format("clv{} -> {};\n", _rchild_clv_index, index);

  index++;
}

void SplitLHGoal::eval(const std::shared_ptr<const Workspace> &ws) {
  SplitReturn ret;
  ret.reserve(ws->restrictedStateCount());

  const auto &parent_clv = ws->CLV(_parent_clv_index);
  const auto &lchild_clv = ws->CLV(_lchild_clv_index);
  const auto &rchild_clv = ws->CLV(_rchild_clv_index);
  const size_t max_areas = ws->maxAreas();

  std::vector<RegionSplit> splits;

  auto inv_dist_map = invert_dist_map(ws->regions(), max_areas);

  size_t index = 0;
  Range dist = 0;
  const Range max_dist = 1ULL << ws->regions();
  for (index = 0, dist = 0; dist < max_dist;
       index++, dist = next_dist(dist, max_areas)) {
    if (_fixed_dist && _fixed_dist.value() != index) { continue; }

    if (!check_excl_dist(index, _excl_area_mask)) { continue; }
    if (!check_incl_dist(index, _incl_area_mask)) { continue; }

    std::vector<AncSplit> anc_split_vec;
    generate_splits(dist, ws->regions(), splits);
    double weight = 1.0 / splits.size();
    for (auto sp : splits) {
      AncSplit anc_split(dist, sp.left, sp.right, weight);

      auto left_index = inv_dist_map[sp.left];
      auto right_index = inv_dist_map[sp.right];

      double lh = parent_clv[index] * lchild_clv[left_index]
                  * rchild_clv[right_index] * weight;
      double loglh = std::log(lh);
      assert(!std::isnan(loglh));
      anc_split.setLikelihood(loglh);
      anc_split_vec.push_back(anc_split);
    }
    ret[dist] = anc_split_vec;
  }
  _result = ret;
}

auto SplitLHGoal::ready(const std::shared_ptr<const Workspace> &ws) const
    -> bool {
  return ws->lastUpdateCLV(_lchild_clv_index) > _last_execution
         && ws->lastUpdateCLV(_rchild_clv_index) > _last_execution
         && ws->lastUpdateCLV(_parent_clv_index) > _last_execution;
}
}  // namespace lagrange
