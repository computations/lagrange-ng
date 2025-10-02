#pragma once

#include "AncSplit.hpp"
#include "Workspace.hpp"

namespace lagrange {

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
  using result_type = std::unique_ptr<LagrangeMatrixBase[]>;

  StateLHGoal(size_t node_id,
              size_t parent_clv,
              size_t lchild_clv,
              size_t rchild_clv) :
      _parent_clv_index{parent_clv},
      _lchild_clv_index{lchild_clv},
      _rchild_clv_index{rchild_clv},
      _node_id{node_id},
      _result{nullptr} {}

  StateLHGoal(const StateLHGoal&);
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

  auto operator=(StateLHGoal&& other) -> StateLHGoal& {
    std::swap(*this, other);

    return *this;
  }

  void eval(const std::shared_ptr<const Workspace>&);

  [[nodiscard]] auto copyResult() const -> result_type {
    std::unique_ptr<LagrangeMatrixBase[]> ret{new LagrangeMatrixBase[_states]};
    for (size_t i = 0; i < _states; i++) { ret[i] = _result[i]; }
    return ret;
  }

  [[nodiscard]] auto result() const -> std::span<LagrangeMatrixBase> {
    return {_result.get(), _result.get() + _states};
  }

  [[nodiscard]] auto result() -> result_type {
    result_type tmp;
    tmp.swap(_result);
    return tmp;
  }

  [[nodiscard]] auto ready(const std::shared_ptr<const Workspace>&) const -> bool;

  void fixDist(Range dist) { _fixed_dist = dist; }

  auto getFixedDist() -> std::optional<Range> { return _fixed_dist; }

  void setInclAreas(Range dist) { _incl_area_mask = dist; }

  void setExclAreas(Range dist) { _excl_area_mask = dist; }

  [[nodiscard]] auto node_id() const -> size_t { return _node_id; }

  void printGraph(std::ostream& os, size_t& index) const;

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;

  size_t _node_id;

  std::optional<Range> _fixed_dist;
  Range _excl_area_mask = 0;
  Range _incl_area_mask = 0;

  ClockTick _last_execution = 0;

  std::unique_ptr<LagrangeMatrixBase[]> _result;
  size_t _states{0};
};

class SplitLHGoal {
 public:
  using result_type = SplitReturn;

  SplitLHGoal(size_t node_id,
              size_t parent_clv,
              size_t lchild_clv,
              size_t rchild_clv) :
      _parent_clv_index{parent_clv},
      _lchild_clv_index{lchild_clv},
      _rchild_clv_index{rchild_clv},
      _node_id{node_id} {}

  void eval(const std::shared_ptr<const Workspace>&);

  auto result() const -> result_type { return _result; }

  /* Make a copy, and then return it to trigger copy elision */
  auto result() -> result_type {
    result_type tmp;
    std::swap(tmp, _result);
    return tmp;
  }

  auto ready(const std::shared_ptr<const Workspace>&) const -> bool;

  void fixDist(Range dist) { _fixed_dist = dist; }

  void setInclAreas(Range dist) { _incl_area_mask = dist; }

  void setExclAreas(Range dist) { _excl_area_mask = dist; }

  auto node_id() const -> size_t { return _node_id; }

 private:
  size_t _parent_clv_index;
  size_t _lchild_clv_index;
  size_t _rchild_clv_index;

  size_t _node_id;

  std::optional<Range> _fixed_dist;
  Range _excl_area_mask = 0;
  Range _incl_area_mask = 0;

  ClockTick _last_execution = 0;

  result_type _result;
};

template <typename T>
concept StreamableGoal =
    std::same_as<T, StateLHGoal> || std::same_as<T, SplitLHGoal>;

template <StreamableGoal G>
class StreamingGoal {
 public:
  StreamingGoal(G g) : _goal{g} {}

  StreamingGoal(const StreamingGoal<G>&) = default;
  StreamingGoal(StreamingGoal<G>&&) = default;

  StreamingGoal() = default;

  StreamingGoal<G>& operator=(const StreamingGoal<G>&) = default;
  StreamingGoal<G>& operator=(StreamingGoal<G>&&) = default;

  auto eval(const std::shared_ptr<const Workspace>& ws) -> G::result_type {
    _goal.eval(ws);
    return _goal.result();
  }

  auto nodeID() const -> size_t { return _goal.node_id(); }

 private:
  G _goal;
};
}  // namespace lagrange
