/*
 * node.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include "Node.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Common.hpp"
#include "Fossil.hpp"
#include "MRCA.hpp"
#include "Operation.hpp"
#include "Utils.hpp"

namespace lagrange {

Node::Node() : _branch_length(0.0), _height(0.0), _id(0) {}

Node::Node(double bl, std::string inname) :
    _branch_length(bl),
    _height(0.0),
    _id(0),
    _label(std::move(inname)) {}

auto Node::isTip() const -> bool { return _children.empty(); }

auto Node::isInternal() const -> bool { return !_children.empty(); }

auto Node::getNumber() const -> size_t { return _id; }

auto Node::getBL() const -> double { return _branch_length; }

void Node::setBL(double bl) { _branch_length = bl; }

auto Node::hasChild(const std::shared_ptr<Node> &test) -> bool {
  return std::any_of(_children.begin(), _children.end(), [&test](auto n) {
    return n == test;
  });
}

auto Node::addChild(const std::shared_ptr<Node> &c) -> bool {
  if (!hasChild(c)) {
    _children.push_back(c);
    return true;
  }
  return false;
}

auto Node::removeChild(const std::shared_ptr<Node> &c) -> bool {
  if (hasChild(c)) {
    for (auto it = _children.begin(); it != _children.end(); it++) {
      if (*it == c) {
        _children.erase(it);
        return true;
      }
    }
  }
  return false;
}

auto Node::getChild(size_t c) const -> std::shared_ptr<Node> {
  return _children.at(c);
}

auto Node::getName() const -> std::string { return _label; }

void Node::setName(const std::string &s) { _label = s; }

auto Node::getNewick() const -> std::string {
  static auto newick_lambda = [](const Node &n) { return n.getName(); };

  return getNewickLambda(newick_lambda);
}

auto Node::getNewickLambda(const std::function<std::string(const Node &)>
                               &newick_lambda) const -> std::string {
  std::ostringstream newick_oss;
  for (size_t i = 0; i < getChildCount(); i++) {
    if (i == 0) { newick_oss << "("; }

    newick_oss << getChild(i)->getNewickLambda(newick_lambda);

    if (i == getChildCount() - 1) {
      newick_oss << ")";
    } else {
      newick_oss << ",";
    }
  }
  newick_oss << newick_lambda(*this);
  return newick_oss.str();
}

auto Node::getChildCount() const -> size_t { return _children.size(); }

auto Node::setHeight() -> double {
  double height = 0.0;
  for (const auto &c : _children) { height = std::max(height, c->setHeight()); }
  _height = height;
  return _height + _branch_length;
}

void Node::setHeightReverse(double h) {
  _height = h + _branch_length;
  for (const auto &c : _children) { c->setHeightReverse(_height); }
}

void getNodesByMRCAEntry(const std::shared_ptr<Node> &current,
                         const std::shared_ptr<MRCAEntry> &mrca,
                         std::vector<std::shared_ptr<Node>> &nodes) {
  /* if it is a tip, then check for membership */
  if (current->isTip()) {
    if (mrca->in(current->getName())) { nodes.push_back(current); }
    return;
  }

  for (const auto &c : current->_children) {
    getNodesByMRCAEntry(c, mrca, nodes);
  }
}

auto getMRCAWithNodes(const std::shared_ptr<Node> &current,
                      const std::vector<std::shared_ptr<Node>> &leaves)
    -> std::shared_ptr<Node> {
  if (current->_children.empty()) {
    if (std::any_of(leaves.begin(), leaves.end(), [&current](auto &n) {
          return n == current;
        })) {
      return current;
    }
    return {nullptr};
  }
  std::shared_ptr<Node> mrca = nullptr;
  for (auto &c : current->_children) {
    auto tmp_mrca = getMRCAWithNodes(c, leaves);
    if (tmp_mrca != nullptr) {
      // if this is our second match, then we can return with this pointer.
      if (mrca != nullptr) { return current; }
      mrca = tmp_mrca;
    }
  }
  return mrca;
}

auto getNodesByMRCALabel(const std::shared_ptr<Node> &current,
                         const MRCALabel &mrca) -> std::shared_ptr<Node> {
  if (current->getMRCALabel() == mrca) { return current; }
  for (const auto &c : current->_children) {
    auto tmp = getNodesByMRCALabel(c, mrca);
    if (tmp) { return tmp; }
  }
  return {nullptr};
}

auto getNodeById(const std::shared_ptr<Node> &current,
                 size_t id) -> std::shared_ptr<Node> {
  if (current->getId() == id) { return current; }
  for (const auto &c : current->_children) {
    auto tmp = getNodeById(c, id);
    if (tmp) { return tmp; }
  }
  return {nullptr};
}

auto Node::findNode(const std::shared_ptr<Node> &n) -> bool {
  if (this == n.get()) { return true; }
  return std::any_of(_children.begin(), _children.end(), [&n](auto &c) {
    return c->findNode(n);
  });
}

auto getParentWithNode(const std::shared_ptr<Node> &current,
                       const std::shared_ptr<Node> &n)
    -> std::shared_ptr<Node> {
  for (const auto &c : current->_children) {
    if (n == c) { return current; }
    auto ret = getParentWithNode(c, n);
    if (ret != nullptr) { return ret; }
  }
  return {nullptr};
}

auto Node::traverseAndGenerateForwardOperations(Workspace &ws,
                                                PeriodRateMatrixMap &pm_map,
                                                BranchProbMatrixMap &bm_map)
    const -> std::pair<std::vector<std::shared_ptr<SplitOperation>>,
                       std::shared_ptr<DispersionOperation>> {
  if (_children.size() != 2 && !_children.empty()) {
    throw std::runtime_error{
        "Tree is not bifircating when generating operations"};
  }

  if (_children.empty()) {
    ws.registerTopCLV(_id);
    return {{}, generateDispersionOperations(ws, pm_map, bm_map)};
  }

  std::vector<std::shared_ptr<SplitOperation>> split_ops;

  auto lchild =
      _children[0]->traverseAndGenerateForwardOperations(ws, pm_map, bm_map);
  auto rchild =
      _children[1]->traverseAndGenerateForwardOperations(ws, pm_map, bm_map);

  split_ops.reserve(lchild.first.size() + rchild.first.size() + 1);
  split_ops.insert(split_ops.end(), lchild.first.begin(), lchild.first.end());
  split_ops.insert(split_ops.end(), rchild.first.begin(), rchild.first.end());

  ws.registerTopCLV(_id);
  ws.registerChildrenCLV(_id);
  lchild.second->terminateTop(ws.getLeftChildCLV(_id));
  rchild.second->terminateTop(ws.getRightChildCLV(_id));

  split_ops.push_back(std::make_shared<SplitOperation>(
      ws.getTopCLV(_id), lchild.second, rchild.second));

  if (_incl_area_mask.hasValue()) {
    split_ops.back()->setInclAreas(_incl_area_mask.get());
  }

  if (_excl_area_mask.hasValue()) {
    split_ops.back()->setExclAreas(_excl_area_mask.get());
  }

  return {split_ops, generateDispersionOperations(ws, pm_map, bm_map)};
}

auto Node::traverseAndGenerateBackwardOperations(Workspace &ws,
                                                 PeriodRateMatrixMap &rm_map,
                                                 BranchProbMatrixMap &pm_map,
                                                 bool root) const
    -> std::pair<std::vector<std::shared_ptr<ReverseSplitOperation>>,
                 std::shared_ptr<DispersionOperation>> {
  if (_children.size() != 2 && !_children.empty()) {
    throw std::runtime_error{
        "Tree is not bifircating when generating operations"};
  }

  if (_children.empty()) { return {{}, {}}; }

  std::vector<std::shared_ptr<ReverseSplitOperation>> rsplit_ops;

  ws.registerTopCLVReverse(_id);

  auto disp_ops = generateDispersionOperationsReverse(ws, rm_map, pm_map);

  if (root && disp_ops == nullptr) {
    ws.registerBot1CLVReverse(_id);
    rsplit_ops.push_back(
        std::make_shared<ReverseSplitOperation>(ws.getBot1CLVReverse(_id),
                                                ws.getTopCLVReverse(_id),
                                                ws.getRightChildCLV(_id)));

    ws.registerBot2CLVReverse(_id);
    rsplit_ops.push_back(
        std::make_shared<ReverseSplitOperation>(ws.getBot2CLVReverse(_id),
                                                ws.getTopCLVReverse(_id),
                                                ws.getLeftChildCLV(_id)));
  } else {
    ws.registerBot1CLVReverse(_id);
    rsplit_ops.push_back(std::make_shared<ReverseSplitOperation>(
        ws.getBot1CLVReverse(_id), ws.getRightChildCLV(_id), disp_ops));

    ws.registerBot2CLVReverse(_id);
    rsplit_ops.push_back(std::make_shared<ReverseSplitOperation>(
        ws.getBot2CLVReverse(_id), ws.getLeftChildCLV(_id), disp_ops));
  }
  if (_incl_area_mask.hasValue()) {
    rsplit_ops[rsplit_ops.size() - 1]->setInclAreas(_incl_area_mask.get());
    rsplit_ops[rsplit_ops.size() - 2]->setInclAreas(_incl_area_mask.get());
  }

  if (_excl_area_mask.hasValue()) {
    rsplit_ops[rsplit_ops.size() - 1]->setExclAreas(_excl_area_mask.get());
    rsplit_ops[rsplit_ops.size() - 2]->setExclAreas(_excl_area_mask.get());
  }

  if (_children[0]->isInternal()) {
    auto child_trav =
        _children[0]->traverseAndGenerateBackwardOperations(ws, rm_map, pm_map);
    auto child_disp_op = child_trav.second;

    child_disp_op->terminateBot(ws.getBot1CLVReverse(_id));
    rsplit_ops.insert(
        rsplit_ops.end(), child_trav.first.begin(), child_trav.first.end());
  }

  if (_children[1]->isInternal()) {
    auto child_trav =
        _children[1]->traverseAndGenerateBackwardOperations(ws, rm_map, pm_map);
    auto child_disp_op = child_trav.second;

    child_disp_op->terminateBot(ws.getBot2CLVReverse(_id));
    rsplit_ops.insert(
        rsplit_ops.end(), child_trav.first.begin(), child_trav.first.end());
  }
  return {rsplit_ops, disp_ops};
}

auto Node::generateDispersionOperations(Workspace &ws,
                                        PeriodRateMatrixMap &rm_map,
                                        BranchProbMatrixMap &pm_map) const
    -> std::shared_ptr<DispersionOperation> {
  std::shared_ptr<DispersionOperation> child_op = nullptr;
  size_t top_clv = ws.getTopCLV(_id);

  for (const auto period : _periods) {
    auto tmp = std::make_shared<DispersionOperation>(
        top_clv, getProbMatrixOperation(ws, rm_map, pm_map, period));
    tmp->setChildOp(child_op);
    tmp->terminateTop(ws.registerGenericCLV());
    top_clv = tmp->topCLVIndex();
    child_op = tmp;
  }

  if (child_op) { child_op->unterminateTop(); }
  return child_op;
}

auto Node::generateDispersionOperationsReverse(
    Workspace &ws,
    PeriodRateMatrixMap &rm_map,
    BranchProbMatrixMap &pm_map) const -> std::shared_ptr<DispersionOperation> {
  std::vector<PeriodSegment> reverse_periods;
  for (auto p : _periods) { reverse_periods.push_back(p); }
  std::reverse(reverse_periods.begin(), reverse_periods.end());

  std::shared_ptr<DispersionOperation> child_op = nullptr;
  for (const auto period : reverse_periods) {
    auto tmp = std::make_shared<DispersionOperation>(
        ws.getTopCLVReverse(_id),
        std::numeric_limits<size_t>::max(),
        getProbMatrixOperation(ws, rm_map, pm_map, period, true));
    tmp->setChildOp(child_op);
    tmp->terminateBot(ws.getTopCLVReverse(_id));
    child_op = tmp;
  };

  if (child_op) { child_op->unterminateBot(); }
  return child_op;
}

void Node::traverseAndGenerateBackwardNodeIdsInternalOnly(
    std::vector<size_t> &ret) const {
  ret.push_back(_id);
  for (const auto &c : _children) {
    if (c->isInternal()) {
      c->traverseAndGenerateBackwardNodeIdsInternalOnly(ret);
    }
  }
}

void Node::traverseAndGeneratePostorderNodeIdsInternalOnly(
    std::vector<size_t> &ret) const {
  for (const auto &c : _children) {
    if (c->isInternal()) {
      c->traverseAndGeneratePostorderNodeIdsInternalOnly(ret);
    }
  }
  ret.push_back(_id);
}

void Node::traverseAndGenerateBackwardNodeNumbersInternalOnly(
    std::vector<size_t> &ret) const {
  ret.push_back(getNumber());
  for (const auto &c : _children) {
    if (c->isInternal()) {
      c->traverseAndGenerateBackwardNodeNumbersInternalOnly(ret);
    }
  }
}

void Node::assignTipData(
    Workspace &ws,
    const std::unordered_map<std::string, Range> &distrib_data) const {
  if (_children.empty()) {
    ws.setTipCLV(ws.getTopCLV(_id), distrib_data.at(_label));
  } else {
    for (const auto &c : _children) { c->assignTipData(ws, distrib_data); }
  }
}

void Node::assignIdRecursive(size_t &id) {
  _id = id++;
  for (auto &c : _children) { c->assignIdRecursive(id); }
}

void Node::assignId() {
  size_t id = 0;
  assignIdRecursive(id);
}

auto Node::getId() const -> size_t { return _id; }

auto Node::checkAlignmentConsistency(const Alignment &align,
                                     size_t count) -> size_t {
  if (isTip()) {
    if (align.data.count(_label) == 0) {
      std::ostringstream oss;
      oss << "Could not find taxa '" << _label
          << "' in alignment file, but taxa is present on the tree";
      throw std::runtime_error{oss.str()};
    }
    return count + 1;
  }
  for (const auto &c : _children) {
    count = c->checkAlignmentConsistency(align, count);
  }
  return count;
}

void Node::assignIncludedAreas(Range incl_area_mask) {
  _incl_area_mask = incl_area_mask;
}

void Node::assignExcludedAreas(Range excl_area_mask) {
  _excl_area_mask = excl_area_mask;
}

void Node::assignFixedDist(Range fixed_dist) { _fixed_dist = fixed_dist; }

void Node::applyPreorderInternalOnly(const std::function<void(Node &)> &func) {
  applyCB(func);
  for (auto &c : _children) {
    if (c->isInternal()) { c->applyPreorderInternalOnly(func); }
  }
}

void Node::applyPreorder(const std::function<void(Node &)> &func) {
  applyCB(func);
  for (auto &c : _children) { c->applyPreorder(func); }
}

auto Node::getFixedDist() const -> Option<Range> { return _fixed_dist; }

auto Node::getIncludedAreas() const -> Option<Range> { return _incl_area_mask; }

auto Node::getExcludedAreas() const -> Option<Range> { return _excl_area_mask; }

void Node::setPeriodSegments(const Periods &periods) {
  _periods = PeriodSpan(periods, _height, _branch_length);
}

auto Node::getRateMatrixOperation(Workspace &ws,
                                  PeriodRateMatrixMap &rm_map,
                                  size_t period) const
    -> std::shared_ptr<MakeRateMatrixOperation> {
  auto it = rm_map.find(period);
  if (it == rm_map.end()) {
    auto rm = std::make_shared<MakeRateMatrixOperation>(
        ws.reserveRateMatrixIndex(period), period);
    rm_map.emplace(period, rm);
    return rm;
  }
  if (it->second == nullptr) {
    throw std::runtime_error{"We got an empty expm"};
  }
  return it->second;
}

auto Node::getProbMatrixOperation(Workspace &ws,
                                  PeriodRateMatrixMap &rm_map,
                                  BranchProbMatrixMap &pm_map,
                                  PeriodSegment period,
                                  bool transpose) const
    -> std::shared_ptr<ExpmOperation> {
  auto key = std::make_pair(period.index, period.duration);
  auto it = pm_map.find(key);
  if (it == pm_map.end()) {
    auto rm = getRateMatrixOperation(ws, rm_map, period.index);
    auto pm = std::make_shared<ExpmOperation>(
        ws.suggestProbMatrixIndex(), _branch_length, rm, transpose);
    pm_map.emplace(key, pm);
    return pm;
  }
  if (it->second == nullptr) {
    throw std::runtime_error{"We got an empty expm"};
  }
  return it->second;
}

void Node::applyCB(const std::function<void(Node &)> &func) { func(*this); }

auto Node::hasResults() const -> bool {
  return hasAncestralState() || hasAncestralSplit();
}

auto Node::hasAncestralState() const -> bool {
  return _ancestral_state.has_value();
}

auto Node::hasAncestralSplit() const -> bool {
  return _ancestral_split.has_value();
}

auto Node::getAncestralState() const
    -> const std::unique_ptr<LagrangeMatrixBase[]> & {
  return _ancestral_state.value();
}

auto Node::getAncestralState() -> std::unique_ptr<LagrangeMatrixBase[]> & {
  return _ancestral_state.value();
}

auto Node::getAncestralSplit() const -> const SplitReturn & {
  return _ancestral_split.value();
}

auto Node::getAncestralSplit() -> SplitReturn & {
  return _ancestral_split.value();
}

auto Node::getTopAncestralSplit() const -> AncSplit {
  if (_ancestral_split) {
    auto best_state = getTopAncestralState();
    auto split = _ancestral_split.value().at(best_state);
    auto best_split = *std::max_element(
        split.begin(), split.end(), [](const auto &a, const auto &b) {
          return a.getLikelihood() < b.getLikelihood();
        });
    return best_split;
  }

  throw std::runtime_error{"No ancestral splits assigned"};
}

auto Node::getTopAncestralState() const -> Range {
  if (_ancestral_state) {
    auto &state = _ancestral_state.value();
    auto period_segment = *_periods.begin();

    Range cur_range;
    size_t cur_index;
    Range best_range = 0;
    double best_lh = -std::numeric_limits<double>::infinity();
    size_t max_range = 1ul << period_segment.regions;
    for (cur_range = 0, cur_index = 0; cur_range < max_range; ++cur_index,
        cur_range = next_dist(cur_range, period_segment.max_areas)) {
      double cur_lwr = state[cur_index];
      if (best_lh < cur_lwr) {
        best_lh = cur_lwr;
        best_range = cur_range;
      }
    }

    return best_range;
  }
  throw std::runtime_error{"No ancestral states assigned"};
}

void Node::assignAncestralState(std::unique_ptr<LagrangeMatrixBase[]> s) {
  _ancestral_state = std::move(s);
}

void Node::assignAncestralSplit(SplitReturn s) { _ancestral_split = s; }

auto Node::getCount() -> size_t { return getCount(0); }

auto Node::getTipCount() -> size_t { return getTipCount(0); }

auto Node::getCount(size_t cur) -> size_t {
  for (const auto &c : _children) { cur = c->getCount(cur); }
  cur += 1;
  return cur;
}

auto Node::getTipCount(size_t cur) -> size_t {
  if (isTip()) { return cur + 1; }
  for (const auto &c : _children) { cur = c->getTipCount(cur); }
  return cur;
}

void Node::setMRCALabel(const MRCALabel &l) { _mrca = l; }

auto Node::getMRCALabel() const -> MRCALabel { return _mrca; }

auto Node::assignFossil(const Fossil &f) -> bool {
  if (getMRCALabel() == f.mrca_name) {
    switch (f.type) {
      case FossilType::NODE:
      case FossilType::INCLUDE:
        assignIncludedAreas(f.area);
        break;
      case FossilType::BRANCH:
      case FossilType::FIXED:
        assignFixedDist(f.area);
        break;
      case FossilType::EXCLUDE:
        assignExcludedAreas(f.area);
        break;
      case FossilType::UNKOWN:
        throw std::runtime_error{
            "unknown fossil type when assigning fossiles to node"};
    }
    return true;
  }
  for (auto &c : _children) {
    if (c->assignFossil(f)) { return true; }
  }
  return false;
}

auto Node::getNodeLabel() const -> std::string {
  if (!_mrca.empty()) { return _mrca; }
  if (isTip()) { return getName(); }
  return std::to_string(getNumber());
}

}  // namespace lagrange
