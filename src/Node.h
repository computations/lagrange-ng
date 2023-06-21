/*
 * node.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef NODE_H
#define NODE_H

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Alignment.h"
#include "Common.h"
#include "Operation.h"
#include "Periods.hpp"
#include "Utils.h"
#include "Workspace.h"

class Node {
 private:
  void assignIdRecursive(size_t &id);

  double _branch_length;  // branch length
  double _height;         // could be from tip or from root

  size_t _number;
  size_t _id;

  PeriodSpan _periods;

  std::string _label;
  std::string _comment;
  std::vector<std::shared_ptr<Node>> _children;

  lagrange_option_t<lagrange_dist_t> _fixed_dist;
  lagrange_option_t<lagrange_dist_t> _incl_area_mask;
  lagrange_option_t<lagrange_dist_t> _excl_area_mask;

  auto getRateMatrixOperation(Workspace &ws, PeriodRateMatrixMap &rm_map,
                              size_t period) const
      -> std::shared_ptr<MakeRateMatrixOperation>;

  auto getProbMatrixOperation(Workspace &ws, PeriodRateMatrixMap &rm_map,
                              BranchProbMatrixMap &pm_map, PeriodSegment period,
                              bool transpose = false) const
      -> std::shared_ptr<ExpmOperation>;

 public:
  Node();
  Node(double bl, size_t number, std::string name);
  ~Node() {
    for (auto &c : _children) { c.reset(); }
  }

  auto isExternal() const -> bool;
  auto isInternal() const -> bool;

  auto getNumber() const -> size_t;
  void setNumber(size_t n);

  auto getId() const -> size_t;

  auto getBL() const -> double;
  void setBL(double bl);

  auto hasChild(const std::shared_ptr<Node> &test) -> bool;
  auto addChild(const std::shared_ptr<Node> &c) -> bool;
  auto removeChild(const std::shared_ptr<Node> &c) -> bool;
  auto getChild(size_t c) const -> std::shared_ptr<Node>;

  auto getName() const -> std::string;
  void setName(const std::string &s);
  void setComment(const std::string &s);

  auto getNewick() const -> std::string;
  auto getNewickLambda(const std::function<std::string(const Node &)> &) const
      -> std::string;

  auto getChildCount() const -> size_t;

  void setSplitStringRecursive(
      const std::vector<size_t> &id_map,
      const std::vector<lagrange_col_vector_t> &dist_lhs, size_t states,
      const std::vector<std::string> &names);
  void setStateStringRecursive(
      const std::vector<size_t> &id_map,
      const std::vector<lagrange_col_vector_t> &dist_lhs, size_t states,
      const std::vector<std::string> &names);

  auto findNode(const std::shared_ptr<Node> &n) -> bool;

  auto setHeight() -> double;
  void setHeightReverse(double height = 0.0);

  friend auto getParentWithNode(const std::shared_ptr<Node> &current,
                                const std::shared_ptr<Node> &n)
      -> std::shared_ptr<Node>;

  friend auto getMRCAWithNode(const std::shared_ptr<Node> &current,
                              const std::vector<std::shared_ptr<Node>> &nodes)
      -> std::shared_ptr<Node>;

  auto traverseAndGenerateForwardOperations(Workspace &ws,
                                            PeriodRateMatrixMap &pm_map,
                                            BranchProbMatrixMap &bm_map) const
      -> std::pair<std::vector<std::shared_ptr<SplitOperation>>,
                   std::shared_ptr<DispersionOperation>>;

  auto traverseAndGenerateBackwardOperations(Workspace &ws,
                                             PeriodRateMatrixMap &rm_map,
                                             BranchProbMatrixMap &pm_map,
                                             bool root = false) const
      -> std::pair<std::vector<std::shared_ptr<ReverseSplitOperation>>,
                   std::shared_ptr<DispersionOperation>>;

  auto generateDispersionOperations(Workspace &ws, PeriodRateMatrixMap &rm_map,
                                    BranchProbMatrixMap &pm_map) const
      -> std::shared_ptr<DispersionOperation>;

  auto generateDispersionOperationsReverse(Workspace &ws,
                                           PeriodRateMatrixMap &rm_map,
                                           BranchProbMatrixMap &pm_map) const
      -> std::shared_ptr<DispersionOperation>;

  auto traverseAndGenerateBackwardOperations(
      Workspace &ws, const std::shared_ptr<MakeRateMatrixOperation> &rm_op,
      bool root = false) const
      -> std::pair<std::vector<ReverseSplitOperation>,
                   std::shared_ptr<DispersionOperation>>;

  void traverseAndGenerateBackwardNodeIdsInternalOnly(
      std::vector<size_t> &) const;

  void traverseAndGenerateBackwardNodeNumbersInternalOnly(
      std::vector<size_t> &ret) const;

  void traverseAndGeneratePostorderNodeIdsInternalOnly(
      std::vector<size_t> &ret) const;

  void applyPreorderInternalOnly(const std::function<void(Node &)> &func);
  void applyPreorder(const std::function<void(Node &)> &func);

  void assignTipData(Workspace &ws,
                     const std::unordered_map<std::string, lagrange_dist_t>
                         &distrib_data) const;

  size_t checkAlignmentConsistency(const Alignment &align, size_t count);

  void assignId();

  void assignInclAreas(lagrange_dist_t fixed_dist);
  void assignFixedDist(lagrange_dist_t fixed_dist);

  lagrange_option_t<lagrange_dist_t> getFixedDist() const;

  lagrange_option_t<lagrange_dist_t> getInclAreas() const;

  void setPeriodSegments(const Periods &periods);
};

auto getMRCAWithNode(const std::shared_ptr<Node> &current,
                     const std::vector<std::shared_ptr<Node>> &nodes)
    -> std::shared_ptr<Node>;

#endif /* NODE_H_ */
