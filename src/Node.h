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

#include "Common.h"
#include "Operation.h"
#include "Utils.h"
#include "Workspace.h"

class Node {
 private:
  void assignIdRecursive(size_t &id);

  double _branch_length;  // branch lengths
  double _height;         // could be from tip or from root
  size_t _number;
  size_t _id;
  size_t _period;
  std::string _label;
  std::string _comment;
  std::string _split_string;
  std::string _state_string;
  std::string _stoch_string;
  std::vector<std::shared_ptr<Node>> _children;
  std::shared_ptr<std::vector<lagrange_dist_t>> _excluded_dists;

  auto getRateMatrixOperation(PeriodRateMatrixMap &rm_map) const
      -> std::shared_ptr<MakeRateMatrixOperation> {
    auto it = rm_map.find(_period);
    if (it == rm_map.end()) {
      auto rm = std::make_shared<MakeRateMatrixOperation>(
          Workspace::suggest_rate_matrix_index());
      rm_map.emplace(_period, rm);
      return rm;
    }
    if (it->second == nullptr) {
      throw std::runtime_error{"We got an empty expm"};
    }
    return it->second;
  }

  auto getProbMatrixOperation(Workspace &ws, PeriodRateMatrixMap &rm_map,
                              BranchProbMatrixMap &pm_map,
                              bool transpose = false) const
      -> std::shared_ptr<ExpmOperation> {
    auto key = std::make_pair(_period, _branch_length);
    auto it = pm_map.find(key);
    if (it == pm_map.end()) {
      auto rm = getRateMatrixOperation(rm_map);
      auto pm = std::make_shared<ExpmOperation>(ws.suggest_prob_matrix_index(),
                                                _branch_length, rm, transpose);
      pm_map.emplace(key, pm);
      return pm;
    }
    if (it->second == nullptr) {
      throw std::runtime_error{"We got an empty expm"};
    }
    return it->second;
  }

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

  void setHeight(double he);

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

  void setSplitString(const std::string &splitstring);
  void setStateString(const std::string &splitstring);
  void setStochString(const std::string &stochstring);
  auto getSplitString() const -> std::string;
  auto getStateString() const -> std::string;
  auto getStochString() const -> std::string;

  void initExclDistVector();

  auto findNode(const std::shared_ptr<Node> &n) -> bool;

  auto getMaxHeightRecursive() const -> double;
  auto getMaxHeight() const -> double;
  void setHeightRecursive(double height);
  void setHeightRecursive();

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
                                             BranchProbMatrixMap &pm_map) const
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
      Workspace &ws,
      const std::shared_ptr<MakeRateMatrixOperation> &rm_op) const
      -> std::pair<std::vector<ReverseSplitOperation>,
                   std::shared_ptr<DispersionOperation>>;

  void traverseAndGenerateBackwardNodeIdsInternalOnly(
      std::vector<size_t> &) const;

  void traverseAndGenerateBackwardNodeNumbersInternalOnly(
      std::vector<size_t> &ret) const;

  void traverseAndGeneratePostorderNodeIdsInternalOnly(
      std::vector<size_t> &ret) const;

  void assignTipData(Workspace &ws,
                     const std::unordered_map<std::string, lagrange_dist_t>
                         &distrib_data) const;

  void assignId();
};

auto getMRCAWithNode(const std::shared_ptr<Node> &current,
                     const std::vector<std::shared_ptr<Node>> &nodes)
    -> std::shared_ptr<Node>;

#endif /* NODE_H_ */
