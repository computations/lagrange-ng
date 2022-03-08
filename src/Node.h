/*
 * node.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#ifndef NODE_H_
#define NODE_H_

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

  std::shared_ptr<MakeRateMatrixOperation> getRateMatrixOperation(
      Workspace &ws, PeriodRateMatrixMap &rm_map) const {
    auto it = rm_map.find(_period);
    if (it == rm_map.end()) {
      auto rm = std::make_shared<MakeRateMatrixOperation>(
          ws.suggest_rate_matrix_index());
      rm_map.emplace(_period, rm);
      return rm;
    } else {
      if (it->second == nullptr) {
        throw std::runtime_error{"We got an empty expm"};
      }
      return it->second;
    }
  }

  std::shared_ptr<ExpmOperation> getProbMatrixOperation(
      Workspace &ws, PeriodRateMatrixMap &rm_map, BranchProbMatrixMap &pm_map,
      bool transpose = false) const {
    auto key = std::make_pair(_period, _branch_length);
    auto it = pm_map.find(key);
    if (it == pm_map.end()) {
      auto rm = getRateMatrixOperation(ws, rm_map);
      auto pm = std::make_shared<ExpmOperation>(ws.suggest_prob_matrix_index(),
                                                _branch_length, rm, transpose);
      pm_map.emplace(key, pm);
      return pm;
    } else {
      if (it->second == nullptr) {
        throw std::runtime_error{"We got an empty expm"};
      }
      return it->second;
    }
  }

 public:
  Node();
  Node(double bl, size_t number, const std::string &name);
  ~Node() {
    for (auto &c : _children) { c.reset(); }
  }

  std::vector<std::shared_ptr<Node>> getChildren();
  bool isExternal() const;
  bool isInternal() const;

  size_t getNumber() const;
  void setNumber(size_t n);

  size_t getId() const;

  double getBL();
  void setBL(double bl);

  double getHeight();
  void setHeight(double he);

  bool hasChild(std::shared_ptr<Node> test);
  bool addChild(std::shared_ptr<Node> c);
  bool removeChild(std::shared_ptr<Node> c);
  std::shared_ptr<Node> getChild(size_t c) const;

  std::string getName() const;
  std::string getComment() const;
  void setName(const std::string &s);
  void setComment(const std::string &s);

  std::string getNewick() const;
  std::string getNewickLambda(
      const std::function<std::string(const Node &)> &) const;

  size_t getChildCount() const;

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
  std::string getSplitString() const;
  std::string getStateString() const;
  std::string getStochString() const;

  void initExclDistVector();
  std::shared_ptr<std::vector<lagrange_dist_t>> &getExclDistVector();

  bool findNode(std::shared_ptr<Node> n);

  double getMaxHeightRecursive() const;
  double getMaxHeight() const;
  void setHeightRecursive(double height);
  void setHeightRecursive();
  void pruneNode(std::shared_ptr<Node> n);

  friend std::shared_ptr<Node> getParentWithNode(
      const std::shared_ptr<Node> &current, const std::shared_ptr<Node> &n);

  friend std::shared_ptr<Node> getMRCAWithNode(
      const std::shared_ptr<Node> &current,
      const std::vector<std::shared_ptr<Node>> &nodes);

  std::pair<std::vector<std::shared_ptr<SplitOperation>>,
            std::shared_ptr<DispersionOperation>>
  traverseAndGenerateForwardOperations(Workspace &ws,
                                       PeriodRateMatrixMap &pm_map,
                                       BranchProbMatrixMap &bm_map) const;

  std::pair<std::vector<std::shared_ptr<ReverseSplitOperation>>,
            std::shared_ptr<DispersionOperation>>
  traverseAndGenerateBackwardOperations(Workspace &ws,
                                        PeriodRateMatrixMap &rm_map,
                                        BranchProbMatrixMap &pm_map) const;

  std::shared_ptr<DispersionOperation> generateDispersionOperations(
      Workspace &ws, PeriodRateMatrixMap &rm_map,
      BranchProbMatrixMap &pm_map) const;

  std::shared_ptr<DispersionOperation> generateDispersionOperationsReverse(
      Workspace &ws, PeriodRateMatrixMap &rm_map,
      BranchProbMatrixMap &pm_map) const;

  std::pair<std::vector<ReverseSplitOperation>,
            std::shared_ptr<DispersionOperation>>
  traverseAndGenerateBackwardOperations(
      Workspace &ws,
      const std::shared_ptr<MakeRateMatrixOperation> &rm_op) const;

  void traverseAndGenerateBackwardNodeIds(std::vector<size_t> &) const;
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

std::shared_ptr<Node> getMRCAWithNode(
    const std::shared_ptr<Node> &current,
    const std::vector<std::shared_ptr<Node>> &nodes);

#endif /* NODE_H_ */
