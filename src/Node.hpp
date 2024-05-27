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
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Alignment.hpp"
#include "AncSplit.hpp"
#include "Common.hpp"
#include "Fossil.hpp"
#include "MRCA.hpp"
#include "Operation.hpp"
#include "Periods.hpp"
#include "Utils.hpp"
#include "Workspace.hpp"

namespace lagrange {

class Node {
 public:
  Node();
  Node(double bl, std::string name);

  ~Node() {
    for (auto &c : _children) { c.reset(); }
  }

  auto isTip() const -> bool;
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

  void setMRCALabel(const MRCALabel &);
  auto getMRCALabel() const -> MRCALabel;

  auto getNewick() const -> std::string;
  auto getNewickLambda(const std::function<std::string(const Node &)> &) const
      -> std::string;

  auto getChildCount() const -> size_t;

  auto getCount() -> size_t;
  auto getTipCount() -> size_t;

  auto findNode(const std::shared_ptr<Node> &n) -> bool;

  auto setHeight() -> double;
  void setHeightReverse(double height = 0.0);

  friend auto getParentWithNode(const std::shared_ptr<Node> &current,
                                const std::shared_ptr<Node> &n)
      -> std::shared_ptr<Node>;

  friend void getNodesByMRCAEntry(const std::shared_ptr<Node> &current,
                                  const std::shared_ptr<MRCAEntry> &mrca,
                                  std::vector<std::shared_ptr<Node>> &nodes);

  friend auto getNodesByMRCALabel(const std::shared_ptr<Node> &current,
                                  const MRCALabel &mrca)
      -> std::shared_ptr<Node>;

  friend auto getMRCAWithNodes(const std::shared_ptr<Node> &current,
                               const std::vector<std::shared_ptr<Node>> &nodes)
      -> std::shared_ptr<Node>;

  friend auto getNodeById(const std::shared_ptr<Node> &current, size_t id)
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

  auto generateDispersionOperations(Workspace &ws,
                                    PeriodRateMatrixMap &rm_map,
                                    BranchProbMatrixMap &pm_map) const
      -> std::shared_ptr<DispersionOperation>;

  auto generateDispersionOperationsReverse(Workspace &ws,
                                           PeriodRateMatrixMap &rm_map,
                                           BranchProbMatrixMap &pm_map) const
      -> std::shared_ptr<DispersionOperation>;

  auto traverseAndGenerateBackwardOperations(
      Workspace &ws,
      const std::shared_ptr<MakeRateMatrixOperation> &rm_op,
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

  void assignTipData(
      Workspace &ws,
      const std::unordered_map<std::string, Dist> &distrib_data) const;

  size_t checkAlignmentConsistency(const Alignment &align, size_t count);

  void assignId();
  std::string getNodeLabel() const;

  void assignIncludedAreas(Dist fixed_dist);
  void assignExcludedAreas(Dist fixed_dist);
  void assignFixedDist(Dist fixed_dist);

  Option<Dist> getFixedDist() const;

  Option<Dist> getIncludedAreas() const;
  Option<Dist> getExcludedAreas() const;

  void setPeriodSegments(const Periods &periods);

  void applyCB(const std::function<void(Node &)> &func);

  bool hasResults() const;
  bool hasAncestralState() const;
  bool hasAncestralSplit() const;

  std::unique_ptr<LagrangeMatrixBase[]> &getAncestralState();
  SplitReturn &getAncestralSplit();

  void assignAncestralState(std::unique_ptr<LagrangeMatrixBase[]>);
  void assignAncestralSplit(SplitReturn);

  bool assignFossil(const Fossil &);

 private:
  auto getCount(size_t) -> size_t;
  auto getTipCount(size_t) -> size_t;
  void assignIdRecursive(size_t &id);

  auto getRateMatrixOperation(Workspace &ws,
                              PeriodRateMatrixMap &rm_map,
                              size_t period) const
      -> std::shared_ptr<MakeRateMatrixOperation>;

  auto getProbMatrixOperation(Workspace &ws,
                              PeriodRateMatrixMap &rm_map,
                              BranchProbMatrixMap &pm_map,
                              PeriodSegment period,
                              bool transpose = false) const
      -> std::shared_ptr<ExpmOperation>;

  double _branch_length;  // branch length
  double _height;         // could be from tip or from root

  size_t _id;

  PeriodSpan _periods;

  std::string _label;
  MRCALabel _mrca;
  std::vector<std::shared_ptr<Node>> _children;

  Option<Dist> _fixed_dist;
  Option<Dist> _incl_area_mask;
  Option<Dist> _excl_area_mask;

  std::optional<std::unique_ptr<LagrangeMatrixBase[]>> _ancestral_state;
  std::optional<SplitReturn> _ancestral_split;
};

auto getMRCAWithNodes(const std::shared_ptr<Node> &current,
                      const std::vector<std::shared_ptr<Node>> &nodes)
    -> std::shared_ptr<Node>;

void getNodesByMRCAEntry(const std::shared_ptr<Node> &current,
                         const std::shared_ptr<MRCAEntry> &mrca,
                         std::vector<std::shared_ptr<Node>> &nodes);

auto getNodesByMRCALabel(const std::shared_ptr<Node> &current,
                         const MRCALabel &mrca) -> std::shared_ptr<Node>;

auto getNodeById(const std::shared_ptr<Node> &current, size_t id)
    -> std::shared_ptr<Node>;
}  // namespace lagrange
#endif /* NODE_H_ */
