#include "IO.hpp"

#include <fstream>
#include <memory>
#include <optional>

#include "AncSplit.hpp"
#include "MRCA.hpp"
#include "Tree.hpp"
#include "Utils.hpp"

#define CONVERT_FLOAT_TO_JSON(y, x) \
  if (std::isfinite(x)) {           \
    y = x;                          \
  } else {                          \
    y = std::to_string(x);          \
  }

namespace lagrange {

auto normalize_split_distribution_by_lwr(SplitReturn &splits) -> void {
  double max_llh = -std::numeric_limits<double>::infinity();
  for (const auto &kv : splits) {
    for (const auto &sp : kv.second) {
      max_llh = std::max(max_llh, sp.getLikelihood());
    }
  }

  double total_llh = 0.0;
  for (const auto &kv : splits) {
    for (const auto &sp : kv.second) {
      total_llh += std::exp(sp.getLikelihood() - max_llh);
    }
  }

  for (auto &kv : splits) {
    for (auto &sp : kv.second) {
      sp.setLWR(std::exp(sp.getLikelihood() - max_llh) / total_llh);
    }
  }
}

auto normalize_state_distribution_by_lwr(
    const std::unique_ptr<LagrangeMatrixBase[]> &states, size_t states_len)
    -> std::unique_ptr<LagrangeMatrixBase[]> {
  std::unique_ptr<LagrangeMatrixBase[]> normalized_states{
      new LagrangeMatrixBase[states_len]};

  for (size_t i = 0; i < states_len; i++) {
    normalized_states.get()[i] = states.get()[i];
  }

  double max_llh = -std::numeric_limits<double>::infinity();

  assert(states_len != 1);

  for (size_t i = 0; i < states_len; ++i) {
    double tmp = normalized_states.get()[i];

    max_llh = std::max(max_llh, tmp);
  }

  double total_llh = 0.0;
  for (size_t i = 1; i < states_len; i++) {
    total_llh += std::exp(normalized_states.get()[i] - max_llh);
  }

  for (size_t i = 0; i < states_len; i++) {
    double tmp = std::exp(normalized_states.get()[i] - max_llh) / total_llh;

    normalized_states.get()[i] = tmp;
  }

  return normalized_states;
}

auto make_state_results_for_node(
    const std::unique_ptr<LagrangeMatrixBase[]> &state_distribution,
    const std::vector<std::string> &region_names,
    size_t states_len,
    size_t max_areas) -> nlohmann::json {
  nlohmann::json node_json;
  auto lwr_distribution =
      normalize_state_distribution_by_lwr(state_distribution, states_len);
  for (size_t dist = 0, dist_index = 0; dist_index < states_len;
       ++dist_index, dist = next_dist(dist, max_areas)) {
    nlohmann::json tmp;
    tmp["distribution"] = dist;
    double llh = state_distribution.get()[dist_index];
    tmp["distribution-string"] =
        lagrange_convert_dist_string(dist, region_names);
    CONVERT_FLOAT_TO_JSON(tmp["llh"], llh);
    double ratio = lwr_distribution.get()[dist_index];
    CONVERT_FLOAT_TO_JSON(tmp["ratio"], llh);
    tmp["regions"] = lagrange_convert_dist_to_list(dist, region_names);
    node_json.push_back(tmp);
  }
  return node_json;
}

auto make_split_results_for_node(SplitReturn &splits,
                                 const std::vector<std::string> &region_names,
                                 size_t states_len,
                                 size_t max_areas) -> nlohmann::json {
  nlohmann::json node_json;
  normalize_split_distribution_by_lwr(splits);
  for (size_t dist = 0, dist_index = 0; dist_index < states_len;
       ++dist_index, dist = next_dist(dist, max_areas)) {
    for (const auto &sp : splits[dist]) {
      nlohmann::json anc_json;
      nlohmann::json left_json;
      nlohmann::json right_json;

      assert(dist == sp.anc_dist);
      anc_json["distribution"] = dist;
      anc_json["distribution-string"] =
          lagrange_convert_dist_string(dist, region_names);
      anc_json["regions"] = lagrange_convert_dist_to_list(dist, region_names);

      left_json["distribution"] = sp.l_dist;
      left_json["distribution-string"] =
          lagrange_convert_dist_string(sp.l_dist, region_names);
      left_json["regions"] =
          lagrange_convert_dist_to_list(sp.l_dist, region_names);

      right_json["distribution"] = sp.r_dist;
      right_json["distribution-string"] =
          lagrange_convert_dist_string(sp.r_dist, region_names);
      right_json["regions"] =
          lagrange_convert_dist_to_list(sp.r_dist, region_names);

      nlohmann::json tmp;
      tmp["anc-dist"] = anc_json;
      tmp["left-dist"] = left_json;
      tmp["right-dist"] = right_json;
      CONVERT_FLOAT_TO_JSON(tmp["llh"], sp.getLikelihood());
      CONVERT_FLOAT_TO_JSON(tmp["ratio"], sp.getLWR());
      node_json.push_back(tmp);
    }
  }

  return node_json;
}

auto make_results_for_nodes(const std::shared_ptr<Tree> &tree,
                            const std::vector<std::string> &region_names,
                            size_t states_len,
                            size_t max_areas) -> nlohmann::json {
  nlohmann::json root_json;

  auto cb = [&](Node &n) {
    nlohmann::json node_json;
    node_json["number"] = n.getNodeLabel();
    if (n.hasAncestralSplit()) {
      node_json["splits"] = make_split_results_for_node(
          n.getAncestralSplit(), region_names, states_len, max_areas);
    }
    if (n.hasAncestralState()) {
      node_json["states"] = make_state_results_for_node(
          n.getAncestralState(), region_names, states_len, max_areas);
    }
    root_json.push_back(node_json);
  };

  tree->applyPreorderInternalOnly(cb);

  return root_json;
}

void write_result_file(const std::shared_ptr<Tree> &tree,
                       const ConfigFile &config,
                       const Context &context) {
  auto root_json = init_json(tree, config);
  nlohmann::json params_json;
  auto params = context.currentParams();
  for (size_t i = 0; i < params.size(); ++i) {
    params_json[i]["dispersion"] = params[i].dispersion_rate;
    params_json[i]["extinction"] = params[i].extinction_rate;
  }
  root_json["params"] = params_json;

  auto states = context.getStateResults();
  auto splits = context.getSplitResults();
  root_json["node-results"] = make_results_for_nodes(
      tree, config.area_names, context.stateCount(), config.maxareas);

  write_json_file(config, root_json);
}

void write_node_tree(const std::shared_ptr<Tree> &tree,
                     const ConfigFile &config) {
  auto node_tree_filename = config.NodeTreeFilename();
  std::cout << "Writing node annotated tree to " << node_tree_filename
            << std::endl;

  std::ofstream node_tree(node_tree_filename);
  node_tree << tree->getNewickLambda([](const Node &n) -> std::string {
    return n.getNodeLabel() + ":" + std::to_string(n.getBL());
  });
}

void write_scaled_tree(const std::shared_ptr<Tree> &tree,
                       const ConfigFile &config) {
  auto scaled_tree_filename = config.scaledTreeFilename();
  std::cout << "Writing scaled tree to " << scaled_tree_filename << std::endl;
  std::ofstream anal_tree(scaled_tree_filename);
  anal_tree << tree->getNewickLambda([](const Node &n) -> std::string {
    return n.getName() + ":" + std::to_string(n.getBL());
  }) << std::endl;
}

auto init_json(const std::shared_ptr<const Tree> &tree,
               const ConfigFile &config) -> nlohmann::json {
  nlohmann::json root_json;
  nlohmann::json attributes_json;
  attributes_json["periods"] =
      !config.periods.empty() ? static_cast<int>(config.periods.size()) : 1;
  attributes_json["regions"] = config.region_count;
  attributes_json["taxa"] = tree->getTipCount();
  attributes_json["nodes-tree"] =
      tree->getNewickLambda([](const Node &n) -> std::string {
        if (n.isInternal()) {
          return n.getNodeLabel() + ":" + std::to_string(n.getBL());
        } else {
          return n.getName() + ":" + std::to_string(n.getBL());
        }
      });
  attributes_json["max-areas"] = config.maxareas;
  attributes_json["state-count"] = lagrange_compute_restricted_state_count(
      config.region_count, config.maxareas);
  root_json["attributes"] = attributes_json;

  return root_json;
}

void write_json_file(const ConfigFile &config,
                     const nlohmann ::json &root_json) {
  auto json_filename = config.resultsFilename();
  std::cout << "Writing results to " << json_filename << std::endl;
  std::ofstream outfile(config.resultsFilename());
  outfile << root_json.dump();
}
}  // namespace lagrange
