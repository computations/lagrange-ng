#include "IO.hpp"

#include <cmath>
#include <fstream>
#include <initializer_list>
#include <logger.hpp>
#include <memory>
#include <sstream>
#include <string>

#undef NDEBUG
#include <cassert>

#include "AncSplit.hpp"
#include "ConfigFile.hpp"
#include "Periods.hpp"
#include "Tree.hpp"
#include "Utils.hpp"

#define CONVERT_FLOAT_TO_JSON(y, x) \
  if (std::isfinite(x)) {           \
    y = x;                          \
  } else {                          \
    y = std::to_string(x);          \
  }

#define STRING(s) #s
#define STRINGIFY(s) STRING(s)
#define GIT_REV_STRING STRINGIFY(LAGRANGE_BUILD_VERSION)

namespace lagrange {

auto normalize_split_distribution_by_lwr(SplitReturn &splits) -> void {
  double max_llh = -std::numeric_limits<double>::infinity();
  for (const auto &kv : splits) {
    for (const auto &sp : kv.second) {
      max_llh = std::max(max_llh, sp.getLikelihood());
    }
  }
  assert(std::isfinite(max_llh));

  double total_llh = 0.0;
  for (const auto &kv : splits) {
    for (const auto &sp : kv.second) {
      total_llh += std::exp(sp.getLikelihood() - max_llh);
    }
  }
  assert(std::isfinite(total_llh));

  for (auto &kv : splits) {
    for (auto &sp : kv.second) {
      sp.setLWR(std::exp(sp.getLikelihood() - max_llh) / total_llh);
    }
  }
}

auto normalize_state_distribution_by_lwr(
    const std::unique_ptr<LagrangeMatrixBase[]> &states,
    size_t states_len) -> std::unique_ptr<LagrangeMatrixBase[]> {
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

  assert(std::isfinite(max_llh));

  double total_llh = 0.0;
  for (size_t i = 1; i < states_len; i++) {
    total_llh += std::exp(normalized_states.get()[i] - max_llh);
  }
  assert(std::isfinite(total_llh));

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
    CONVERT_FLOAT_TO_JSON(tmp["llh"], llh);
    double ratio = lwr_distribution.get()[dist_index];
    CONVERT_FLOAT_TO_JSON(tmp["ratio"], ratio);
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
      anc_json["regions"] = lagrange_convert_dist_to_list(dist, region_names);

      left_json["distribution"] = sp.l_dist;
      left_json["regions"] =
          lagrange_convert_dist_to_list(sp.l_dist, region_names);

      right_json["distribution"] = sp.r_dist;
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
    if (!n.hasResults()) { return; }
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

nlohmann::json produce_json_file(const std::shared_ptr<Tree> &tree,
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

  root_json["node-results"] = make_results_for_nodes(
      tree, config.area_names(), context.stateCount(), config.max_areas());
  return root_json;
}

inline std::string make_csv_header(
    const std::initializer_list<const char *> &fields) {
  return std::accumulate(std::next(fields.begin()),
                         fields.end(),
                         std::string(*fields.begin()),
                         [](const std::string &acc, const std::string &entry)
                             -> std::string { return acc + ", " + entry; })
         + "\n";
}

inline std::string make_csv_row(
    const std::initializer_list<std::string> &entries) {
  return std::accumulate(std::next(entries.begin()),
                         entries.end(),
                         std::string(*entries.begin()),
                         [](const std::string &acc, const std::string &entry)
                             -> std::string { return acc + ", " + entry; })
         + "\n";
}

inline std::ofstream init_csv(
    const std::filesystem::path &filename,
    const std::initializer_list<const char *> &fields) {
  std::ofstream csv_file(filename);
  csv_file << make_csv_header(fields);
  return csv_file;
}

void write_csv_state_file(const std::shared_ptr<Tree> &tree,
                          const ConfigFile &config,
                          const Context &context) {
  LOG_INFO("Writing ancestral states to {}",
           config.statesCSVResultsFilename().string());
  auto fields = {"node", "dist", "llh", "ratio"};
  auto outfile = init_csv(config.statesCSVResultsFilename(), fields);

  auto max_areas = config.max_areas();
  auto total_states = context.stateCount();

  auto cb = [&outfile, max_areas, total_states](const Node &n) {
    if (!n.hasAncestralState()) { return; }
    auto node_label = n.getNodeLabel();
    const auto &states = n.getAncestralState();

    auto lwr_distribution =
        normalize_state_distribution_by_lwr(states, total_states);
    Range dist = 0;
    size_t dist_index = 0;
    Range incl_areas =
        n.getIncludedAreas().hasValue() ? n.getIncludedAreas().get() : 0;
    Range excl_areas =
        n.getExcludedAreas().hasValue() ? n.getExcludedAreas().get() : 0;
    while (true) {
      {
        auto tmp =
            next_dist(dist, max_areas, dist_index, excl_areas, incl_areas);
        dist = tmp.first;
        dist_index = tmp.second;
      }
      if (dist_index >= total_states) { break; }
      auto lwr = lwr_distribution.get()[dist_index];
      if (lwr < 1e-6) { continue; }
      outfile << make_csv_row(
          {node_label,
           std::to_string(dist),
           std::to_string(states.get()[dist_index]),
           std::to_string(lwr_distribution.get()[dist_index])});
    }
  };
  tree->applyPreorderInternalOnly(cb);
}

void write_csv_split_file(const std::shared_ptr<Tree> &tree,
                          const ConfigFile &config) {
  LOG_INFO("Writing ancestral splits to {}",
           config.splitsCSVResultsFilename().string());
  auto fields = {"node", "anc-dist", "left-dist", "right-dist", "llh", "ratio"};
  auto outfile = init_csv(config.splitsCSVResultsFilename(), fields);

  auto cb = [&outfile](const Node &n) {
    if (!n.hasAncestralSplit()) { return; }
    auto node_label = n.getNodeLabel();
    auto splits = n.getAncestralSplit();

    normalize_split_distribution_by_lwr(splits);
    for (const auto &kv : splits) {
      for (const auto &sp : kv.second) {
        if (sp.getLWR() < 1e-6) { continue; }
        outfile << make_csv_row({node_label,
                                 std::to_string(sp.anc_dist),
                                 std::to_string(sp.l_dist),
                                 std::to_string(sp.r_dist),
                                 std::to_string(sp.getLikelihood()),
                                 std::to_string(sp.getLWR())});
      }
    }
  };
  tree->applyPreorderInternalOnly(cb);
}

void write_csv_periods_file(const ConfigFile &config, const Context &context) {
  LOG_INFO("Writing period parameters to {}",
           config.periodsCSVResultsFilename().string());
  auto fields = {"from", "to", "dispersion", "extinction"};
  auto outfile = init_csv(config.periodsCSVResultsFilename(), fields);

  auto p = PeriodSpan(config.periods()).begin();
  double last = 0.0;
  for (const auto &period : context.currentParams()) {
    auto seg = *p;
    outfile << make_csv_row({std::to_string(last),
                             std::to_string(last + seg.duration),
                             std::to_string(period.dispersion_rate),
                             std::to_string(period.extinction_rate)});
    last += seg.duration;
    ++p;
  }
}

void write_csv_distribution_file(const ConfigFile &config) {
  LOG_INFO("Writing distribution labels to {}",
           config.distributionsCSVResultsFilename().string());
  auto fields = {"dist", "list"};
  auto outfile = init_csv(config.distributionsCSVResultsFilename(), fields);

  size_t total_states = 1ul << config.region_count();

  for (Range dist = 0; dist < total_states; ++dist) {
    outfile << dist << ",";
    auto region_list = lagrange_convert_dist_to_list(dist, config.area_names());

    outfile << "\"[";
    if (!region_list.empty()) {
      auto region_list_string = std::accumulate(
          std::next(region_list.begin()),
          region_list.end(),
          "'" + *region_list.begin() + "'",
          [](const std::string &acc, const std::string &entry) -> std::string {
            return acc + ", '" + entry + "'";
          });
      outfile << region_list_string;
    }
    outfile << "]\"\n";
  }
}

void write_csv_node_info_file(const std::shared_ptr<Tree> &tree,
                              const ConfigFile &config) {
  LOG_INFO("Writing node info to {}",
           config.distributionsCSVResultsFilename().string());
  auto fields = {"node", "left-child", "right-child"};
  auto outfile = init_csv(config.nodeInfoCSVResultsFilename(), fields);

  auto cb = [&outfile](const Node &n) {
    outfile << make_csv_row({n.getNodeLabel(),
                             n.getChild(0)->getNodeLabel(),
                             n.getChild(1)->getNodeLabel()});
  };

  tree->applyPreorderInternalOnly(cb);
}

void write_csv_result_files(const std::shared_ptr<Tree> &tree,
                            const ConfigFile &config,
                            const Context &context) {
  write_csv_distribution_file(config);
  write_csv_periods_file(config, context);
  write_csv_node_info_file(tree, config);
  if (config.computeStates()) { write_csv_state_file(tree, config, context); }
  if (config.computeSplits()) { write_csv_split_file(tree, config); }
}

void write_result_files(const std::shared_ptr<Tree> &tree,
                        const ConfigFile &config,
                        const Context &context) {
  if (config.output_file_type() == OutputType::JSON) {
    auto root_json = produce_json_file(tree, config, context);
    write_json_file(config, root_json);
    return;
  }
  if (config.output_file_type() == OutputType::CSV) {
    write_csv_result_files(tree, config, context);
  }
}

inline auto make_annotated_node_newick_lambda() {
  return [](const Node &n) -> std::string {
    std::ostringstream oss;
    oss << n.getNodeLabel();
    if (n.isInternal()) {
      oss << "[&&NHX:left-child=" + n.getChild(0)->getNodeLabel()
                 + ":right-child=" + n.getChild(1)->getNodeLabel() + "]";
    }
    oss << ":" + std::to_string(n.getBL());
    return oss.str();
  };
}

inline auto make_annotated_splits_newick_lambda(
    const std::vector<std::string> &region_names) {
  return [&](const Node &n) -> std::string {
    std::ostringstream oss;
    if (n.isTip()) { oss << n.getNodeLabel(); }
    if (n.isInternal() && n.hasAncestralSplit()) {
      auto anc_split = n.getTopAncestralSplit();
      oss << lagrange_convert_dist_string(anc_split.l_dist, region_names);
      oss << "|";
      oss << lagrange_convert_dist_string(anc_split.r_dist, region_names);
    }
    oss << ":" + std::to_string(n.getBL());
    return oss.str();
  };
}

inline auto make_annotated_states_newick_lambda(
    const std::vector<std::string> &region_names) {
  return [&](const Node &n) -> std::string {
    std::ostringstream oss;
    if (n.isTip()) { oss << n.getNodeLabel(); }
    if (n.isInternal() && n.hasAncestralState()) {
      oss << lagrange_convert_dist_string(n.getTopAncestralState(),
                                          region_names);
    }
    oss << ":" + std::to_string(n.getBL());
    return oss.str();
  };
}

void write_node_tree(const std::shared_ptr<Tree> &tree,
                     const ConfigFile &config) {
  auto node_tree_filename = config.nodeTreeFilename();
  LOG(INFO, "Writing node annotated tree to {}", node_tree_filename.string());

  std::ofstream node_tree(node_tree_filename);
  node_tree << tree->getNewickLambda(make_annotated_node_newick_lambda())
            << std::endl;
}

void write_scaled_tree(const std::shared_ptr<Tree> &tree,
                       const ConfigFile &config) {
  auto scaled_tree_filename = config.scaledTreeFilename();
  LOG(INFO, "Writing scaled tree to {}", scaled_tree_filename.string());

  std::ofstream anal_tree(scaled_tree_filename);
  anal_tree << tree->getNewickLambda(make_annotated_node_newick_lambda())
            << std::endl;
}

void write_states_tree(const std::shared_ptr<Tree> &tree,
                       const ConfigFile &config) {
  auto states_tree_filename = config.statesTreeFilename();
  LOG(INFO, "Writing states tree to {}", states_tree_filename.string());

  std::ofstream states_tree(states_tree_filename);
  states_tree << tree->getNewickLambda(
      make_annotated_states_newick_lambda(config.area_names()))
              << std::endl;
}

void write_splits_tree(const std::shared_ptr<Tree> &tree,
                       const ConfigFile &config) {
  auto splits_tree_filename = config.splitsTreeFilename();
  LOG(INFO, "Writing splits tree to {}", splits_tree_filename.string());

  std::ofstream splits_tree(splits_tree_filename);
  splits_tree << tree->getNewickLambda(
      make_annotated_splits_newick_lambda(config.area_names()))
              << std::endl;
}

auto init_json(const std::shared_ptr<const Tree> &tree,
               const ConfigFile &config) -> nlohmann::json {
  nlohmann::json root_json;
  nlohmann::json attributes_json;
  attributes_json["periods"] =
      !config.periods().empty() ? static_cast<int>(config.periods().size()) : 1;
  attributes_json["regions"] = config.region_count();
  attributes_json["taxa"] = tree->getTipCount();
  attributes_json["nodes-tree"] =
      tree->getNewickLambda(make_annotated_node_newick_lambda());
  attributes_json["max-areas"] = config.max_areas();
  attributes_json["state-count"] = lagrange_compute_restricted_state_count(
      config.region_count(), config.max_areas());
  root_json["attributes"] = attributes_json;

  return root_json;
}

void write_json_file(const ConfigFile &config,
                     const nlohmann ::json &root_json) {
  auto json_filename = config.jsonResultsFilename();

  LOG(INFO, "Writing results to {}", json_filename.string());

  std::ofstream outfile(config.jsonResultsFilename());
  outfile << root_json.dump();
}

void print_run_header(const ConfigFile &config) {
  LOG(INFO, "Lagrange-ng version: {}", GIT_REV_STRING)
  LOG(INFO, "Using tree: {}", config.tree_filename().c_str());
  LOG(INFO, "Using data: {}", config.data_filename().c_str());
}
}  // namespace lagrange
