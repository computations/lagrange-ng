#pragma once
#include <memory>
#include <nlohmann/thirdparty/hedley/hedley.hpp>

#include "AncSplit.hpp"
#include "Common.hpp"
#include "ConfigFile.hpp"
#include "Context.hpp"
#include "nlohmann/json.hpp"

namespace lagrange {


auto normalize_split_distribution_by_lwr(SplitReturn &splits) -> void;

auto normalize_state_distribution_by_lwr(
    const std::unique_ptr<LagrangeMatrixBase[]> &states, size_t states_len)
    -> std::unique_ptr<LagrangeMatrixBase[]>;

auto make_state_results_for_node(
    const std::unique_ptr<LagrangeMatrixBase[]> &state_distribution,
    const std::vector<std::string> &region_names,
    size_t states_len,
    size_t max_areas) -> nlohmann::json;

auto make_split_results_for_node(SplitReturn &splits,
                                 const std::vector<std::string> &region_names,
                                 size_t states_len,
                                 size_t max_areas) -> nlohmann::json;

auto make_results_for_node(
    const std::vector<std::unique_ptr<LagrangeMatrixBase[]>> &states,
    SplitReturnList &splits,
    const std::vector<size_t> &state_id_map,
    const std::vector<std::string> &region_names,
    size_t states_len,
    size_t max_areas) -> nlohmann::json;

void write_result_files(const std::shared_ptr<Tree> &tree,
                       const ConfigFile &config,
                       const Context &context);

void write_node_tree(const std::shared_ptr<Tree> &tree,
                     const ConfigFile &config);

void write_scaled_tree(const std::shared_ptr<Tree> &tree,
                       const ConfigFile &config);

auto init_json(const std::shared_ptr<const Tree> &tree,
               const ConfigFile &config) -> nlohmann::json;

void write_json_file(const ConfigFile &config, const nlohmann::json &root_json);

nlohmann::json produce_json_file(const std::shared_ptr<Tree> &tree,
                                 const ConfigFile &config,
                                 const Context &context);

}  // namespace lagrange
