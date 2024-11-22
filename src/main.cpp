/*
 * main.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: Stephen A. Smith
 *   Last Edit: 27 Oct 2020
 *      Author: Ben Bettisworth
 */

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <logger.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "Alignment.hpp"
#include "Common.hpp"
#include "ConfigFile.hpp"
#include "Context.hpp"
#include "IO.hpp"
#include "TreeReader.hpp"
#include "Utils.hpp"
#include "WorkerState.hpp"

using namespace lagrange;

static void set_expm_mode(Context &context, const ConfigFile &config) {
  if (config.has_expm_mode()) {
    switch (config.expm_mode()) {
      case LagrangeEXPMComputationMode::ADAPTIVE:
        if (config.region_count() > KRYLOV_RANGE_COUNT_THRESHOLD) {
          LOG(INFO, "Enabling adaptive EXPM mode");
          context.useArnoldi();
        } else {
          LOG(INFO, "Using Pade's method for EXPM computation");
          context.useArnoldi(false, false);
        }
        break;
      case LagrangeEXPMComputationMode::PADE:
        LOG(INFO, "Using Pade's method for EXPM computation");
        context.useArnoldi(false, false);
        break;
      case LagrangeEXPMComputationMode::KRYLOV:
        LOG(INFO,
                "Using Krylov subspace based method for EXPM computation");
        context.useArnoldi(true, false);
        break;
      default:
        throw std::runtime_error{"Unknown Expm Mode"};
    }
  } else {
    if (config.region_count() > KRYLOV_RANGE_COUNT_THRESHOLD) {
      LOG(INFO, "Enabling adaptive EXPM mode");
      context.useArnoldi();
    } else {
      LOG(INFO, "Using Pade's method for EXPM computation");
      context.useArnoldi(false, false);
    }
  }
}

static void setup_tree(const std::shared_ptr<Tree> &tree,
                       const ConfigFile &config) {
  tree->setHeightBottomUp();
  tree->setPeriods(config.periods());
  tree->assignMCRALabels(config.mrcas());
  tree->assignFossils(config.fossils());
}

static void assign_results_to_tree(std::shared_ptr<Tree> &tree,
                                   const ConfigFile &config,
                                   const Context &context) {
  if (config.computeStates()) {
    auto states = context.getStateResults();
    auto cb = [&](Node &n) {
      if (auto result = states.find(n.getId()); result != states.end()) {
        n.assignAncestralState(std::move(result->second));
      }
    };
    tree->applyPreorderInternalOnly(cb);
  }

  if (config.computeSplits()) {
    auto splits = context.getSplitResults();
    auto cb = [&](Node &n) {
      if (auto result = splits.find(n.getId()); result != splits.end()) {
        n.assignAncestralSplit(std::move(result->second));
      }
    };
    tree->applyPreorderInternalOnly(cb);
  }
}

static void handle_tree(std::shared_ptr<Tree> &tree,
                        const Alignment &data,
                        const ConfigFile &config) {
  auto root_json = init_json(tree, config);

  setup_tree(tree, config);

  Context context(tree, config.region_count(), config.max_areas());
  context.registerLHGoal();
  if (config.compute_all_states() || config.compute_all_splits()) {
    context.registerStateLHGoal();
  } else if (config.computeStates()) {
    auto tmp_nodes = config.state_nodes();
    for (const auto &split_node : config.split_nodes()) {
      tmp_nodes.insert(split_node);
    }
    context.registerStateLHGoal(tmp_nodes, config.mrcas());
  }
  if (config.compute_all_splits()) {
    context.registerSplitLHGoal();
  } else if (config.computeSplits()) {
    context.registerSplitLHGoal(config.split_nodes(), config.mrcas());
  }
  context.set_opt_method(config.opt_method());
  context.init();
  context.updateRates({context.getPeriodCount(), config.period_params()});
  context.registerTipClvs(data.data);
  context.set_lh_epsilon(config.lh_epsilon());
  set_expm_mode(context, config);

  std::vector<WorkerState> worker_states;
  worker_states.reserve(config.workers());
  std::vector<std::thread> threads;
  LOG(INFO, "Starting Workers");
  for (size_t i = 0; i < config.workers(); i++) {
    LOG(INFO, "Making Worker #{}", i + 1);
    worker_states.emplace_back();
    worker_states.back().setAssignedThreads(config.threads_per_worker());
    threads.emplace_back(&Context::optimizeAndComputeValues,
                         std::ref(context),
                         std::ref(worker_states[i]),
                         config.computeStates(),
                         config.computeSplits(),
                         std::cref(config.run_mode()));
  }

  LOG(INFO, "Waiting for Workers to finish");
  for (auto &t : threads) { t.join(); }

  assign_results_to_tree(tree, config, context);
  write_result_files(tree, config, context);
  write_node_tree(tree, config);
  if (config.computeStatesStrict()) { write_states_tree(tree, config); }
  if (config.computeSplits()) { write_splits_tree(tree, config); }
}

static ConfigFile read_config_file(const std::string &filename) {
  std::ifstream infile(filename);
  return ConfigFile{infile};
}

static std::vector<std::shared_ptr<Tree>> read_tree_file_line_by_line(
    const std::filesystem::path &filename) {
  std::vector<std::shared_ptr<Tree>> ret;

  if (!std::filesystem::exists(filename)) {
    throw std::runtime_error{"Failed to find the tree file "
                             + filename.string()};
  }

  std::ifstream ifs(filename);
  std::string temp;
  size_t count = 1;
  while (getline(ifs, temp)) {
    if (temp.size() > 1) {
      auto intree = TreeReader::readTree(temp);
      LOG(INFO, "Tree number {} has {} leaves", count, intree->getTipCount());
      ret.push_back(intree);
      count++;
    }
  }
  return ret;
}

bool check_alignment_against_trees(
    const Alignment &align, const std::vector<std::shared_ptr<Tree>> &trees) {
  bool ok = true;
  for (const auto &t : trees) { ok &= t->checkAlignmentConsistency(align); }
  return ok;
}

auto main(int argc, char *argv[]) -> int {
#if MKL_ENABLED
  mkl_set_num_threads(1);
#else
  openblas_set_num_threads(1);
#endif

  logger::get_log_states().add_stream(
      stdout, INFO | IMPORTANT | PROGRESS | WARNING | ERROR);

  auto start_time = std::chrono::high_resolution_clock::now();
  if (argc != 2) {
    LOG(ERROR, "Lagrange-ng needs a config file");
    LOG(ERROR, "Usage: lagrange configfile");
    return 1;
  } else {
    std::string config_filename(argv[1]);
    auto config = read_config_file(config_filename);
    print_run_header(config);

    LOG(INFO, "Reading tree...");
    std::vector<std::shared_ptr<Tree>> intrees =
        read_tree_file_line_by_line(config.tree_filename());

    LOG(INFO, "Reading data...");
    Alignment data =
        read_alignment(config.data_filename(), config.alignment_file_type());

    LOG(INFO, "Checking data...");
    check_alignment_against_trees(data, intrees);

    if (data.region_count == 0) {
      throw std::runtime_error{"Region count cannot be zero"};
    }

    config.region_count(data.region_count);
    if (!config.has_max_areas()) { config.max_areas(config.region_count()); }

    config.finalize_periods();

    LOG(INFO, "Running analysis...");
    for (auto &intree : intrees) { handle_tree(intree, data, config); }
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end_time - start_time;
  LOG(INFO, "Analysis took {:.2}s", duration.count());
  return 0;
}
