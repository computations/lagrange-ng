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
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "Alignment.h"
#include "Common.h"
#include "ConfigFile.h"
#include "Context.h"
#include "IO.hpp"
#include "TreeReader.h"
#include "Utils.h"
#include "WorkerState.h"

static void set_expm_mode(Context &context, const ConfigFile &config) {
  auto &expm_mode = config.expm_mode;
  if (expm_mode.has_value()) {
    switch (expm_mode.get()) {
      case lagrange_expm_computation_mode::ADAPTIVE:
        if (config.region_count > KRYLOV_RANGE_COUNT_THRESHOLD) {
          std::cout << "Enabling adaptive expm computation" << std::endl;
          context.useArnoldi();
        } else {
          std::cout << "Using Pade's method expm computation" << std::endl;
          context.useArnoldi(false, false);
        }
        break;
      case lagrange_expm_computation_mode::PADE:
        std::cout << "Using Pade's method expm computation" << std::endl;
        context.useArnoldi(false, false);
        break;
      case lagrange_expm_computation_mode::KRYLOV:
        std::cout << "Using Krylov subspace based method for expm computation"
                  << std::endl;
        context.useArnoldi(true, false);
        break;
      default:
        throw std::runtime_error{"Unknown Expm Mode"};
    }
  } else {
    if (config.region_count > KRYLOV_RANGE_COUNT_THRESHOLD) {
      std::cout << "Enabling adaptive expm computation" << std::endl;
      context.useArnoldi();
    } else {
      std::cout << "Using Pade's method expm computation" << std::endl;
    }
  }
}

static void setup_tree(const std::shared_ptr<Tree> &tree,
                       const ConfigFile &config) {
  tree->setHeightBottomUp();
  tree->setPeriods(config.periods);
  tree->assignFossils(config.fossils);
}

static void handle_tree(const std::shared_ptr<Tree> &tree,
                        const Alignment &data, const ConfigFile &config) {
  auto root_json = init_json(tree, config);

  setup_tree(tree, config);

  Context context(tree, config.region_count, config.maxareas);
  context.registerLHGoal();
  if (config.states) { context.registerStateLHGoal(); }
  if (config.splits) { context.registerSplitLHGoal(); }
  context.init();
  context.updateRates(
      {context.getPeriodCount(), {config.dispersal, config.extinction}});
  context.registerTipClvs(data.data);
  context.set_lh_epsilon(config.lh_epsilon);
  set_expm_mode(context, config);

  std::vector<WorkerState> worker_states;
  worker_states.reserve(config.workers.get());
  std::vector<std::thread> threads;
  std::cout << "Starting Workers" << std::endl;
  for (size_t i = 0; i < config.workers.get(); i++) {
    std::cout << "Making Worker #" << i + 1 << std::endl;
    worker_states.emplace_back();
    worker_states.back().set_assigned_threads(config.threads_per_worker.get());
    threads.emplace_back(&Context::optimizeAndComputeValues, std::ref(context),
                         std::ref(worker_states[i]), config.states,
                         config.splits, true, std::cref(config.mode.get()));
  }

  std::cout << "Waiting for workers to finish" << std::endl;
  for (auto &t : threads) { t.join(); }

  write_result_file(tree, config, context);
  write_node_tree(tree, config);
  write_scaled_tree(tree, config);
}

static void setThreads(ConfigFile &config) {
  if (!config.workers.has_value()) { config.workers = 1; }
  if (!config.threads_per_worker.has_value()) { config.threads_per_worker = 1; }
}

static ConfigFile read_config_file(const std::string &filename) {
  std::ifstream infile(filename);
  return parse_config_file(infile);
}

static std::vector<std::shared_ptr<Tree>> read_tree_file_line_by_line(
    const std::string &filename) {
  std::vector<std::shared_ptr<Tree>> ret;
  std::ifstream ifs(filename);
  std::string temp;
  int count = 1;
  while (getline(ifs, temp)) {
    if (temp.size() > 1) {
      auto intree = TreeReader::readTree(temp);
      std::cout << "Tree " << count << " has " << intree->getExternalNodeCount()
                << " leaves." << std::endl;
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

bool validate_and_make_prefix(const std::filesystem::path &prefix) {
  bool ok = true;
  try {
    if (!prefix.parent_path().empty()) {
      std::filesystem::create_directories(prefix.parent_path());
    }
  } catch (const std::filesystem::filesystem_error &err) {
    std::cerr << err.what() << std::endl;
    ok = false;
  }

  return ok;
}

auto main(int argc, char *argv[]) -> int {
#if MKL_ENABLED
  mkl_set_num_threads(1);
#else
  openblas_set_num_threads(1);
#endif
  auto start_time = std::chrono::high_resolution_clock::now();
  if (argc != 2) {
    std::cout << "you need more arguments." << std::endl;
    std::cout << "usage: lagrange configfile" << std::endl;
    return 0;
  } else {
    std::string config_filename(argv[1]);
    auto config = read_config_file(config_filename);
    validate_and_make_prefix(config.prefix);
    setThreads(config);

    std::cout << "reading tree..." << std::endl;
    std::vector<std::shared_ptr<Tree>> intrees =
        read_tree_file_line_by_line(config.treefile);

    std::cout << "reading data..." << std::endl;
    Alignment data =
        read_alignment(config.datafile, config.alignment_file_type);

    std::cout << "checking data..." << std::endl;
    check_alignment_against_trees(data, intrees);

    if (data.region_count == 0) {
      throw std::runtime_error{"Region count cannot be zero"};
    }

    config.region_count = data.region_count;
    if (config.maxareas == 0) { config.maxareas = config.region_count; }

    std::cout << "running analysis..." << std::endl;
    for (auto &intree : intrees) { handle_tree(intree, data, config); }
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end_time - start_time;
  std::cout << "Analysis took: " << duration.count() << "s" << std::endl;
  return 0;
}
