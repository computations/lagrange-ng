#pragma once

#include <filesystem>
#include <fstream>
#include <optional>
#include <ranges>

#include "CSV.hpp"
#include "Periods.hpp"

/**
 * I am pretty sure that we just need to store the current paramters, and then
 * load them back when we restart the program. When doing BFGS, we will loose
 * the estimated hessian, which will be difficult to compute, but this will
 * still be faster than not starting from a good set of parameters.
 *
 * Additionally, if we mark the program as being done, we can just evaluate the
 * model at the given parameter set and output results.
 *
 * I am going to serialize the checkpoint as a append only CSV log.
 */
class Checkpoint {
 public:
  Checkpoint(const std::filesystem::path &checkpoint_filename) :
      _filename{checkpoint_filename} {}

  void write_checkpoint(const lagrange::PeriodParameterRange auto &parameters) {
    if (!_checkpoint_file) { init_checkpoint(parameters); }
    std::vector<double> parameter_list;
    for (auto &p : parameters) {
      parameter_list.push_back(p.dispersion_rate);
      parameter_list.push_back(p.extinction_rate);
      if (p.hasAdjustmentMatrix()) {
        parameter_list.push_back(p.distance_penalty);
      }
    }
    write_csv_row(*_checkpoint_file, parameter_list);
  }

 private:
  void init_checkpoint(const lagrange::PeriodParameterRange auto &parameters) {
    _checkpoint_file = std::ofstream(_filename);
    std::vector<std::string> header;
    for (const auto &[i, param] : std::views::enumerate(parameters)) {
      header.push_back(std::format("d_{}", i));
      header.push_back(std::format("e_{}", i));
      if (param.hasAdjustmentMatrix()) {
        header.push_back(std::format("p_{}", i));
      }
    }

    write_csv_row(*_checkpoint_file, header);
  }

  std::filesystem::path _filename;
  std::optional<std::ofstream> _checkpoint_file;
};
