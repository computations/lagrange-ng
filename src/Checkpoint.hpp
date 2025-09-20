#pragma once

#include <filesystem>
#include <fstream>
#include <optional>
#include <ranges>

#include "CSV.hpp"
#include "ConfigFileParser.hpp"
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

  Checkpoint(Checkpoint &&) = default;
  Checkpoint &operator=(Checkpoint &&) = default;

  void logCheckpoint(const lagrange::PeriodParameterRange auto &parameters,
                     size_t iter,
                     bool final) {
    if (!_checkpoint_file) { initCheckpoint(parameters); }
    std::vector<std::string> parameter_list;
    parameter_list.push_back(std::format("{}", iter));
    for (auto &p : parameters) {
      parameter_list.push_back(std::format("{}", p.dispersion_rate));
      parameter_list.push_back(std::format("{}", p.extinction_rate));
      if (p.hasAdjustmentMatrix()) {
        parameter_list.push_back(std::format("{}", p.distance_penalty));
      }
    }
    parameter_list.push_back(std::format("{:d}", final));
    write_csv_row(*_checkpoint_file, parameter_list);
  }

  bool existingCheckpoint() const { return std::filesystem::exists(_filename); }

  bool isFinalized() const { return _final; }

  std::vector<double> loadCheckpoint() {
    LOG_INFO("Loading checkpoint {}", _filename.string());
    CSVReader reader{_filename};
    parseFields(reader.header());
    auto rows = reader.read_rows();
    if (rows.size() == 0) { return {}; }

    std::vector<double> loaded_params;

    auto last_row = rows.back();
    for (size_t i = 0; i < _fields.size(); ++i) {
      auto [_, value] = last_row.get<std::string>(i + 1);
      loaded_params.push_back(std::stod(value));
    }

    _final = last_row.get<bool>("final");

    _checkpoint_file = std::ofstream{_filename, std::ios::app};
    return loaded_params;
  }

 private:
  void parseFields(const CSVHeaderType &header) {
    for (auto &h : header) {
      if (h == "iteration" || h == "final") { continue; }
      _fields.push_back(FieldInfo::parse(h));
    }
  }

  void makeFields(const lagrange::PeriodParameterRange auto &parameters) {
    for (auto param : parameters) {
      _fields.push_back({.period_name = param.name,
                         .type = FieldInfo::parameter_type::dispersion});
      _fields.push_back({.period_name = param.name,
                         .type = FieldInfo::parameter_type::extinction});
      if (param.hasAdjustmentMatrix()) {
        _fields.push_back({.period_name = param.name,
                           .type = FieldInfo::parameter_type::distance});
      }
    }
  }

  void initCheckpoint(const lagrange::PeriodParameterRange auto &parameters) {
    makeFields(parameters);

    _checkpoint_file = std::ofstream(_filename);
    std::vector<std::string> header;

    header.push_back("iteration");
    for (auto &f : _fields) { header.push_back(f.to_string()); }
    header.push_back("final");

    write_csv_row(*_checkpoint_file, header);
  }

  struct FieldInfo {
    std::string period_name;
    enum class parameter_type {
      dispersion,
      extinction,
      distance,
    } type;

    std::string to_string() const {
      return std::format("{}:{}", period_name, describeParameter(type));
    }

    static std::string_view describeParameter(parameter_type type) {
      using namespace std::string_view_literals;
      switch (type) {
        case parameter_type::dispersion:
          return "d"sv;
        case parameter_type::extinction:
          return "e"sv;
        case parameter_type::distance:
          return "p"sv;
      }
      throw std::runtime_error{"Failed to describe parameter type"};
    }

    static parameter_type parseType(const std::string_view &key) {
      if (key == "d") { return parameter_type::dispersion; }
      if (key == "e") { return parameter_type::extinction; }
      if (key == "p") { return parameter_type::distance; }
      throw std::runtime_error{"Failed to parse parameter type"};
    }

    static FieldInfo parse(std::string_view header) {
      size_t split_index = header.find_last_of(":");
      auto type = parseType(header.substr(split_index + 1));
      return {.period_name = std::string(header.substr(0, split_index)),
              .type = type};
    }
  };

  std::filesystem::path _filename;
  std::optional<std::ofstream> _checkpoint_file;
  std::vector<FieldInfo> _fields;
  bool _final = false;
};
