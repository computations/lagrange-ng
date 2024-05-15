#include "ConfigFile.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include "Fossil.h"

enum class config_lexeme_type_t { VALUE, EQUALS_SIGN, END };

void set_mrcas_for_fossils(ConfigFile &config) {
  for (auto &f : config.fossils) { f.clade = config.mrcas[f.mrca_name]; }
}

void check_prefix(ConfigFile &config) {
  if (config.prefix.empty()) { config.prefix = config.treefile; }
  if (config.prefix.filename().string()[0] == '.') {
    std::cout << "Warning, the current prefix starts with a dot, results will "
                 "be hidden on linux systems"
              << std::endl;
  }
}

void finalize(ConfigFile &config) {
  set_mrcas_for_fossils(config);
  check_prefix(config);
}

class config_lexer_t {
 public:
  explicit config_lexer_t(std::string input)
      : _input{std::move(input)}, _current_index{0} {};

  auto consume() -> config_lexeme_type_t {
    auto token = peak();
    if (token == config_lexeme_type_t::VALUE) {
      std::stringstream builder;
      while (char tmp = _input[_current_index]) {
        if (is_punct(tmp) || std::isspace(tmp)) { break; }
        builder << tmp;
        _current_index++;
      }

      _value = builder.str();
      skip_whitespace();
      return token;
    }
    _current_index++;
    skip_whitespace();
    return token;
  }

  auto peak() -> config_lexeme_type_t {
    size_t tmp_index = _current_index;
    char current_char = _input[tmp_index++];

    if (_current_index == _input.size()) { return config_lexeme_type_t::END; }
    if (current_char == '=') { return config_lexeme_type_t::EQUALS_SIGN; }
    return config_lexeme_type_t::VALUE;
  }

  auto consume_value_as_string() -> std::string {
    std::string tmp;
    std::swap(tmp, _value);
    return tmp;
  }

  auto consume_value_as_float() -> double {
    auto f_str = consume_value_as_string();
    size_t pos = 0;
    double val = std::stod(f_str, &pos);
    if (pos != f_str.size()) {
      throw std::runtime_error{std::string("Float conversion failed around") +
                               describe_position()};
    }
    return val;
  }

  auto consume_value_as_size_t() -> size_t {
    auto f_str = consume_value_as_string();
    size_t pos = 0;
    size_t val = std::stoull(f_str, &pos);
    if (pos != f_str.size()) {
      throw std::runtime_error{std::string("Float conversion failed around") +
                               describe_position()};
    }
    return val;
  }

  auto describe_position() const -> std::string {
    std::stringstream builder;
    builder << "position " << _current_index;
    return builder.str();
  }

  void expect(config_lexeme_type_t token_type) {
    auto ret = consume_token_pos();
    if (ret.first != token_type) {
      throw std::runtime_error{
          std::string("Got the wrong token type at position ") +
          std::to_string(ret.second + 1) + " was expecting " +
          describe_token(token_type)};
    }
  }

  void consume_until(config_lexeme_type_t token_type) {
    while (token_type != consume()) {}
  }

  auto at_end() -> bool { return _input.size() == _current_index; }

 private:
  bool is_punct(char c) { return c == '='; }

  auto consume_token_pos() -> std::pair<config_lexeme_type_t, size_t> {
    auto start_index = _current_index;
    auto token = consume();
    return {token, start_index};
  }

  void skip_whitespace() {
    while (_current_index < _input.size()) {
      char c = _input[_current_index];
      if (std::isspace(c) == 0) { break; }
      _current_index++;
    }
  }

  static auto describe_token(config_lexeme_type_t token_type) -> std::string {
    switch (token_type) {
      case config_lexeme_type_t::VALUE:
        return {"either a number or a string"};
      case config_lexeme_type_t::EQUALS_SIGN:
        return {"equals sign"};
      case config_lexeme_type_t::END:
        return {"end of line"};
      default:
        return {"unknown token"};
    }
  }

  std::string _input;
  std::string _value;
  size_t _current_index;
};

std::string parse_filename(config_lexer_t &lexer) {
  lexer.expect(config_lexeme_type_t::VALUE);
  auto filename_value = lexer.consume_value_as_string();
  return filename_value;
}

std::vector<std::string> parse_list(config_lexer_t &lexer) {
  std::vector<std::string> values;
  while (lexer.peak() != config_lexeme_type_t::END) {
    lexer.expect(config_lexeme_type_t::VALUE);
    values.push_back(lexer.consume_value_as_string());
  }
  return values;
}

double parse_double(config_lexer_t &lexer) {
  lexer.expect(config_lexeme_type_t::VALUE);

  return lexer.consume_value_as_float();
}

double parse_size_t(config_lexer_t &lexer) {
  lexer.expect(config_lexeme_type_t::VALUE);

  return lexer.consume_value_as_size_t();
}

Fossil parse_node_fossil(config_lexer_t &lexer) {
  lexer.expect(config_lexeme_type_t::VALUE);
  auto fossile_mrca = lexer.consume_value_as_string();

  lexer.expect(config_lexeme_type_t::VALUE);
  auto fossil_area = lagrange_convert_dist_binary_string_to_dist(
      lexer.consume_value_as_string());

  return {.mrca_name = fossile_mrca,
          .clade = {},
          .age = 0.0,
          .area = fossil_area,
          .type = fossil_type::node};
}

Fossil parse_branch_fossil(config_lexer_t &lexer) {
  lexer.expect(config_lexeme_type_t::VALUE);
  auto fossile_mrca = lexer.consume_value_as_string();

  lexer.expect(config_lexeme_type_t::VALUE);
  auto fossil_area = lagrange_convert_dist_binary_string_to_dist(
      lexer.consume_value_as_string());

  lexer.expect(config_lexeme_type_t::VALUE);
  double fossil_age = lexer.consume_value_as_float();

  return {.mrca_name = fossile_mrca,
          .clade = {},
          .age = fossil_age,
          .area = fossil_area,
          .type = fossil_type::branch};
};

Fossil parse_fixed_fossil(config_lexer_t &lexer) {
  lexer.expect(config_lexeme_type_t::VALUE);
  auto fossile_mrca = lexer.consume_value_as_string();

  lexer.expect(config_lexeme_type_t::VALUE);
  auto fossil_area = lagrange_convert_dist_binary_string_to_dist(
      lexer.consume_value_as_string());

  return {.mrca_name = fossile_mrca,
          .clade = {},
          .age = 0,
          .area = fossil_area,
          .type = fossil_type::fixed};
}

Fossil parse_fossil(config_lexer_t &lexer) {
  lexer.expect(config_lexeme_type_t::VALUE);
  auto fossil_type_string = lexer.consume_value_as_string();

  std::transform(fossil_type_string.cbegin(), fossil_type_string.cend(),
                 fossil_type_string.begin(),
                 [](char c) -> char { return std::tolower(c); });

  if (fossil_type_string == "n" || fossil_type_string == "node") {
    return parse_node_fossil(lexer);
  } else if (fossil_type_string == "b" || fossil_type_string == "branch") {
    return parse_branch_fossil(lexer);
  } else if (fossil_type_string == "f" || fossil_type_string == "fixed") {
    return parse_fixed_fossil(lexer);
  }
  throw std::runtime_error{"Unknown fossile type"};
}

ConfigFile parse_config_file(std::istream &instream) {
  ConfigFile config;
  std::string line;
  size_t line_number = 1;
  while (getline(instream, line)) {
    config_lexer_t lexer(line);
    try {
      if (lexer.peak() == config_lexeme_type_t::END) { continue; }

      lexer.expect(config_lexeme_type_t::VALUE);
      auto config_value = lexer.consume_value_as_string();

      if (config_value == "treefile") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        config.treefile = parse_filename(lexer);
      } else if (config_value == "datafile") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        config.datafile = parse_filename(lexer);
      } else if (config_value == "ratematrix") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        config.ratematrixfile = parse_filename(lexer);
      } else if (config_value == "areanames") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        config.areaNames = parse_list(lexer);
      } else if (config_value == "prefix") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        config.prefix = parse_filename(lexer);
      } else if (config_value == "periods") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        auto tmp_values = parse_list(lexer);

        std::vector<double> time_points;
        time_points.resize(tmp_values.size());

        std::transform(
            tmp_values.begin(), tmp_values.end(), time_points.begin(),
            [](const std::string &s) -> double { return std::stof(s); });

        config.periods = Periods(time_points);
      } else if (config_value == "mrca") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        lexer.expect(config_lexeme_type_t::VALUE);
        auto mrca_name = lexer.consume_value_as_string();

        auto mrca_entries = parse_list(lexer);
        config.mrcas[mrca_name] = std::shared_ptr<MRCAEntry>(
            new MRCAEntry{.clade = std::move(mrca_entries)});
      } else if (config_value == "ancstate") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        config.ancstates = parse_list(lexer);
      } else if (config_value == "fossil") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);

        config.fossils.push_back(parse_fossil(lexer));
      } else if (config_value == "calctype") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        lexer.expect(config_lexeme_type_t::VALUE);

        auto calc_type_value = lexer.consume_value_as_string();
        if (std::tolower(calc_type_value[0]) == 'm') {
          config.marginal = true;
        } else {
          config.marginal = false;
        }
      } else if (config_value == "report") {
        lexer.expect(config_lexeme_type_t::VALUE);
        auto report_value = lexer.consume_value_as_string();
        if (report_value != "split") { config.splits = false; }
      } else if (config_value == "splits") {
        config.splits = true;
      } else if (config_value == "states") {
        config.states = true;
      } else if (config_value == "dispersal") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
      } else if (config_value == "extinction") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);

        config.extinction = parse_double(lexer);
      } else if (config_value == "lh-epsilon") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);

        config.extinction = parse_double(lexer);
      } else if (config_value == "workers") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);

        config.workers = parse_size_t(lexer);
      } else if (config_value == "threads-per-worker") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);

        config.threads_per_worker = parse_size_t(lexer);
      } else if (config_value == "maxareas") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);

        config.maxareas = parse_size_t(lexer);
      } else if (config_value == "expm-mode") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        lexer.expect(config_lexeme_type_t::VALUE);

        auto expm_type = lexer.consume_value_as_string();

        if (expm_type == "krylov") {
          config.expm_mode = lagrange_expm_computation_mode::KRYLOV;
        } else if (expm_type == "pade") {
          config.expm_mode = lagrange_expm_computation_mode::PADE;
        } else if (expm_type == "adaptive") {
          config.expm_mode = lagrange_expm_computation_mode::ADAPTIVE;
        }
      } else if (config_value == "mode") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        lexer.expect(config_lexeme_type_t::VALUE);

        auto mode_type = lexer.consume_value_as_string();
        if (mode_type == "optimize") {
          config.mode = lagrange_operation_mode::OPTIMIZE;
        } else if (mode_type == "evaluate") {
          config.mode = lagrange_operation_mode::EVALUATE;
        }
      } else if (config_value == "alignment-type") {
        lexer.expect(config_lexeme_type_t::EQUALS_SIGN);
        lexer.expect(config_lexeme_type_t::VALUE);

        auto align_type_str = lexer.consume_value_as_string();
        if (align_type_str == "fasta") {
          config.alignment_file_type = AlignmentFileType::fasta;
        }
        if (align_type_str == "phylip") {
          config.alignment_file_type = AlignmentFileType::phylip;
        }
      } else {
        std::stringstream oss;
        oss << "Option '" << config_value << "' on line " << line_number
            << " was not recognized";
        throw std::runtime_error{oss.str()};
      }

      lexer.expect(config_lexeme_type_t::END);
    } catch (const std::exception &e) {
      std::ostringstream oss;
      oss << "There was a problem parsing line " << line_number
          << " of the config file:\n"
          << e.what();
      throw std::runtime_error{oss.str()};
    }

    line_number++;
  }

  finalize(config);

  return config;
}

std::filesystem::path ConfigFile::get_results_filename() {
  auto results_filename = prefix;
  results_filename += ".results.json";
  return results_filename;
}

std::filesystem::path ConfigFile::get_node_tree_filename() {
  auto node_tree_filename = prefix;
  node_tree_filename += ".nodes.tre";
  return node_tree_filename;
}

std::filesystem::path ConfigFile::get_scaled_tree_filename() {
  auto scaled_tree_filename = prefix;
  scaled_tree_filename += ".scaled.tre";
  return scaled_tree_filename;
}
