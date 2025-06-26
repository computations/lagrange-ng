#include "CSV.hpp"

std::vector<std::string> read_row(std::string_view row) {
  return row | std::views::chunk_by([](auto, auto c) -> bool {
           if (c == ',') { return false; }
           return true;
         })
         | std::ranges::to<std::vector<std::string>>();
}

template <>
uint64_t CSVRow::get(const std::string_view& key) {
  return std::stoul(get_data(key));
}

template <>
double CSVRow::get(const std::string_view& key) {
  return std::stod(get_data(key));
}
