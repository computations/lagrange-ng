#include "CSV.hpp"

std::vector<std::string> read_row(std::string_view row) {
  return row | std::views::chunk_by([](auto, auto c) -> bool {
           if (c == ',') { return false; }
           return true;
         })
         | std::views::transform([](const auto& a) -> auto {
             auto beg_itr = a.begin();
             auto end_itr = std::prev(a.end());
             while (*beg_itr == ',' || std::isspace(*beg_itr)) {
               beg_itr = std::next(beg_itr);
             }
             if (*end_itr == ' ') { end_itr = std::prev(end_itr); }
             end_itr = std::next(end_itr);
             if (end_itr - beg_itr < 0) {
               throw CSVValueError{"Expected a value"};
             }
             return std::ranges::subrange{beg_itr, end_itr};
           })
         | std::ranges::to<std::vector<std::string>>();
}

template <>
uint64_t CSVRow::get(const std::string_view& key) const {
  return std::stoul(get_data(key));
}

template <>
double CSVRow::get(const std::string_view& key) const {
  return std::stod(get_data(key));
}

template <>
long CSVRow::get(const std::string_view& key) const {
  return std::stol(get_data(key));
}
