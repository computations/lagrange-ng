#include "CSV.hpp"

#include <numeric>

auto make_csv_row(const std::initializer_list<std::string_view> &fields)
    -> std::string {
  return std::accumulate(std::next(fields.begin()),
                         fields.end(),
                         std::string(*fields.begin()),
                         [](const std::string &acc,
                            const std::string_view &entry) -> std::string {
                           return std::format("{}, {}", acc, entry);
                         })
         + "\n";
}

std::vector<std::string> read_row(std::string_view row) {
  return row | std::views::chunk_by([](auto, auto c) -> bool {
           if (c == ',') { return false; }
           return true;
         })
         | std::views::transform([](const auto &a) -> auto {
             auto beg_itr = a.begin();
             auto end_itr = std::prev(a.end());

             while (
                 beg_itr != a.end()
                 && (*beg_itr == ','
                     || std::isspace(static_cast<unsigned char>(*beg_itr)))) {
               beg_itr = std::next(beg_itr);
             }

             while (end_itr != a.begin()
                    && std::isspace(static_cast<unsigned char>(*end_itr))) {
               end_itr = std::prev(end_itr);
             }
             end_itr = std::next(end_itr);

             if (beg_itr >= end_itr) {
               throw CSVValueError{"Expected a value"};
             }
             return std::ranges::subrange{beg_itr, end_itr};
           })
         | std::ranges::to<std::vector<std::string>>();
}

template <>
uint64_t CSVRow::get(const std::string_view &key) const {
  return std::stoul(get_data(key));
}

template <>
double CSVRow::get(const std::string_view &key) const {
  return std::stod(get_data(key));
}

template <>
long CSVRow::get(const std::string_view &key) const {
  return std::stol(get_data(key));
}

template <>
bool CSVRow::get(const std::string_view &key) const {
  return std::stoi(get_data(key));
}
