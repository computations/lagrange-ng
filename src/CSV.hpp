#ifndef CSV_READER_H
#define CSV_READER_H

#include <concepts>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

// #ifndef __cpp_lib_ranges_join_with
// #endif

template <std::ranges::range R>
auto make_csv_row(const R& entries) -> auto
  requires std::convertible_to<std::ranges::range_value_t<R>, std::string>
           || std::convertible_to<std::ranges::range_value_t<R>,
                                  std::string_view>
{
#ifdef __cpp_lib_ranges_join_with
  return entries | std::ranges::views::join_with(',');
#else
  using namespace std::string_literals;
  return std::views::zip(entries, std::views::repeat(","s))
         | std::views::transform([](auto t) -> std::vector<std::string> {
             return {std::string(std::get<0>(t)), std::get<1>(t)};
           })
         | std::views::join | std::views::take(2 * entries.size() - 1);

#endif
}

template <std::ranges::range R>
auto make_csv_row(const R& entries) -> auto
  requires std::floating_point<std::ranges::range_value_t<R>>
{
#ifdef __cpp_lib_ranges_join_with
  return entries | std::views::transform([](const auto& val) -> std::string {
           return std::to_string(val);
         })
         | std::ranges::views::join_with(',');
#else
  using namespace std::string_literals;
  return std::views::zip(entries, std::views::repeat(","s))
         | std::views::transform([](auto t) -> std::vector<std::string> {
             return {std::to_string(std::get<0>(t)), std::get<1>(t)};
           })
         | std::views::join | std::views::take(2 * entries.size() - 1);
#endif
}

template <std::ranges::range R>
auto write_csv_row(std::ostream& os, const R& entries) {
  for (auto e : make_csv_row(entries)) { os << e; }
  os << "\n";
  os.flush();
}

using CSVHeaderType = std::vector<std::string>;
using CSVHeaderPointer = std::shared_ptr<CSVHeaderType>;

class CSVValueError : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

class CSVRow {
 public:
  CSVRow() = default;

  CSVRow(CSVHeaderPointer header, std::vector<std::string>&& data) :
      _header{header},
      _data{data} {}

  template <typename T>
  T get(const std::string_view& key) const {
    return get_data(key);
  }

  template <typename T>
  std::pair<std::string_view, T> get(size_t index) const {
    std::string_view header_val = (*_header)[index];
    return {header_val, get_data(header_val)};
  }

 private:
  std::string get_data(const std::string_view& key) const {
#ifdef __cpp_lib_ranges_zip
    for (auto [h_key, data] : std::views::zip(*_header, _data)) {
      if (h_key == key) { return data; }
    }
#else
    auto header = _header->begin();
    auto data = _data.begin();
    for (; header != _header->end() && data != _data.end(); ++header, ++data) {
      if (*header == key) { return *data; }
    }
#endif
    return {};
  }

  CSVHeaderPointer _header;
  std::vector<std::string> _data;
};

template <>
uint64_t CSVRow::get(const std::string_view& key) const;

template <>
double CSVRow::get(const std::string_view& key) const;

template <>
long CSVRow::get(const std::string_view& key) const;

template <>
bool CSVRow::get(const std::string_view& key) const;

std::vector<std::string> read_row(std::string_view row);

class CSVIter {
 public:
  using value_type = CSVRow;
  using difference_type = size_t;
  typedef CSVRow* pointer;
  typedef CSVRow& reference;

  CSVIter(std::shared_ptr<std::istream> csv_file, CSVHeaderPointer header) :
      _csv_file{csv_file},
      _header{header} {
    get_row();
  }

  CSVIter() : _csv_file{nullptr}, _header{nullptr} {}

  const CSVRow& operator*() const { return _row; }

  const CSVRow* operator->() const { return &_row; }

  CSVIter& operator++() {
    if (_csv_file == nullptr) { return *this; }

    std::string row_buffer;
    std::getline(*_csv_file, row_buffer);

    if (row_buffer.empty()) {
      _csv_file = nullptr;
      return *this;
    }

    _row = CSVRow{_header, read_row(row_buffer)};
    return *this;
  }

  CSVIter operator++(int) {
    auto tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator==(const CSVIter& other) const {
    return (this == &other)
           || ((_csv_file == nullptr) && (other._csv_file == nullptr));
  }

  bool operator!=(const CSVIter& other) const { return !(*this == other); }

 private:
  void get_row() {
    if (_csv_file == nullptr) { return; }
    std::string row_buffer;
    std::getline(*_csv_file, row_buffer);

    if (row_buffer.empty()) {
      _csv_file = nullptr;
      return;
    }

    _row = CSVRow{_header, read_row(row_buffer)};
  }

  std::shared_ptr<std::istream> _csv_file;
  CSVRow _row;
  CSVHeaderPointer _header;
};

class CSVReader {
 public:
  CSVReader(const std::filesystem::path& csv_filename) :
      _csv_file{new std::ifstream{csv_filename}},
      _header{new CSVHeaderType} {
    read_header();
  }

  CSVReader(std::istringstream&& istream) :
      _csv_file{new std::istringstream(std::move(istream))},
      _header{new CSVHeaderType} {
    read_header();
  }

  CSVIter begin() { return {_csv_file, _header}; }

  CSVIter end() { return {nullptr, nullptr}; }

  CSVHeaderType const& header() const { return *_header; }

  std::vector<CSVRow> read_rows() {
    std::vector<CSVRow> rows;
    for (auto row : *this) { rows.push_back(row); }
    return rows;
  }

 private:
  void read_header() {
    std::string header_buffer;
    std::getline(*_csv_file, header_buffer);
    *_header = read_row(header_buffer);
  }

  std::shared_ptr<std::istream> _csv_file;
  CSVHeaderPointer _header;
};

#endif
