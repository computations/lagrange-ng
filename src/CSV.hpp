#ifndef CSV_READER_H
#define CSV_READER_H

#include <filesystem>
#include <fstream>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

using CSVHeaderType = std::vector<std::string>;
using CSVHeaderPointer = std::shared_ptr<CSVHeaderType>;

class CSVRow {
 public:
  CSVRow() = default;

  CSVRow(CSVHeaderPointer header, std::vector<std::string>&& data) :
      _header{header},
      _data{data} {}

  template <typename T>
  T get(const std::string_view& key) {
    return get_data(key);
  }

 private:
  std::string get_data(const std::string_view& key) {
    for (auto [h_key, data] : std::views::zip(*_header, _data)) {
      if (h_key == key) { return data; }
    }
    return {};
  }

  CSVHeaderPointer _header;
  std::vector<std::string> _data;
};

template <>
uint64_t CSVRow::get(const std::string_view& key);

template <>
double CSVRow::get(const std::string_view& key);

std::vector<std::string> read_row(std::string_view row);

class CSVIter {
 public:
  using value_type = CSVRow;
  using difference_type = size_t;
  typedef CSVRow* pointer;
  typedef CSVRow& reference;

  CSVIter(std::shared_ptr<std::ifstream> csv_file, CSVHeaderPointer header) :
      _csv_file{csv_file},
      _header{header} {}

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
    std::string row_buffer;
    std::getline(*_csv_file, row_buffer);

    if (row_buffer.empty()) {
      _csv_file = nullptr;
      return;
    }

    _row = CSVRow{_header, read_row(row_buffer)};
  }

  std::shared_ptr<std::ifstream> _csv_file;
  CSVRow _row;
  CSVHeaderPointer _header;
};

class CSVReader {
  CSVReader(const std::filesystem::path& csv_filename) :
      _csv_file{new std::ifstream{csv_filename}},
      _header{new std::vector<std::string>} {
    std::string header_buffer;
    std::getline(*_csv_file, header_buffer);
    *_header = read_row(header_buffer);
  }

  CSVIter begin() { return {_csv_file, _header}; }

  CSVIter end() { return {nullptr, nullptr}; }

 private:
  std::shared_ptr<std::ifstream> _csv_file;
  CSVHeaderPointer _header;
};

#endif
