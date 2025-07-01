#include "AdjustmentMatrix.hpp"

#include <algorithm>
#include <functional>
#include <ranges>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "CSV.hpp"

namespace lagrange {

using namespace std::string_view_literals;

static bool check_duplicates(std::span<const AdjustmentArc> arcs) {
  auto duplicates = arcs | std::views::adjacent<2>
                    | std::views::transform([](const auto& t) -> bool {
                        auto& [a, b] = t;
                        auto [a_from, a_to, a_dist] = a;
                        auto [b_from, b_to, b_dist] = b;
                        return a_from == b_from && a_to == b_to;
                      });
  return std::ranges::none_of(duplicates, std::identity{});
}

/* assume that the matrix is sorted */
/* need to check that:
 * - No duplicates
 * - All arcs are found (both symmetric and non)
 */
static bool is_valid(std::span<const AdjustmentArc> arcs) {
  auto runs = arcs | std::views::chunk_by([](const auto& a, const auto& b) {
                auto [a_from, a_to, a_dist] = a;
                auto [b_from, b_to, b_dist] = b;
                return a_from == b_from;
              })
              | std::views::adjacent<2>;

  auto decreasing = std::ranges::all_of(runs, [](const auto t) -> bool {
    auto& [a, b] = t;
    return a.size() > b.size();
  });

  auto equal = std::ranges::all_of(runs, [](const auto t) -> bool {
    auto& [a, b] = t;
    return a.size() == b.size();
  });

  return check_duplicates(arcs) && (decreasing != equal);
}

void AdjustmentMatrix::read_arcs(CSVReader& reader,
                                 const std::vector<std::string>& area_names) {
  _region_count = area_names.size();
  auto reverse_area_name_map =
      area_names | std::views::enumerate
      | std::views::transform([](const auto& a) -> auto {
          return std::make_pair(std::get<1>(a),
                                static_cast<size_t>(std::get<0>(a)));
        })
      | std::ranges::to<std::unordered_map>();

  constexpr auto from_key = "from"sv;
  constexpr auto to_key = "to"sv;
  constexpr auto dist_key = "dist"sv;

  for (auto& row : reader) {
    auto from = row.get<std::string>(from_key);
    auto to = row.get<std::string>(to_key);
    _arcs.push_back({.from = reverse_area_name_map.at(from),
                     .to = reverse_area_name_map.at(to),
                     .dist = row.get<double>(dist_key)});
  }

  std::sort(
      _arcs.begin(), _arcs.end(), [](const auto& a, const auto& b) -> bool {
        return a.from < b.from || (a.from == b.from && a.to < b.to);
      });
}

AdjustmentMatrix::AdjustmentMatrix(const std::filesystem::path& filename,
                                   const std::vector<std::string>& area_names) {
  CSVReader reader(filename);
  read_arcs(reader, area_names);

  _type = determine_matrix_symmetry(_arcs, area_names.size());
  if (_type == AdjustmentMatrixType::invalid) {
    throw std::runtime_error{"Adjustment matrix is invalid"};
  }
}

/*
 * Tasks:
 *  - Load the CSV file
 *  - For each row, make an "arc", a tuple of (from, to, distance)
 *  - After loading the csv, check if it is symmetric or not
 *  - setup the _real_ matrix based on that determination
 */
AdjustmentMatrix::AdjustmentMatrix(std::istringstream&& matstream,
                                   const std::vector<std::string>& area_names) {
  CSVReader reader(std::move(matstream));
  read_arcs(reader, area_names);

  _type = determine_matrix_symmetry(_arcs, _region_count);
  if (_type == AdjustmentMatrixType::invalid) {
    throw std::runtime_error{"Adjustment matrix is invalid"};
  }
}

std::shared_ptr<double[]> AdjustmentMatrix::to_matrix() const {
  auto ret = std::shared_ptr<double[]>{new double[compute_size()]};

  for (auto arc : _arcs) {
    ret[compute_index(arc.from, arc.to)] = arc.dist;
    if (_type == AdjustmentMatrixType::symmetric) {
      ret[compute_index(arc.to, arc.from)] = arc.dist;
    }
  }
  for (size_t i = 0; i < _region_count; ++i) { ret[compute_index(i, i)] = 0.0; }
  return ret;
}

size_t AdjustmentMatrix::compute_size() const {
  return _region_count * _region_count;
}

size_t AdjustmentMatrix::compute_index(size_t from, size_t to) const {
  return from * _region_count + to;
}

AdjustmentMatrixType AdjustmentMatrix::determine_matrix_symmetry(
    const std::vector<AdjustmentArc>& arcs, size_t region_count) {
  AdjustmentMatrixType type = AdjustmentMatrixType::invalid;
  if (arcs.size() == ((region_count * (region_count - 1)) / 2)) {
    type = AdjustmentMatrixType::symmetric;
  } else if (arcs.size() == (region_count * (region_count - 1))) {
    type = AdjustmentMatrixType::nonsymmetric;
  }
  if (type == AdjustmentMatrixType::invalid || !is_valid(arcs)) {
    return AdjustmentMatrixType::invalid;
  }
  return type;
}
}  // namespace lagrange
