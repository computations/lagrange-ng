#include "AdjustmentMatrix.hpp"

#include <algorithm>
#include <functional>
#include <ranges>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "CSV.hpp"
#include "Utils.hpp"

namespace lagrange {

using namespace std::string_view_literals;

#if defined(__cpp_lib_ranges_zip) \
    && (!defined(__clang_major__) || __clang_major__ > 18)
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
#endif

/* assume that the matrix is sorted */
/* need to check that:
 * - No duplicates
 * - All arcs are found (both symmetric and non)
 */
static bool is_valid(std::span<const AdjustmentArc> arcs) {
#if defined(__cpp_lib_ranges_zip) \
    && (!defined(__clang_major__) || __clang_major__ > 18)
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
#else
  std::unordered_set<
      std::pair<decltype(AdjustmentArc::from), decltype(AdjustmentArc::to)> >
      arc_set;
  for (auto a : arcs) { arc_set.insert({a.from, a.to}); }

  auto expected_sym_size = (arcs.size() - 1) * (arcs.size()) / 2;
  auto expected_nonsym_size = arcs.size() * arcs.size();

  return arc_set.size() == expected_sym_size
         || arc_set.size() == expected_nonsym_size;
#endif
}

void AdjustmentMatrix::read_arcs(CSVReader& reader,
                                 const std::vector<std::string>& area_names) {
  _region_count = area_names.size();
#ifdef __cpp_lib_ranges_zip
  auto reverse_area_name_map =
      std::views::zip(area_names, std::views::iota(0ul, area_names.size()))
      | std::ranges::to<std::unordered_map>();
#else
  std::unordered_map<std::string, size_t> reverse_area_name_map;
  for (size_t i = 0; i < area_names.size(); ++i) {
    reverse_area_name_map[area_names[i]] = i;
  }
#endif

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
