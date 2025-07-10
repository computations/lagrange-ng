#ifndef PERIODS_H
#define PERIODS_H

#include <cmath>
#include <format>
#include <limits>
#include <memory>
#include <ranges>
#include <vector>

#include "Common.hpp"
#include "Utils.hpp"

namespace lagrange {

struct PeriodParams {
  double dispersion_rate;
  double extinction_rate;
  double distance_penalty;
  std::shared_ptr<double[]> adjustment_matrix = nullptr;
  size_t regions;
  RangeMask include_area_mask = 0;
  RangeMask exclude_area_mask = 0;

  auto toString() const -> std::string {
    return std::format("(disp: {}, ext: {})", dispersion_rate, extinction_rate);
  }

  auto getDispersionRate(size_t from, size_t to) const -> double {
    if (exclude_area_mask
        && (lagrange_bextr(exclude_area_mask, from)
            || lagrange_bextr(exclude_area_mask, to))) {
      return 0.0;
    }
    return dispersion_rate
           * (adjustment_matrix.get() != nullptr
                  ? std::pow(adjustment_matrix[from * regions + to],
                             distance_penalty)
                  : 1.0);
  }

  inline auto getExtinctionRate(size_t to) const -> double {
    if (include_area_mask && lagrange_bextr(include_area_mask, to)) {
      return 0.0;
    }
    return extinction_rate;
  }

  void applyParameters(const std::vector<double> &x, size_t &index) {
    dispersion_rate = x[index];
    index += 1;
    extinction_rate = x[index];
    index += 1;

    if (adjustment_matrix != nullptr) {
      distance_penalty = x[index];
      index += 1;
    }
  }

  void applyParameters(const PeriodParams &other) {
    dispersion_rate = other.dispersion_rate;
    extinction_rate = other.extinction_rate;
    distance_penalty = other.distance_penalty;
  }

  size_t paramCount() const {
    return 2 + (adjustment_matrix == nullptr ? 0.0 : 1.0);
  }

  bool hasAdjustmentMatrix() const { return adjustment_matrix != nullptr; }
};

struct PeriodSegment {
  size_t index;
  double duration;
  size_t regions;
  size_t max_areas;
};

class PeriodTimes {
 public:
  using const_iterator = std::vector<double>::const_iterator;

  template <typename R>
  PeriodTimes(const R &periods)
    requires(requires(R r) {
      requires std::ranges::range<R>;
      { *r.begin() } -> std::same_as<double>;
    })
      : _end_points{periods.begin(), periods.end()} {
    terminate();
  }

  PeriodTimes() { terminate(); }

  PeriodTimes(const PeriodTimes &) = default;
  PeriodTimes(PeriodTimes &&) = default;

  auto operator=(const PeriodTimes &) -> PeriodTimes & = default;
  auto operator=(PeriodTimes &&) -> PeriodTimes & = default;

  [[nodiscard]] auto begin() const -> const_iterator {
    return _end_points.begin();
  }

  [[nodiscard]] auto end() const -> const_iterator { return _end_points.end(); }

  [[nodiscard]] auto min() const -> double { return 0.0; }

  [[nodiscard]] auto max() const -> double { return _end_points.back(); }

  [[nodiscard]] auto empty() const -> bool { return _end_points.empty(); }

  [[nodiscard]] auto size() const -> size_t { return _end_points.size(); }

  void setMaxAreas(size_t);
  void setRegionCount(size_t);

  [[nodiscard]] auto regions() const -> size_t;
  [[nodiscard]] auto maxAreas() const -> size_t;

 private:
  void terminate() {
    _end_points.emplace_back(std::numeric_limits<double>::infinity());
  }

  size_t _regions = 0;
  size_t _max_areas = 0;

  std::vector<double> _end_points;
};

class PeriodSpan {
 public:
  struct Iterator {
   public:
    using value_type = PeriodSegment;

    Iterator(double start,
             double len,
             size_t index,
             size_t regions,
             size_t max_areas,
             PeriodTimes::const_iterator period) :
        _time{start},
        _length{len},
        _index{index},
        _regions{regions},
        _max_areas{max_areas},
        _period{period} {}

    auto operator*() const -> value_type {
      return {.index = _index,
              .duration = segmentLength(),
              .regions = _regions,
              .max_areas = _max_areas};
    }

    auto operator++() -> Iterator & {
      _index += 1;
      _length -= std::max(segmentLength(), 0.0);
      _time = *_period;
      _period++;
      return *this;
    }

    friend auto operator==(const Iterator &a, const Iterator &b) -> bool {
      return a._length == b._length || a._period == b._period;
    }

    friend auto operator!=(const Iterator &a, const Iterator &b) -> bool {
      return !(a == b);
    }

   private:
    [[nodiscard]] auto segmentLength() const -> double {
      return std::min(*_period - _time, _length);
    }

    double _time;
    double _length;
    size_t _index;
    size_t _regions;
    size_t _max_areas;
    PeriodTimes::const_iterator _period;
  };

  PeriodSpan() = default;

  PeriodSpan(const PeriodTimes &periods) :
      PeriodSpan(periods, 0.0, periods.max()) {}

  PeriodSpan(const PeriodTimes &periods, double start, double length) :
      _start_time{start},
      _length{length},
      _max_areas{periods.maxAreas()},
      _regions{periods.regions()} {
    _begin = periods.begin();
    _end = periods.end();

    _start_index = 0;

    while (*_begin <= _start_time) {
      _begin++;
      _start_index++;
    }

    if (_length == 0.0) { _length = std::numeric_limits<double>::epsilon(); }
  }

  [[nodiscard]] auto begin() const -> Iterator {
    return {_start_time, _length, _start_index, _regions, _max_areas, _begin};
  }

  [[nodiscard]] auto end() const -> Iterator {
    return {0.0, 0.0, 0, _regions, _max_areas, _end};
  }

 private:
  size_t _start_index;
  double _start_time;
  double _length;
  size_t _max_areas;
  size_t _regions;

  PeriodTimes::const_iterator _begin;
  PeriodTimes::const_iterator _end;
};

}  // namespace lagrange

#endif
