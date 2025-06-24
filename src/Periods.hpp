#ifndef PERIODS_H
#define PERIODS_H

#include <limits>
#include <ranges>
#include <vector>

namespace lagrange {

struct PeriodSegment {
  size_t index;
  double duration;
  size_t regions;
  size_t max_areas;
};

class Periods {
 public:
  using const_iterator = std::vector<double>::const_iterator;

  Periods(const std::ranges::range auto &periods) :
      _periods{periods.begin(), periods.end()} {
    terminate();
  }

  Periods() { terminate(); }

  Periods(const Periods &) = default;
  Periods(Periods &&) = default;

  auto operator=(const Periods &) -> Periods & = default;
  auto operator=(Periods &&) -> Periods & = default;

  [[nodiscard]] auto begin() const -> const_iterator {
    return _periods.begin();
  }

  [[nodiscard]] auto end() const -> const_iterator { return _periods.end(); }

  [[nodiscard]] auto min() const -> double { return _periods.front(); }

  [[nodiscard]] auto max() const -> double { return _periods.back(); }

  [[nodiscard]] auto empty() const -> bool { return _periods.empty(); }

  [[nodiscard]] auto size() const -> size_t { return _periods.size(); }

  void setMaxAreas(size_t);
  void setRegionCount(size_t);

  [[nodiscard]] auto regions() const -> size_t;
  [[nodiscard]] auto maxAreas() const -> size_t;

 private:
  void terminate() {
    _periods.emplace_back(std::numeric_limits<double>::infinity());
  }

  size_t _regions = 0;
  size_t _max_areas = 0;

  std::vector<double> _periods;
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
             Periods::const_iterator period) :
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
    Periods::const_iterator _period;
  };

  PeriodSpan() = default;

  PeriodSpan(const Periods &periods) :
      PeriodSpan(periods, 0.0, periods.max()) {}

  PeriodSpan(const Periods &periods, double start, double length) :
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

  Periods::const_iterator _begin;
  Periods::const_iterator _end;
};

}  // namespace lagrange

#endif
