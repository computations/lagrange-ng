#ifndef PERIODS_H
#define PERIODS_H

#include <format>
#include <limits>
#include <cmath>
#include <memory>
#include <ranges>
#include <vector>

namespace lagrange {

struct PeriodDerivative {
  double d_dispersion;
  double d_extinction;

  auto norm() const -> double {
    return d_dispersion * d_dispersion + d_extinction * d_extinction;
  }
};

struct PeriodParams {
  double dispersion_rate;
  double extinction_rate;
  double distance_penalty;
  std::shared_ptr<double[]> adjustment_matrix = nullptr;
  size_t regions;

  void applyDerivative(const PeriodDerivative &d) {
    dispersion_rate += d.d_dispersion;
    extinction_rate += d.d_extinction;
    if (dispersion_rate < 0) { dispersion_rate = 0.0; }
    if (extinction_rate < 0) { extinction_rate = 0.0; }
  }

  auto toString() const -> std::string {
    return std::format("(disp: {}, ext: {})", dispersion_rate, extinction_rate);
  }

  inline auto getDispersionRate(size_t from, size_t to) const -> double {
    return dispersion_rate
           * (adjustment_matrix != nullptr
                  ? std::pow(adjustment_matrix[from * regions + to],
                             distance_penalty)
                  : 1.0);
  }

  inline auto getExtinctionRate() const -> double { return extinction_rate; }

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

  void applyParameters(const PeriodParams& other){
    dispersion_rate = other.dispersion_rate;
    extinction_rate = other.extinction_rate;
    distance_penalty = other.distance_penalty;
  }

  size_t paramCount() const {
    return 2 + (adjustment_matrix == nullptr ? 0.0 : 1.0);
  }
};

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
