#ifndef PERIODS_H
#define PERIODS_H

#include <iterator>
#include <limits>
#include <vector>

struct PeriodSegment {
  size_t index;
  double duration;
};

class Periods {
 public:
  using const_iterator = std::vector<double>::const_iterator;
  Periods(const std::vector<double> &periods) : _periods{periods} {
    terminate();
  }

  Periods() : _periods{} { terminate(); }
  Periods(const Periods &) = default;
  Periods(Periods &&) = default;

  Periods &operator=(const Periods &) = default;
  Periods &operator=(Periods &&) = default;

  const_iterator begin() const { return _periods.begin(); }
  const_iterator end() const { return _periods.end(); }

  double min() const { return _periods.front(); }
  double max() const { return _periods.back(); }

  bool empty() const { return _periods.empty(); }

  size_t size() const { return _periods.size(); }

 private:
  void terminate() {
    _periods.emplace_back(std::numeric_limits<double>::infinity());
  }

  std::vector<double> _periods;
};

class PeriodSpan {
 public:
  struct Iterator {
   public:
    using value_type = PeriodSegment;

    Iterator(double start, double len, size_t index,
             Periods::const_iterator period)
        : _time{start}, _length{len}, _index{index}, _period{period} {}

    value_type operator*() const {
      return {.index = _index, .duration = segment_length()};
    }

    Iterator &operator++() {
      _index += 1;
      _length -= std::max(segment_length(), 0.0);
      _time = *_period;
      _period++;
      return *this;
    }

    friend bool operator==(const Iterator &a, const Iterator &b) {
      return a._length == b._length || a._period == b._period;
    }
    friend bool operator!=(const Iterator &a, const Iterator &b) {
      return !(a == b);
    }

   private:
    double segment_length() const {
      return std::min(*_period - _time, _length);
    }

    double _time;
    double _length;
    size_t _index;
    Periods::const_iterator _period;
  };
  PeriodSpan() = default;

  PeriodSpan(const Periods &periods, double start, double length)
      : _start_time{start}, _length{length} {
    _begin = periods.begin();
    _end = periods.end();

    while (*_begin <= _start_time) {
      _begin++;
      _start_index++;
    }
  }

  Iterator begin() const {
    return Iterator(_start_time, _length, _start_index, _begin);
  }

  Iterator end() const { return Iterator(0.0, 0.0, 0, _end); }

 private:
  size_t _start_index;
  double _start_time;
  double _length;

  Periods::const_iterator _begin;
  Periods::const_iterator _end;
};

#endif
