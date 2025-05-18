#pragma once
#include <assert.h>

#include <memory>
#include <numeric>
#include <sstream>
#include <vector>

struct CladogenesisParams {
  double copy;
  double sympatry;
  double allopatry;
  double jump;
};

struct PeriodParams {
  double dispersion_rate;
  double extinction_rate;
  std::shared_ptr<std::vector<std::vector<double>>> adjustment_matrix = nullptr;

  auto toString() const -> std::string {
    std::ostringstream os;
    os << "(disp: " << dispersion_rate << ", ext: " << extinction_rate << ")";
    return os.str();
  }

  inline auto getDispersionRate(size_t from, size_t to) const -> double {
    return dispersion_rate
           * (adjustment_matrix != nullptr ? (*adjustment_matrix)[from][to]
                                           : 1.0);
  }

  inline auto getExtinctionRate() const -> double { return extinction_rate; }

  double operator[](size_t i) const {
    if (i == 0) { return dispersion_rate; }
    if (i == 1) { return extinction_rate; }
    return std::numeric_limits<double>::quiet_NaN();
  }

  double &operator[](size_t i) {
    assert(i < 2);
    if (i == 0) { return dispersion_rate; }
    return extinction_rate;
  }
};
