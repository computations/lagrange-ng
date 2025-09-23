#pragma once
#include <memory>

#include "CSV.hpp"
#include "Operation.hpp"
#include "Workspace.hpp"

namespace lagrange {
template <StreamableGoal T>
class StreamingGoal {
  void eval(const std::shared_ptr<Workspace>& ws) {
    _goal->eval(ws);
    std::lock_guard<std::mutex> lock{*_lock};
    for (auto a : _goal->result()) { write_csv_row(*_outfile, a); }
  }

 private:
  std::shared_ptr<T> _goal;
  std::shared_ptr<std::mutex> _lock;
  std::shared_ptr<std::ofstream> _outfile;
};
}  // namespace lagrange
