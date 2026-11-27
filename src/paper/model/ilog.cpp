#include "paper/model/ilog.hpp"

#include "sssp/result.hpp"

namespace sssp::paper::model {

IlogTable::IlogTable(std::size_t range) {
  if (range == 0) return;
  lg_.assign(range + 1, 0);
  SSSP_COUNT(++build_steps_);  // lg[1] = 0
  // Sec 2, l.303-304: lg[1] = 0 and lg[i] = lg[i >> 1] + 1 for i >= 2; one step per index.
  for (std::size_t i = 2; i <= range; ++i) {
    SSSP_COUNT(++build_steps_);
    lg_[i] = static_cast<std::uint8_t>(lg_[i >> 1] + 1);
  }
}

}  // namespace sssp::paper::model
