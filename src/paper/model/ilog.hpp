#pragma once
// Integer logarithms in the word model of Sec 2 ("Operations on lengths and indices",
// arXiv 2611.34911v1, txt l.301-306).
//
// Only integers (sizes, indices, counters) are passed here; real lengths never are (Sec 2,
// l.300: no real register supplies a logarithm). The word versions use shifts and comparisons
// only (constant words, at most 64 steps); the table realises the paper's O(K) tabulation.

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace sssp::paper::model {

/// floor(log2 x) for x >= 1, by repeated halving (shift) of the word. Throws std::domain_error
/// for x = 0.
constexpr unsigned floor_log2(std::uint64_t x) {
  if (x == 0) throw std::domain_error("floor_log2(0) is undefined");
  unsigned r = 0;
  while (x > 1) {
    x >>= 1;
    ++r;
  }
  return r;
}

/// ceil(log2 x) for x >= 1, i.e. the least c with x <= 2^c. Throws std::domain_error for x = 0.
/// Used as the bit count of a uniform draw from 0..h-1 (Sec 2, l.308).
constexpr unsigned ceil_log2(std::uint64_t x) {
  if (x == 0) throw std::domain_error("ceil_log2(0) is undefined");
  return x == 1 ? 0 : floor_log2(x - 1) + 1;
}

/// Table of floor(log2 i) for 1 <= i <= K, built in O(K) word operations by the recurrence
/// lg[1] = 0, lg[i] = lg[i >> 1] + 1 (Sec 2, l.303-304: "tabulated ... under halving, starting
/// with i = 1").
class IlogTable {
 public:
  explicit IlogTable(std::size_t range);

  /// K, the largest tabulated index (0 means the table is empty).
  [[nodiscard]] std::size_t range() const noexcept { return lg_.empty() ? 0 : lg_.size() - 1; }

  /// floor(log2 i); throws std::out_of_range unless 1 <= i <= range().
  [[nodiscard]] unsigned floor_log2(std::uint64_t i) const {
    check(i);
    return lg_[i];
  }

  /// ceil(log2 i) = floor(log2(i-1)) + 1 for i >= 2, and 0 for i = 1; same domain as floor_log2.
  [[nodiscard]] unsigned ceil_log2(std::uint64_t i) const {
    check(i);
    return i == 1 ? 0u : static_cast<unsigned>(lg_[i - 1]) + 1u;
  }

 private:
  void check(std::uint64_t i) const {
    if (i == 0 || i > range()) throw std::out_of_range("IlogTable index outside 1..K");
  }

  std::vector<std::uint8_t> lg_;  ///< lg_[0] is unused
};

}  // namespace sssp::paper::model
