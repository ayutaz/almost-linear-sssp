#pragma once
// Block lookups in the word model of Sec 2 (arXiv 2611.34911v1, txt l.305-306): "For fixed-size
// blocks of any array, a scan stores the block index and offset for every member. No integer
// division is needed for these lookups."
//
// Integer division and multiplication are not among the constant-time word operations listed
// in l.302-303 (array access, shifts, bit operations, addition, subtraction, comparison).
// BlockLayout is built by one scan that only adds and compares; afterwards every lookup is one
// array access (decision (b) of Q-51; resolves D-18).

#include <cassert>
#include <cstdint>
#include <vector>

namespace sssp::paper::model {

/// Consecutive blocks of `block_size` members over the index range [0, k); the final block may
/// be short. Stores, for every member, its block index and its offset in the block, and for
/// every block its first member (plus the sentinel block_start(num_blocks()) = k).
class BlockLayout {
 public:
  BlockLayout() = default;

  /// One scan over [0, k) (Sec 2 l.305). Throws std::invalid_argument for block_size = 0.
  BlockLayout(std::uint32_t k, std::uint32_t block_size);

  [[nodiscard]] std::uint32_t size() const noexcept {
    return static_cast<std::uint32_t>(block_of_.size());
  }
  [[nodiscard]] std::uint32_t block_size() const noexcept { return block_size_; }
  [[nodiscard]] std::uint32_t num_blocks() const noexcept {
    return start_.empty() ? 0 : static_cast<std::uint32_t>(start_.size() - 1);
  }

  /// Block index of member i (= i / B), one array access. Precondition: i < size().
  [[nodiscard]] std::uint32_t block_of(std::uint32_t i) const noexcept {
    assert(i < size());
    return block_of_[i];
  }
  /// Offset of member i in its block (= i mod B), one array access. Precondition: i < size().
  [[nodiscard]] std::uint32_t offset_of(std::uint32_t i) const noexcept {
    assert(i < size());
    return offset_of_[i];
  }
  /// First member of block b (= b * B); block_start(num_blocks()) = size().
  [[nodiscard]] std::uint32_t block_start(std::uint32_t b) const noexcept {
    assert(b < start_.size());
    return start_[b];
  }
  /// One past the last member of block b (= min(size(), (b + 1) * B)).
  [[nodiscard]] std::uint32_t block_end(std::uint32_t b) const noexcept {
    assert(b + 1 < start_.size());
    return start_[b + 1];
  }

  /// Members visited by the construction scan (diagnostic; 0 when SSSP_DIAGNOSTICS=0).
  [[nodiscard]] std::uint64_t build_steps() const noexcept { return build_steps_; }

 private:
  std::uint32_t block_size_ = 0;
  std::vector<std::uint32_t> block_of_;
  std::vector<std::uint32_t> offset_of_;
  std::vector<std::uint32_t> start_;  ///< num_blocks() + 1 entries; empty for the default layout
  std::uint64_t build_steps_ = 0;
};

}  // namespace sssp::paper::model
