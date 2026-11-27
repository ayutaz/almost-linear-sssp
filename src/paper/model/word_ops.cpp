#include "paper/model/word_ops.hpp"

#include <stdexcept>

#include "sssp/result.hpp"

namespace sssp::paper::model {

BlockLayout::BlockLayout(std::uint32_t k, std::uint32_t block_size) : block_size_(block_size) {
  if (block_size == 0) throw std::invalid_argument("BlockLayout: block size must be >= 1");
  block_of_.resize(k);
  offset_of_.resize(k);
  start_.push_back(0);
  // Sec 2 l.305: one scan; the offset advances by one and wraps to 0 at the block size, the
  // block index advances at every wrap. Only addition and comparison are used.
  std::uint32_t block = 0;
  std::uint32_t off = 0;
  for (std::uint32_t i = 0; i < k; ++i) {
    SSSP_COUNT(++build_steps_);
    if (off == block_size) {
      off = 0;
      ++block;
      start_.push_back(i);
    }
    block_of_[i] = block;
    offset_of_[i] = off;
    ++off;
  }
  if (k != 0) start_.push_back(k);
}

}  // namespace sssp::paper::model
