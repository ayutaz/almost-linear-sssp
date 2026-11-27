#pragma once
// Binary min-heap of (key, value) entries without decrease-key ("lazy deletion": callers push a
// new entry and discard stale ones on pop). Same semantics as std::priority_queue with
// std::greater, but implemented here so that storage can be reserved up front (no allocation
// during a timed run) and key comparisons can be counted.

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "sssp/result.hpp"

namespace sssp::baseline::heaps {

template <class K, class V = std::uint32_t>
class BinaryHeap {
 public:
  struct Entry {
    K key;
    V value;
  };

  BinaryHeap() = default;
  explicit BinaryHeap(std::size_t capacity) { reserve(capacity); }

  void reserve(std::size_t capacity) { data_.reserve(capacity); }
  [[nodiscard]] std::size_t capacity() const noexcept { return data_.capacity(); }
  [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }
  [[nodiscard]] bool empty() const noexcept { return data_.empty(); }
  /// Removes all entries; keeps the reserved storage.
  void clear() noexcept { data_.clear(); }
  /// Number of key comparisons performed so far (0 when SSSP_DIAGNOSTICS=0).
  [[nodiscard]] std::uint64_t comparisons() const noexcept { return comparisons_; }

  void push(K key, V value) {
    data_.push_back(Entry{key, value});
    sift_up(data_.size() - 1);
  }

  [[nodiscard]] const Entry& top() const {
    if (data_.empty()) throw std::out_of_range("BinaryHeap::top on empty heap");
    return data_.front();
  }

  Entry pop() {
    if (data_.empty()) throw std::out_of_range("BinaryHeap::pop on empty heap");
    const Entry result = data_.front();
    const Entry last = data_.back();
    data_.pop_back();
    if (!data_.empty()) sift_down(0, last);
    return result;
  }

 private:
  bool less(const K& a, const K& b) noexcept {
    SSSP_COUNT(++comparisons_);
    return a < b;
  }

  void sift_up(std::size_t i) {
    const Entry x = data_[i];
    while (i > 0) {
      const std::size_t p = (i - 1) / 2;
      if (!less(x.key, data_[p].key)) break;
      data_[i] = data_[p];
      i = p;
    }
    data_[i] = x;
  }

  /// Places x into the hole at index i, moving smaller children up.
  void sift_down(std::size_t i, const Entry& x) {
    const std::size_t n = data_.size();
    for (;;) {
      std::size_t c = 2 * i + 1;
      if (c >= n) break;
      if (c + 1 < n && less(data_[c + 1].key, data_[c].key)) ++c;
      if (!less(data_[c].key, x.key)) break;
      data_[i] = data_[c];
      i = c;
    }
    data_[i] = x;
  }

  std::vector<Entry> data_;
  std::uint64_t comparisons_ = 0;
};

}  // namespace sssp::baseline::heaps
