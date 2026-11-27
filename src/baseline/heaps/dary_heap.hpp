#pragma once
// Indexed d-ary min-heap with decrease-key over ids [0, n). All storage (heap array, position
// map, key array) is allocated once in the constructor, so a Dijkstra run never allocates.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "sssp/result.hpp"

namespace sssp::baseline::heaps {

template <class K, unsigned D>
class IndexedDaryHeap {
  static_assert(D >= 2, "arity must be at least 2");

 public:
  using key_type = K;
  using id_type = std::uint32_t;
  static constexpr id_type kNone = std::numeric_limits<id_type>::max();

  struct Entry {
    K key;
    id_type id;
  };

  explicit IndexedDaryHeap(std::size_t n) : pos_(n, kNone), key_(n) { heap_.reserve(n); }

  [[nodiscard]] std::size_t capacity() const noexcept { return pos_.size(); }
  [[nodiscard]] std::size_t size() const noexcept { return heap_.size(); }
  [[nodiscard]] bool empty() const noexcept { return heap_.empty(); }
  [[nodiscard]] bool contains(id_type id) const noexcept {
    return id < pos_.size() && pos_[id] != kNone;
  }
  /// Current key of an id (the last key it had if it was popped).
  [[nodiscard]] K key(id_type id) const { return key_.at(id); }
  [[nodiscard]] std::uint64_t comparisons() const noexcept { return comparisons_; }

  void push(id_type id, K key) {
    check_id(id);
    if (pos_[id] != kNone) throw std::invalid_argument("IndexedDaryHeap::push: id already present");
    key_[id] = key;
    heap_.push_back(id);
    sift_up(heap_.size() - 1);
  }

  void decrease_key(id_type id, K key) {
    check_id(id);
    if (pos_[id] == kNone) throw std::invalid_argument("IndexedDaryHeap::decrease_key: id absent");
    if (key_[id] < key) throw std::invalid_argument("IndexedDaryHeap::decrease_key: key increase");
    key_[id] = key;
    sift_up(pos_[id]);
  }

  Entry pop() {
    if (heap_.empty()) throw std::out_of_range("IndexedDaryHeap::pop on empty heap");
    const id_type top = heap_.front();
    const id_type last = heap_.back();
    heap_.pop_back();
    pos_[top] = kNone;
    if (!heap_.empty()) sift_down(0, last);
    return Entry{key_[top], top};
  }

  /// Checks heap order and the position map (test hook; O(n)).
  [[nodiscard]] bool validate() const {
    std::size_t present = 0;
    for (id_type v = 0; v < pos_.size(); ++v) {
      if (pos_[v] == kNone) continue;
      ++present;
      if (pos_[v] >= heap_.size() || heap_[pos_[v]] != v) return false;
    }
    if (present != heap_.size()) return false;
    for (std::size_t i = 1; i < heap_.size(); ++i) {
      if (key_[heap_[i]] < key_[heap_[(i - 1) / D]]) return false;
    }
    return true;
  }

 private:
  void check_id(id_type id) const {
    if (id >= pos_.size()) throw std::out_of_range("IndexedDaryHeap: id out of range");
  }

  bool less(id_type a, id_type b) noexcept {
    SSSP_COUNT(++comparisons_);
    return key_[a] < key_[b];
  }

  void place(std::size_t i, id_type id) noexcept {
    heap_[i] = id;
    pos_[id] = static_cast<id_type>(i);
  }

  void sift_up(std::size_t i) {
    const id_type x = heap_[i];
    while (i > 0) {
      const std::size_t p = (i - 1) / D;
      if (!less(x, heap_[p])) break;
      place(i, heap_[p]);
      i = p;
    }
    place(i, x);
  }

  /// Places x into the hole at index i, moving the smallest child up while it beats x.
  void sift_down(std::size_t i, id_type x) {
    const std::size_t n = heap_.size();
    for (;;) {
      const std::size_t first = D * i + 1;
      if (first >= n) break;
      const std::size_t end = first + D < n ? first + D : n;
      std::size_t best = first;
      for (std::size_t c = first + 1; c < end; ++c) {
        if (less(heap_[c], heap_[best])) best = c;
      }
      if (!less(heap_[best], x)) break;
      place(i, heap_[best]);
      i = best;
    }
    place(i, x);
  }

  std::vector<id_type> heap_;
  std::vector<id_type> pos_;
  std::vector<K> key_;
  std::uint64_t comparisons_ = 0;
};

}  // namespace sssp::baseline::heaps
