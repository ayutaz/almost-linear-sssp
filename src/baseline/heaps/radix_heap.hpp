#pragma once
// Indexed radix heap (Ahuja, Mehlhorn, Orlin & Tarjan 1990) for monotone 64-bit unsigned keys:
// every key ever inserted must be >= the last popped key. Item x lives in bucket
// bit_width(key(x) XOR last) in [0, 64]; bucket 0 holds keys equal to `last`.
// Buckets are intrusive doubly-linked lists over a node pool allocated in the constructor,
// so no operation allocates. Decrease-key moves the item to its new bucket in O(1).
//
// pop: if bucket 0 is empty, the first non-empty bucket i is scanned for its minimum, which
// becomes the new `last`; every item of bucket i then moves to a bucket < i. Each item moves at
// most 64 times, so a Dijkstra run costs O(m + n * 64).

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "sssp/result.hpp"

namespace sssp::baseline::heaps {

class RadixHeap {
 public:
  using key_type = std::uint64_t;
  using id_type = std::uint32_t;
  static constexpr id_type kNone = std::numeric_limits<id_type>::max();
  static constexpr std::size_t kBuckets = 65;

  struct Entry {
    key_type key;
    id_type id;
  };

  explicit RadixHeap(std::size_t n) : nodes_(n) { heads_.fill(kNone); }

  [[nodiscard]] std::size_t capacity() const noexcept { return nodes_.size(); }
  [[nodiscard]] std::size_t size() const noexcept { return size_; }
  [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
  [[nodiscard]] bool contains(id_type id) const noexcept {
    return id < nodes_.size() && nodes_[id].bucket != kAbsent;
  }
  [[nodiscard]] key_type key(id_type id) const { return nodes_.at(id).key; }
  /// The last popped key (lower bound for every future key).
  [[nodiscard]] key_type last() const noexcept { return last_; }
  [[nodiscard]] std::uint64_t comparisons() const noexcept { return comparisons_; }

  void push(id_type id, key_type key) {
    check_id(id);
    if (nodes_[id].bucket != kAbsent) throw std::invalid_argument("RadixHeap::push: id present");
    check_monotone(key);
    nodes_[id].key = key;
    insert(id, bucket_of(key));
    ++size_;
  }

  void decrease_key(id_type id, key_type key) {
    check_id(id);
    Node& x = nodes_[id];
    if (x.bucket == kAbsent) throw std::invalid_argument("RadixHeap::decrease_key: id absent");
    if (x.key < key) throw std::invalid_argument("RadixHeap::decrease_key: key increase");
    check_monotone(key);
    x.key = key;
    const auto b = bucket_of(key);
    if (b != x.bucket) {
      remove(id);
      insert(id, b);
    }
  }

  Entry pop() {
    if (size_ == 0) throw std::out_of_range("RadixHeap::pop on empty heap");
    if (heads_[0] == kNone) redistribute();
    const id_type x = heads_[0];
    remove(x);
    nodes_[x].bucket = kAbsent;
    --size_;
    return Entry{nodes_[x].key, x};
  }

  /// Checks bucket membership, list links and size (test hook; O(n)).
  [[nodiscard]] bool validate() const {
    std::size_t counted = 0;
    for (std::size_t b = 0; b < kBuckets; ++b) {
      id_type prev = kNone;
      for (id_type x = heads_[b]; x != kNone; x = nodes_[x].next) {
        const Node& xn = nodes_[x];
        if (xn.bucket != b || xn.prev != prev || xn.key < last_ || bucket_of(xn.key) != b) {
          return false;
        }
        prev = x;
        if (++counted > size_) return false;
      }
    }
    return counted == size_;
  }

 private:
  static constexpr std::uint8_t kAbsent = 0xFF;

  struct Node {
    key_type key = 0;
    id_type prev = kNone;
    id_type next = kNone;
    std::uint8_t bucket = kAbsent;
  };

  void check_id(id_type id) const {
    if (id >= nodes_.size()) throw std::out_of_range("RadixHeap: id out of range");
  }
  void check_monotone(key_type key) const {
    if (key < last_) throw std::invalid_argument("RadixHeap: key below the last popped key");
  }

  [[nodiscard]] std::uint8_t bucket_of(key_type key) const noexcept {
    return static_cast<std::uint8_t>(std::bit_width(key ^ last_));
  }

  void insert(id_type id, std::uint8_t b) noexcept {
    Node& x = nodes_[id];
    x.bucket = b;
    x.prev = kNone;
    x.next = heads_[b];
    if (heads_[b] != kNone) nodes_[heads_[b]].prev = id;
    heads_[b] = id;
  }

  void remove(id_type id) noexcept {
    Node& x = nodes_[id];
    if (x.prev != kNone) {
      nodes_[x.prev].next = x.next;
    } else {
      heads_[x.bucket] = x.next;
    }
    if (x.next != kNone) nodes_[x.next].prev = x.prev;
    x.prev = x.next = kNone;
  }

  /// Precondition: bucket 0 is empty and the heap is not.
  void redistribute() {
    std::size_t i = 1;
    while (heads_[i] == kNone) ++i;
    id_type x = heads_[i];
    key_type m = nodes_[x].key;
    for (x = nodes_[x].next; x != kNone; x = nodes_[x].next) {
      SSSP_COUNT(++comparisons_);
      if (nodes_[x].key < m) m = nodes_[x].key;
    }
    last_ = m;
    x = heads_[i];
    heads_[i] = kNone;
    while (x != kNone) {
      const id_type next = nodes_[x].next;
      insert(x, bucket_of(nodes_[x].key));  // always a bucket < i
      x = next;
    }
  }

  std::vector<Node> nodes_;
  std::array<id_type, kBuckets> heads_{};
  key_type last_ = 0;
  std::size_t size_ = 0;
  std::uint64_t comparisons_ = 0;
};

}  // namespace sssp::baseline::heaps
