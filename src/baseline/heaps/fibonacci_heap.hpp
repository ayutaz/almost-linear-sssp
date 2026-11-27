#pragma once
// Indexed Fibonacci min-heap (Fredman & Tarjan 1987; CLRS ch. 19) over ids [0, n).
// Every id owns one node of a pool allocated in the constructor; the consolidation buffers are
// also allocated once, so no operation allocates.
//
// Amortised costs: push / decrease-key O(1), pop O(log n).

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "sssp/result.hpp"

namespace sssp::baseline::heaps {

template <class K>
class FibonacciHeap {
 public:
  using key_type = K;
  using id_type = std::uint32_t;
  static constexpr id_type kNone = std::numeric_limits<id_type>::max();

  struct Entry {
    K key;
    id_type id;
  };

  explicit FibonacciHeap(std::size_t n) : nodes_(n) { roots_scratch_.reserve(n); }

  [[nodiscard]] std::size_t capacity() const noexcept { return nodes_.size(); }
  [[nodiscard]] std::size_t size() const noexcept { return size_; }
  [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
  [[nodiscard]] bool contains(id_type id) const noexcept {
    return id < nodes_.size() && nodes_[id].in_heap;
  }
  [[nodiscard]] K key(id_type id) const { return nodes_.at(id).key; }
  [[nodiscard]] std::uint64_t comparisons() const noexcept { return comparisons_; }

  // Structural accessors (tests and diagnostics).
  [[nodiscard]] id_type parent_of(id_type id) const { return nodes_.at(id).parent; }
  [[nodiscard]] std::uint32_t degree(id_type id) const { return nodes_.at(id).degree; }
  [[nodiscard]] bool is_marked(id_type id) const { return nodes_.at(id).mark; }
  [[nodiscard]] std::uint64_t cascading_cuts() const noexcept { return cascading_cuts_; }
  [[nodiscard]] std::size_t root_count() const noexcept {
    std::size_t c = 0;
    if (min_ == kNone) return 0;
    id_type r = min_;
    do {
      ++c;
      r = nodes_[r].right;
    } while (r != min_);
    return c;
  }

  void push(id_type id, K key) {
    check_id(id);
    Node& x = nodes_[id];
    if (x.in_heap) throw std::invalid_argument("FibonacciHeap::push: id already present");
    x.key = key;
    x.parent = kNone;
    x.child = kNone;
    x.degree = 0;
    x.mark = false;
    x.in_heap = true;
    add_root(id);
    ++size_;
  }

  void decrease_key(id_type id, K key) {
    check_id(id);
    Node& x = nodes_[id];
    if (!x.in_heap) throw std::invalid_argument("FibonacciHeap::decrease_key: id absent");
    if (x.key < key) throw std::invalid_argument("FibonacciHeap::decrease_key: key increase");
    x.key = key;
    const id_type p = x.parent;
    if (p != kNone && less(id, p)) {
      cut(id, p);
      cascading_cut(p);
    }
    if (id != min_ && less(id, min_)) min_ = id;
  }

  Entry pop() {
    if (size_ == 0) throw std::out_of_range("FibonacciHeap::pop on empty heap");
    const id_type z = min_;
    Node& zn = nodes_[z];
    // Move every child of z to the root list.
    if (zn.child != kNone) {
      id_type c = zn.child;
      do {
        nodes_[c].parent = kNone;
        nodes_[c].mark = false;
        c = nodes_[c].right;
      } while (c != zn.child);
      splice_into_roots(zn.child);
      zn.child = kNone;
      zn.degree = 0;
    }
    const id_type next = zn.right;
    unlink(z);
    zn.in_heap = false;
    --size_;
    if (next == z) {
      min_ = kNone;
    } else {
      min_ = next;
      consolidate();
    }
    return Entry{zn.key, z};
  }

  /// Checks heap order, parent/child/sibling links, degrees, marks and size (test hook).
  [[nodiscard]] bool validate() const {
    std::size_t counted = 0;
    if (min_ == kNone) return size_ == 0;
    if (!nodes_[min_].in_heap || nodes_[min_].parent != kNone) return false;
    id_type r = min_;
    do {
      if (nodes_[r].parent != kNone || nodes_[r].key < nodes_[min_].key) return false;
      if (!validate_tree(r, counted)) return false;
      r = nodes_[r].right;
    } while (r != min_);
    return counted == size_;
  }

 private:
  struct Node {
    K key{};
    id_type parent = kNone;
    id_type child = kNone;
    id_type left = kNone;
    id_type right = kNone;
    std::uint32_t degree = 0;
    bool mark = false;
    bool in_heap = false;
  };

  /// Degrees are at most log_phi(n) < 1.4405 log2(n) + 2 < 64 for n < 2^32 ids.
  static constexpr std::size_t kMaxDegree = 64;

  void check_id(id_type id) const {
    if (id >= nodes_.size()) throw std::out_of_range("FibonacciHeap: id out of range");
  }

  bool less(id_type a, id_type b) noexcept {
    SSSP_COUNT(++comparisons_);
    return nodes_[a].key < nodes_[b].key;
  }

  /// Inserts a single node into the root list and updates the minimum.
  void add_root(id_type id) {
    if (min_ == kNone) {
      nodes_[id].left = nodes_[id].right = id;
      min_ = id;
      return;
    }
    insert_after(min_, id);
    if (less(id, min_)) min_ = id;
  }

  void insert_after(id_type at, id_type id) noexcept {
    const id_type r = nodes_[at].right;
    nodes_[id].left = at;
    nodes_[id].right = r;
    nodes_[r].left = id;
    nodes_[at].right = id;
  }

  void unlink(id_type id) noexcept {
    const id_type l = nodes_[id].left;
    const id_type r = nodes_[id].right;
    nodes_[l].right = r;
    nodes_[r].left = l;
    nodes_[id].left = nodes_[id].right = id;
  }

  /// Splices a whole circular sibling list (given by any member) into the root list.
  void splice_into_roots(id_type list) noexcept {
    if (min_ == kNone) {
      min_ = list;
      return;
    }
    const id_type a = min_;
    const id_type a_right = nodes_[a].right;
    const id_type b_left = nodes_[list].left;
    nodes_[a].right = list;
    nodes_[list].left = a;
    nodes_[b_left].right = a_right;
    nodes_[a_right].left = b_left;
  }

  /// Makes root y a child of root x.
  void link(id_type y, id_type x) noexcept {
    unlink(y);
    Node& xn = nodes_[x];
    if (xn.child == kNone) {
      xn.child = y;
    } else {
      insert_after(xn.child, y);
    }
    nodes_[y].parent = x;
    nodes_[y].mark = false;
    ++xn.degree;
  }

  void consolidate() {
    roots_scratch_.clear();
    id_type r = min_;
    do {
      roots_scratch_.push_back(r);
      r = nodes_[r].right;
    } while (r != min_);
    by_degree_.fill(kNone);
    for (id_type x : roots_scratch_) {
      std::uint32_t d = nodes_[x].degree;
      while (by_degree_[d] != kNone) {
        id_type y = by_degree_[d];
        if (less(y, x)) {
          const id_type t = x;
          x = y;
          y = t;
        }
        link(y, x);
        by_degree_[d] = kNone;
        ++d;
      }
      by_degree_[d] = x;
    }
    min_ = kNone;
    for (const id_type x : by_degree_) {
      if (x == kNone) continue;
      if (min_ == kNone || less(x, min_)) min_ = x;
    }
  }

  /// Moves x (a child of p) to the root list.
  void cut(id_type x, id_type p) noexcept {
    Node& pn = nodes_[p];
    if (pn.child == x) pn.child = nodes_[x].right == x ? kNone : nodes_[x].right;
    unlink(x);
    --pn.degree;
    nodes_[x].parent = kNone;
    nodes_[x].mark = false;
    insert_after(min_, x);
  }

  void cascading_cut(id_type y) noexcept {
    for (id_type z = nodes_[y].parent; z != kNone; z = nodes_[y].parent) {
      if (!nodes_[y].mark) {
        nodes_[y].mark = true;
        return;
      }
      cut(y, z);
      ++cascading_cuts_;
      y = z;
    }
  }

  bool validate_tree(id_type x, std::size_t& counted) const {
    const Node& xn = nodes_[x];
    if (!xn.in_heap) return false;
    ++counted;
    if (xn.child == kNone) return xn.degree == 0;
    std::uint32_t deg = 0;
    id_type c = xn.child;
    do {
      const Node& cn = nodes_[c];
      if (cn.parent != x || cn.key < xn.key) return false;
      if (nodes_[cn.right].left != c) return false;
      if (!validate_tree(c, counted)) return false;
      ++deg;
      c = cn.right;
    } while (c != xn.child && deg <= nodes_.size());
    return deg == xn.degree;
  }

  std::vector<Node> nodes_;
  std::vector<id_type> roots_scratch_;
  std::array<id_type, kMaxDegree> by_degree_{};
  id_type min_ = kNone;
  std::size_t size_ = 0;
  std::uint64_t comparisons_ = 0;
  std::uint64_t cascading_cuts_ = 0;
};

}  // namespace sssp::baseline::heaps
