#pragma once
// Indexed Fibonacci min-heap used as the boundary heap of GrowthPass
// (Sec 3.2, Algorithm 1 line 5; operation bounds of Lemma 3.3, after Fredman–Tarjan).
//
// Why not baseline::heaps::FibonacciHeap: Algorithm 1 line 5 asks for a *fresh* heap and fresh
// membership identifiers for every growth. A pass may run up to k' growths, so re-allocating a
// k'-node pool per growth costs O(k'^2), and draining the leftover boundary entries by
// extraction costs O(log s) per entry instead of the O(1) per node that the proof of Lemma 3.3
// charges for discarding a heap. This heap keeps one node pool per pass and makes `discard()`
// O(1) by bumping a membership generation.
//
// Keys are compared with `operator<` only (comparison model, Sec 2); ties are allowed.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "sssp/result.hpp"

namespace sssp::paper {

template <class Key>
class FibHeap {
 public:
  using key_type = Key;
  using id_type = std::uint32_t;
  static constexpr id_type kNone = std::numeric_limits<id_type>::max();

  struct Entry {
    Key key;
    id_type id;
  };

  /// Pool for ids [0, capacity). The degree array has the explicit logarithmic size
  /// 3 + 2*ceil(log2(s + 1)) of the proof of Lemma 3.3, with s = capacity >= any heap size.
  explicit FibHeap(std::size_t capacity)
      : nodes_(capacity), by_degree_(3 + 2 * ceil_log2(capacity + 1), kNone) {
    if (capacity >= static_cast<std::size_t>(kNone)) {
      throw std::length_error("FibHeap: too many ids");
    }
    roots_scratch_.reserve(capacity);
    used_degrees_.reserve(by_degree_.size());
  }

  [[nodiscard]] std::size_t capacity() const noexcept { return nodes_.size(); }
  [[nodiscard]] std::size_t size() const noexcept { return size_; }
  [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
  [[nodiscard]] std::size_t degree_array_size() const noexcept { return by_degree_.size(); }
  [[nodiscard]] std::uint64_t comparisons() const noexcept { return comparisons_; }
  [[nodiscard]] std::uint64_t cascading_cuts() const noexcept { return cascading_cuts_; }

  [[nodiscard]] bool contains(id_type id) const noexcept {
    return id < nodes_.size() && nodes_[id].member == generation_;
  }
  [[nodiscard]] const Key& key(id_type id) const {
    require_present(id, "FibHeap::key");
    return nodes_[id].key;
  }

  // Structural accessors (tests and diagnostics); meaningful only for present ids.
  [[nodiscard]] id_type parent_of(id_type id) const {
    require_present(id, "FibHeap::parent_of");
    return nodes_[id].parent;
  }
  [[nodiscard]] std::uint32_t degree(id_type id) const {
    require_present(id, "FibHeap::degree");
    return nodes_[id].degree;
  }
  [[nodiscard]] bool is_marked(id_type id) const {
    require_present(id, "FibHeap::is_marked");
    return nodes_[id].mark;
  }
  [[nodiscard]] std::size_t root_count() const noexcept {
    if (min_ == kNone) return 0;
    std::size_t c = 0;
    id_type r = min_;
    do {
      ++c;
      r = nodes_[r].right;
    } while (r != min_);
    return c;
  }

  /// Lemma 3.3: a new singleton root; updates the minimum pointer. O(1).
  void insert(id_type id, const Key& key) {
    if (id >= nodes_.size()) throw std::out_of_range("FibHeap::insert: id out of range");
    if (contains(id)) throw std::invalid_argument("FibHeap::insert: id already present");
    Node& x = nodes_[id];
    x.key = key;
    x.parent = kNone;
    x.child = kNone;
    x.degree = 0;
    x.mark = false;
    x.member = generation_;
    add_root(id);
    ++size_;
  }

  /// Lemma 3.3: a decrease violating heap order cuts the node to the root list, followed by
  /// cascading cuts of marked ancestors. Amortised O(1). Equal keys are accepted.
  void decrease(id_type id, const Key& key) {
    require_present(id, "FibHeap::decrease");
    Node& x = nodes_[id];
    if (less_keys(x.key, key)) throw std::invalid_argument("FibHeap::decrease: key increase");
    x.key = key;
    const id_type p = x.parent;
    if (p != kNone && less(id, p)) {
      cut(id, p);
      cascading_cut(p);
    }
    if (id != min_ && less(id, min_)) min_ = id;
  }

  /// Lemma 3.3: remove the minimum, promote its children to roots, link roots of equal degree
  /// (larger key becomes the child) through the degree-indexed array, and rebuild the minimum
  /// pointer from the remaining roots. Amortised O(1 + log(s + 1)).
  Entry extract_min() {
    if (size_ == 0) throw std::out_of_range("FibHeap::extract_min on empty heap");
    const id_type z = min_;
    Node& zn = nodes_[z];
    if (zn.child != kNone) {
      id_type c = zn.child;
      do {
        nodes_[c].parent = kNone;
        nodes_[c].mark = false;  // roots are always unmarked
        c = nodes_[c].right;
      } while (c != zn.child);
      splice_into_roots(zn.child);
      zn.child = kNone;
      zn.degree = 0;
    }
    const id_type next = zn.right;
    unlink(z);
    zn.member = 0;
    --size_;
    if (next == z) {
      min_ = kNone;
    } else {
      min_ = next;
      consolidate();
    }
    return Entry{zn.key, z};
  }

  /// Releases every entry in O(1): the next generation treats all nodes as absent
  /// (Algorithm 1 line 5, "a fresh heap and membership identifiers for each growth").
  void discard() noexcept {
    min_ = kNone;
    size_ = 0;
    ++generation_;
    if (generation_ == 0) {  // wrapped after 2^32 discards: clear the stamps once
      for (Node& n : nodes_) n.member = 0;
      generation_ = 1;
    }
  }

  /// Checks root/child links, parent pointers, heap order, degrees, marks (roots unmarked),
  /// the minimum pointer, the size, the degree-array bound and the Fibonacci subtree bound
  /// (a node of degree d roots at least F_{d+2} nodes; proof of Lemma 3.3). Test hook.
  [[nodiscard]] bool validate() const {
    if (min_ == kNone) return size_ == 0;
    if (!contains(min_) || nodes_[min_].parent != kNone) return false;
    std::size_t counted = 0;
    id_type r = min_;
    do {
      const Node& rn = nodes_[r];
      if (!contains(r) || rn.parent != kNone || rn.mark) return false;
      if (less_keys(rn.key, nodes_[min_].key)) return false;
      if (nodes_[rn.right].left != r) return false;
      std::size_t sub = 0;
      if (!validate_tree(r, sub)) return false;
      counted += sub;
      if (counted > size_) return false;
      r = rn.right;
    } while (r != min_);
    return counted == size_;
  }

 private:
  struct Node {
    Key key{};
    id_type parent = kNone;
    id_type child = kNone;
    id_type left = kNone;
    id_type right = kNone;
    std::uint32_t degree = 0;
    std::uint32_t member = 0;  ///< == generation_ iff present
    bool mark = false;
  };

  static std::size_t ceil_log2(std::size_t x) noexcept {
    std::size_t c = 0;
    while ((std::size_t{1} << c) < x) ++c;
    return c;
  }

  void require_present(id_type id, const char* what) const {
    if (!contains(id)) throw std::invalid_argument(std::string(what) + ": id absent");
  }

  bool less_keys(const Key& a, const Key& b) const noexcept {
    SSSP_COUNT(++comparisons_);
    return a < b;
  }
  bool less(id_type a, id_type b) const noexcept { return less_keys(nodes_[a].key, nodes_[b].key); }

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

  /// Splices a circular sibling list (given by any member) into the root list.
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

  /// Makes root y a child of root x; every newly linked child is unmarked (Lemma 3.3).
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
    for (id_type x : roots_scratch_) {
      std::uint32_t d = nodes_[x].degree;
      while (true) {
        if (d >= by_degree_.size()) {
          throw std::logic_error("FibHeap: degree exceeds the Lemma 3.3 bound");
        }
        if (by_degree_[d] == kNone) break;
        id_type y = by_degree_[d];
        if (less(y, x)) std::swap(x, y);
        link(y, x);  // larger-keyed root becomes the child
        by_degree_[d] = kNone;
        ++d;
      }
      if (by_degree_[d] == kNone) used_degrees_.push_back(d);
      by_degree_[d] = x;
    }
    // Rebuild the minimum pointer while collecting the remaining roots.
    min_ = kNone;
    for (const std::uint32_t d : used_degrees_) {
      const id_type x = by_degree_[d];
      by_degree_[d] = kNone;
      if (x == kNone) continue;
      if (min_ == kNone || less(x, min_)) min_ = x;
    }
    used_degrees_.clear();
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

  /// Lemma 3.3: a nonroot is marked on its first child loss and cut on its second, continuing
  /// upwards.
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

  /// Validates the subtree of x and returns its size in `sub`.
  bool validate_tree(id_type x, std::size_t& sub) const {
    const Node& xn = nodes_[x];
    if (!contains(x) || xn.degree >= by_degree_.size()) return false;
    sub = 1;
    std::uint32_t deg = 0;
    if (xn.child != kNone) {
      id_type c = xn.child;
      do {
        const Node& cn = nodes_[c];
        if (cn.parent != x || less_keys(cn.key, xn.key)) return false;
        if (nodes_[cn.right].left != c) return false;
        std::size_t child_size = 0;
        if (!validate_tree(c, child_size)) return false;
        sub += child_size;
        ++deg;
        c = cn.right;
      } while (c != xn.child && deg <= nodes_.size());
    }
    if (deg != xn.degree) return false;
    // Fibonacci bound: size >= F_{d+2} with F_0 = 0, F_1 = 1.
    std::size_t f0 = 0;
    std::size_t f1 = 1;
    for (std::uint32_t i = 0; i < deg + 1; ++i) {
      const std::size_t f2 = f0 + f1;
      f0 = f1;
      f1 = f2;
    }
    return sub >= f1;
  }

  std::vector<Node> nodes_;
  std::vector<id_type> by_degree_;
  std::vector<std::uint32_t> used_degrees_;
  std::vector<id_type> roots_scratch_;
  id_type min_ = kNone;
  std::size_t size_ = 0;
  std::uint32_t generation_ = 1;
  mutable std::uint64_t comparisons_ = 0;
  std::uint64_t cascading_cuts_ = 0;
};

}  // namespace sssp::paper
