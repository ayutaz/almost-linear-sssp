#pragma once
// Sec 3.3 "Disjoint Sets for Ordered Merge Events", Lemma 3.6.
//
// Union by rank with full path compression on k items. Merge events are processed in the order
// the caller issues them; ordering them is the caller's cost (Sec 3.3). Only integer word
// operations are used (Sec 2).
//
// Lemma 3.6: initialisation plus Q unions/finds cost O(Q + k*tau(k)). The proof charges at most
// three parent-pointer traversals to each find (last edge + two band crossings) and the rest to
// items inside the three rank bands [0,b1), [b1,b2), [b2,b3). find_steps() counts every
// traversal so that callers and property tests can check this bound by counted operations.

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "sssp/result.hpp"

namespace sssp::paper {

class DisjointSets {
 public:
  using Item = std::uint32_t;

  /// Sec 3.3: initially each of the k items is a root of rank zero.
  explicit DisjointSets(std::size_t k) : parent_(k), rank_(k, 0), sets_(k) {
    if (k >= static_cast<std::size_t>(std::numeric_limits<Item>::max())) {
      throw std::length_error("DisjointSets: too many items");
    }
    for (std::size_t i = 0; i < k; ++i) parent_[i] = static_cast<Item>(i);
  }

  /// Sec 3.3: return the root of x's parent chain and reset the parent of every nonroot on the
  /// chain to that root.
  Item find(Item x) {
    assert(x < parent_.size());
    SSSP_COUNT(++find_calls_);
    Item root = x;
    while (parent_[root] != root) {
      root = parent_[root];
      SSSP_COUNT(++find_steps_);
    }
    while (parent_[x] != root) {
      const Item next = parent_[x];
      parent_[x] = root;
      x = next;
    }
    return root;
  }

  /// Sec 3.3: find both roots; do nothing if they agree. Otherwise attach the root of smaller
  /// rank to that of larger rank; on a tie attach b's root below a's root and raise a's rank.
  /// Returns true iff two different sets were merged.
  bool unite(Item a, Item b) {
    Item ra = find(a);
    Item rb = find(b);
    if (ra == rb) return false;
    if (rank_[ra] < rank_[rb]) {
      const Item t = ra;
      ra = rb;
      rb = t;
    }
    parent_[rb] = ra;
    if (rank_[ra] == rank_[rb]) ++rank_[ra];
    --sets_;
    SSSP_COUNT(++unions_);
    return true;
  }

  /// Number of items k.
  [[nodiscard]] std::size_t size() const noexcept { return parent_.size(); }
  /// Number of disjoint sets currently represented.
  [[nodiscard]] std::size_t num_sets() const noexcept { return sets_; }

  /// Inspection (tests, invariant checks): current parent pointer and rank. No compression.
  [[nodiscard]] Item parent(Item x) const { return parent_[x]; }
  [[nodiscard]] std::uint32_t rank(Item x) const { return rank_[x]; }

  /// Diagnostic counters (zero when SSSP_DIAGNOSTICS=0).
  [[nodiscard]] std::uint64_t find_calls() const noexcept { return find_calls_; }
  [[nodiscard]] std::uint64_t find_steps() const noexcept { return find_steps_; }
  [[nodiscard]] std::uint64_t unions() const noexcept { return unions_; }

  /// Adds "dsu.find_calls", "dsu.find_steps" and "dsu.unions" to diag.counters.
  void add_counters(Diagnostics& diag) const {
    diag.counters["dsu.find_calls"] += find_calls_;
    diag.counters["dsu.find_steps"] += find_steps_;
    diag.counters["dsu.unions"] += unions_;
  }

 private:
  std::vector<Item> parent_;
  std::vector<std::uint32_t> rank_;  ///< <= floor(log2 k) (proof of Lemma 3.6)
  std::size_t sets_ = 0;
  std::uint64_t find_calls_ = 0;
  std::uint64_t find_steps_ = 0;
  std::uint64_t unions_ = 0;
};

}  // namespace sssp::paper
