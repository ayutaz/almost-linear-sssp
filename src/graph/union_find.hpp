#pragma once
// Internal union-find shared by src/graph (P2-UNIONFIND-DEDUP); not part of the public API.
//
// Used by connected_components (edge_list.cpp) and by count_components, connect_components and
// the road_like connectivity repair (generators.cpp). Those callers depend only on the
// partition: on the value returned by unite() and on whether two find() results are equal.
// They never use which element is a root (components are numbered by scanning the vertices in
// increasing order), so the root rule below does not change any of their outputs.
//
// Deliberately separate from (do not merge with):
// - src/paper/primitives/dsu.hpp DisjointSets: union by rank with full path compression and
//   operation counters, as required by Sec 3.3, Lemma 3.6 (Mode A);
// - the union-finds inside the paper's reference/oracle code and the tests, which stay
//   independent of the code they check.

#include <concepts>
#include <cstddef>
#include <utility>
#include <vector>

#include "sssp/graph.hpp"

namespace sssp::detail {

/// Disjoint sets over the elements [0, n): union by size with path halving, so a sequence of
/// operations costs amortised O(alpha(n)) each. Deterministic: on unite(a, b) the root of the
/// larger set becomes the root of the union; on equal sizes, the root of a's set.
/// Index must be able to represent every element (Vertex for vertex sets of a graph).
template <std::unsigned_integral Index = Vertex>
class UnionFind {
 public:
  explicit UnionFind(std::size_t n) : parent_(n), size_(n, 1) {
    for (std::size_t i = 0; i < n; ++i) parent_[i] = static_cast<Index>(i);
  }

  /// The root of x's set (x < n). Path halving: every visited element skips its parent.
  Index find(Index x) {
    while (parent_[x] != x) {
      parent_[x] = parent_[parent_[x]];
      x = parent_[x];
    }
    return x;
  }

  /// Merges the sets of a and b; returns true iff they were different sets.
  bool unite(Index a, Index b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;
    if (size_[a] < size_[b]) std::swap(a, b);
    parent_[b] = a;
    size_[a] += size_[b];
    return true;
  }

 private:
  std::vector<Index> parent_;
  std::vector<std::size_t> size_;
};

}  // namespace sssp::detail
