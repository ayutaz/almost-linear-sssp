#pragma once
// Kruskal merge forest with equal weights merged simultaneously (Mode A, paper-faithful):
//   Sec 2 l.288-290  local base case: merge equal-weight MST edges simultaneously; the tree
//                    edges connect the merge children and Lemma 3.5 excludes lighter edges
//                    between them;
//   Alg 4 l.4        the same step inside LocalHierarchy;
//   Sec 5.4 l.1198   Y trees: sort the edges of the local MST, then process equal-weight edges
//                    simultaneously;
//   Sec 3.3          the merge events are processed in the given order by the disjoint sets of
//                    Lemma 3.6 (DisjointSets); ordering them is charged to the caller.
//
// Node heights are edge weights; lengths are only copied and compared (Sec 2).
//
// Lemma 3.5 (b_T(x, y) <= w(xy), and b_T is the minimax path value) gives: for every edge xy
// whose endpoints lie below two different children of a node C, w(xy) >= h(C); the height of
// lca(x, y) is the minimax value between x and y.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// Node of a merge forest: leaves 0..n-1 are the vertices, internal nodes follow.
using MergeNodeId = std::uint32_t;
inline constexpr MergeNodeId kNoMergeNode = std::numeric_limits<MergeNodeId>::max();

template <class W>
struct MergeNode {
  W height{};                         ///< the common weight of the class that created it
  std::vector<MergeNodeId> children;  ///< >= 2, increasing ids
  std::vector<EdgeId> edges;          ///< its forest edges (weight == height), children - 1,
                                      ///< in (weight, id) order; they connect the children
};

template <class W>
struct MergeForest {
  std::size_t num_leaves = 0;
  /// Internal node num_leaves + i is internal[i]. Created class by class in increasing weight;
  /// inside one class in the order of their first forest edge. Children precede parents.
  std::vector<MergeNode<W>> internal;
  std::vector<MergeNodeId> parent;   ///< every node; kNoMergeNode for roots
  std::vector<MergeNodeId> roots;    ///< parentless nodes, increasing (one per component)
  std::vector<EdgeId> forest_edges;  ///< minimum spanning forest under (weight, id), that order

  [[nodiscard]] std::size_t num_nodes() const noexcept { return num_leaves + internal.size(); }
  [[nodiscard]] bool is_leaf(MergeNodeId x) const noexcept { return x < num_leaves; }
};

/// Kruskal on the candidate edges (ids into `edges`; loops and duplicates are harmless) by the
/// (weight, id) order, each equal-weight class merged simultaneously: every component formed in
/// a class that joins at least two components from before the class becomes one node of that
/// height. Input need not be connected (one root per component).
/// Throws like sssp::validate_edges, and std::out_of_range for a candidate id >= edges.size().
/// Counters (if `diag` is non-null): merge_events.{calls,nodes,forest_edges,weight_classes}
/// and the DisjointSets counters dsu.{find_calls,find_steps,unions}.
template <class W>
MergeForest<W> build_merge_events(std::size_t n, std::span<const Edge<W>> edges,
                                  std::span<const EdgeId> candidates, Diagnostics* diag = nullptr);

/// The same with every edge as a candidate.
template <class W>
MergeForest<W> build_merge_events(std::size_t n, std::span<const Edge<W>> edges,
                                  Diagnostics* diag = nullptr);

}  // namespace sssp::paper
