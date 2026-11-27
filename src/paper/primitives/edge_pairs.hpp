#pragma once
// The lightest edge per vertex pair (Mode A, paper-faithful):
//   Alg 4 l.4 (local base case), Alg 5 l.4 (NameRoutes), Alg 8 l.5 (outer base case), and the
//   size-threshold scan of Sec 2 l.286-288 all "retain a lightest edge per vertex pair".
// Grouping by the unordered pair uses two stable counting-sort passes (Lemma 3.1), O(m + n).
// Lengths are only compared (comparison model, Sec 2): one comparison per further edge of a pair
// group, counted into an optional OpCounters sink (SSSP_COUNT, compiled out when
// SSSP_DIAGNOSTICS=0).
//
// [Supplement] Ties: the paper allows any lightest edge; we keep the least input index, i.e.
// the survivor is the least edge under the (length, index) order of Sec 3.2. This makes the
// result unique and keeps the MST under that order unchanged.

#include <cstddef>
#include <span>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// Survivors of the lightest-edge-per-pair scan.
struct LightestEdgePairs {
  /// One edge per unordered pair {u, v} with u != v that has an edge: the (length, index)-least
  /// one. Input indices in increasing order.
  std::vector<EdgeId> kept;
  /// For every input edge, the kept edge joining the same unordered pair (kept edges point to
  /// themselves); kNoEdge for loops, which are dropped.
  std::vector<EdgeId> representative;
};

/// Throws like sssp::validate_edges (std::out_of_range for bad endpoints, std::invalid_argument
/// for invalid weights, std::length_error for too many vertices or edges). If ops is given, the
/// length comparisons (non-loop edges minus kept edges) are added to ops->comparisons.
template <class W>
LightestEdgePairs lightest_edge_pairs(std::size_t n, std::span<const Edge<W>> edges,
                                      OpCounters* ops = nullptr);

/// Only the kept indices of lightest_edge_pairs (increasing).
template <class W>
std::vector<EdgeId> lightest_edge_per_pair(std::size_t n, std::span<const Edge<W>> edges,
                                           OpCounters* ops = nullptr);

}  // namespace sssp::paper
