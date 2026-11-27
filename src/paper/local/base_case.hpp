#pragma once
// Base case of Algorithm 4 (LocalHierarchy, l.2-5) as fixed by Sec 2 (l.286-290):
//   l.3  keep one native entry state per vertex, with its physical image, and both
//        orientations of every original nonloop edge as state arcs;
//   l.4  find an MST after retaining a lightest edge per vertex pair, and merge its
//        equal-weight edges simultaneously to form the hierarchy, using the MST edges as
//        physical connectors (links);
//   l.5  return the states, arcs, maps and hierarchy.
// A node created by the edges of weight w has height w and the components of the MST edges
// of weight < w as children; its links are its weight-w MST edges, so c0 = 1. Lemma 3.5
// (b_T(x, y) <= w(xy)) excludes arcs lighter than h(C) between distinct children, so there are
// no light arcs at all and H3 holds trivially.
//
// Mode A scope: this file builds only the base case. The size test k <= 2 L_k^4 of Alg 4 l.2
// belongs to the Alg 4 entry point, which must also push the FallbackRecord required by the
// project's fallback policy when it takes this branch; this function only bumps counters.

#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// Output of Algorithm 4 (Lemma 5.10): the ordered hierarchy (states, arcs, projection, forest,
/// c0), the native entry state of every vertex, and for every state arc the id of the original
/// edge it came from (used by later stages to recover physical paths).
template <class W>
struct LocalHierarchyResult {
  OrderedHierarchy<W> hierarchy;
  std::vector<StateId> native_entry;  ///< vertex -> native entry state
  std::vector<EdgeId> arc_edge;       ///< state arc -> original edge id
};

/// Alg 4 l.3-5 on `g`. Alg 4 assumes a connected input; on a disconnected one this returns one
/// tree per connected component (several roots, no inter-root arcs, so H4 holds) [supplement].
/// Throws std::invalid_argument if an edge (including a self-loop) has weight 0 or, for
/// integer weights, the infinity sentinel: Alg 4 requires positive weights and zero-weight
/// components are contracted by the caller (Sec 5.1 l.962, Alg 8 l.2).
/// Counters (if `diag` is non-null): local_base_case_calls, local_base_case_states,
/// local_base_case_merge_nodes.
template <class W>
LocalHierarchyResult<W> local_hierarchy_base_case(const CsrGraph<W>& g,
                                                  Diagnostics* diag = nullptr);

}  // namespace sssp::paper
