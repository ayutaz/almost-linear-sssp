#pragma once
// Entry and exit of Algorithm 8 (arXiv 2611.34911v1, Sec 8.3, txt l.2145-2185):
//
// - contract_zero_edges: Alg 8 l.2. Contract every zero-edge component, discard all resulting
//   loops, retain the source component and keep original edge identifiers. The result has only
//   positive edges (Proof of Theorem 1, l.2158-2159), which Sec 5 and the base case require.
//   Parallel edges between two classes are KEPT: l.2 does not deduplicate; "retain a lightest
//   edge per vertex pair" is the separate step l.5 (primitives/edge_pairs).
// - restore_predecessor_tree: Lemma 8.5 as used by Alg 8 l.6 and l.16. From exact distances of
//   the contracted source component, give every original vertex the label of its class (infinity
//   outside the source component or for an infinite label), traverse the finite tight original
//   arcs d(u) + w(uv) = d(v) < inf from s, keep the first-discovery arcs as the tree, and output
//   the sums of original weights down that tree.
//
// Both are single graph traversals, O(m + n) (Sec 8.3, proof of Theorem 1, step 1 and step 3).

#include <cstdint>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::outer {

/// Output of Alg 8 l.2 on (G, s).
template <class W>
struct ZeroContraction {
  /// Contracted multigraph on N = number of zero-edge classes of the source component: positive
  /// weights, no loops, parallel edges kept. Edge i is original edge original_edge[i].
  CsrGraph<W> graph;
  /// Contracted source, the class of s.
  Vertex source = kNoVertex;
  /// Original vertex -> contracted vertex; kNoVertex outside the source component.
  std::vector<Vertex> vertex_map;
  /// Contracted vertex -> original representative (the first member reached by the traversal;
  /// the representative of the source class is s itself).
  std::vector<Vertex> representative;
  /// Contracted edge id -> original edge id, strictly increasing.
  std::vector<EdgeId> original_edge;
  /// Zero-weight spanning tree of every class (the traversal forest of its zero edges), rooted
  /// at the representative: parent vertex and original edge id. kNoVertex / kNoEdge for
  /// representatives and for vertices outside the source component.
  std::vector<Vertex> zero_parent;
  std::vector<EdgeId> zero_parent_edge;
};

/// Alg 8 l.2. Throws std::out_of_range if source is not a vertex of g. If diag is non-null,
/// adds the counters contract_{input,retained}_{vertices,edges}, contract_merged_vertices,
/// contract_unreached_vertices, contract_removed_vertices (= merged + unreached),
/// contract_dropped_loops, contract_unreached_edges, contract_removed_edges (= loops +
/// unreached) and contract_zero_edges (zero-weight edges inside the source component), and
/// counts vertex/edge scans and comparisons into diag->ops.
template <class W>
ZeroContraction<W> contract_zero_edges(const CsrGraph<W>& g, Vertex source,
                                       Diagnostics* diag = nullptr);

/// Lemma 8.5. contracted_dist[x] must be the exact distance of contracted vertex x from
/// c.source in c.graph (WeightTraits<W>::infinity() for unreachable). Throws
/// std::invalid_argument if c does not belong to (g, source), if the length is wrong, if the
/// source label is not 0, or if some finite label is not attained by a tight original arc (the
/// labels are then not exact distances). The result's diag holds the counters
/// restore_tree_vertices and restore_unreachable_vertices and the scans, comparisons and
/// additions performed.
template <class W>
SsspResult<W> restore_predecessor_tree(const CsrGraph<W>& g, Vertex source,
                                       const ZeroContraction<W>& c,
                                       const std::vector<W>& contracted_dist);

extern template struct ZeroContraction<std::uint64_t>;
extern template struct ZeroContraction<double>;
extern template ZeroContraction<std::uint64_t> contract_zero_edges(const CsrGraph<std::uint64_t>&,
                                                                   Vertex, Diagnostics*);
extern template ZeroContraction<double> contract_zero_edges(const CsrGraph<double>&, Vertex,
                                                            Diagnostics*);
extern template SsspResult<std::uint64_t> restore_predecessor_tree(
    const CsrGraph<std::uint64_t>&, Vertex, const ZeroContraction<std::uint64_t>&,
    const std::vector<std::uint64_t>&);
extern template SsspResult<double> restore_predecessor_tree(const CsrGraph<double>&, Vertex,
                                                            const ZeroContraction<double>&,
                                                            const std::vector<double>&);

}  // namespace sssp::paper::outer
