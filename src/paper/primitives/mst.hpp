#pragma once
// Minimum spanning trees of Sec 3.2 (Mode A, paper-faithful):
//   - the edge order (length, input index)               Sec 3.2, first paragraph
//   - quotient graphs whose edges keep their input index  Sec 3.2, first paragraph
//   - GrowthPass                                          Algorithm 1, Lemma 3.4
//   - the three-pass MST and its bottleneck property      Lemma 3.5
//   - the lightest edge per vertex pair                   Algorithm 4 line 4
//   - a Kruskal reference (sort + union-find)             test oracle, not from the paper
//
// Lengths are only copied and compared (comparison model, Sec 2); no length is added.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// Total order on edges: (length, input index). With unique indices every graph has exactly one
/// minimum spanning forest under this order (Sec 3.2, proof of Lemma 3.4).
template <class W>
struct MstEdgeKey {
  W weight;
  EdgeId id;

  friend bool operator<(const MstEdgeKey& a, const MstEdgeKey& b) noexcept {
    if (a.weight < b.weight) return true;
    if (b.weight < a.weight) return false;
    return a.id < b.id;
  }
};

/// Growth limit B = +infinity of Algorithm 1 (the third pass of Lemma 3.5).
inline constexpr std::size_t kUnboundedGrowth = std::numeric_limits<std::size_t>::max();

/// An edge between two quotient vertices; `id` is the original input index.
template <class W>
struct QuotientEdge {
  Vertex a;
  Vertex b;
  W weight;
  EdgeId id;
};

/// Quotient vertices are connected sets of original vertices (Sec 3.2); only the edges between
/// them are stored. Parallel edges are kept; loops may be present on input to growth_pass and
/// are never produced by it.
template <class W>
struct QuotientGraph {
  std::size_t num_vertices = 0;
  std::vector<QuotientEdge<W>> edges;
};

/// The quotient in which every vertex is its own set: loops dropped, indices kept.
/// Throws std::out_of_range / std::invalid_argument like sssp::validate_edges.
template <class W>
QuotientGraph<W> singleton_quotient(std::size_t n, std::span<const Edge<W>> edges);

/// Per-pass counters (Lemma 3.4); used by tests of the vertex bound and by Diagnostics.
struct GrowthPassStats {
  std::size_t limit = 0;         ///< B (kUnboundedGrowth for +infinity)
  std::size_t vertices_in = 0;   ///< k'
  std::size_t edges_in = 0;      ///< l' (non-loop quotient edges)
  std::size_t vertices_out = 0;  ///< quotient vertices after contracting F
  std::size_t edges_out = 0;
  std::size_t growths = 0;                ///< Alg 1 l.4 iterations
  std::size_t stopped_empty = 0;          ///< Alg 1 l.9: empty heap
  std::size_t stopped_large = 0;          ///< Alg 1 l.9: heap size > B
  std::size_t stopped_attached = 0;       ///< Alg 1 l.12: reached an earlier growth
  std::size_t extractions = 0;            ///< Alg 1 l.11
  std::size_t max_extract_heap_size = 0;  ///< largest heap size at an extraction (<= min(k', B))
  std::uint64_t heap_inserts = 0;
  std::uint64_t heap_decreases = 0;
  /// Edge-key comparisons (heap plus boundary-improvement tests); 0 if SSSP_DIAGNOSTICS=0.
  std::uint64_t comparisons = 0;
};

template <class W>
struct GrowthPassResult {
  std::vector<EdgeId> selected;      ///< F, original indices in selection order
  std::vector<Vertex> component_of;  ///< input quotient vertex -> output quotient vertex
  QuotientGraph<W> quotient;         ///< loop-free quotient by the components of F
  GrowthPassStats stats;
};

/// Algorithm 1 (GrowthPass) with limit B >= 1 (kUnboundedGrowth = +infinity).
/// Throws std::invalid_argument for B = 0 and std::out_of_range for bad endpoints.
/// The input need not be connected (each component then ends with one empty-heap growth).
template <class W>
GrowthPassResult<W> growth_pass(const QuotientGraph<W>& h, std::size_t limit);

/// Finite growth limits of Lemma 3.5 for an input with k vertices:
/// B1 = max(2, ceil(log2 L_k)), B2 = L_k, with L_k = ceil(log2(k + 2)) (Sec 2), evaluated in
/// integer arithmetic without overflow. The third pass uses B = +infinity.
struct GrowthLimits {
  std::size_t b1;
  std::size_t b2;
};
GrowthLimits three_pass_limits(std::size_t k) noexcept;

struct SpanningForestResult {
  std::vector<EdgeId> edges;            ///< selected input indices, increasing
  std::vector<GrowthPassStats> passes;  ///< passes actually run (at most three)
};

/// Lemma 3.5: minimum spanning tree of a connected multigraph with k >= 1 vertices by three
/// GrowthPass rounds with limits B1, B2 (from the original k) and +infinity. The result is the
/// unique MST under the (length, index) order. Throws std::invalid_argument if k = 0 or the
/// graph is disconnected (Lemma 3.5 assumes connectivity; see spanning_forest_three_pass).
/// If `diag` is given, adds per-pass counters "mst.pass<i>.{limit,vertices_in,edges_in,
/// vertices_out,edges_out,growths,extractions}" (limit only for finite B), "mst.passes", and
/// the heap work to `diag->ops`.
template <class W>
SpanningForestResult mst_three_pass(std::size_t n, std::span<const Edge<W>> edges,
                                    Diagnostics* diag = nullptr);
template <class W>
SpanningForestResult mst_three_pass(const CsrGraph<W>& g, Diagnostics* diag = nullptr);

/// The same three passes on a possibly disconnected graph: returns the minimum spanning forest.
/// Not covered by Lemma 3.5 (each component ends its growths on an empty heap, so the vertex
/// bound of Lemma 3.4 becomes c + 2l'/B for c components).
template <class W>
SpanningForestResult spanning_forest_three_pass(std::size_t n, std::span<const Edge<W>> edges,
                                                Diagnostics* diag = nullptr);
template <class W>
SpanningForestResult spanning_forest_three_pass(const CsrGraph<W>& g, Diagnostics* diag = nullptr);

/// For every unordered pair {u, v} with u != v, the index of the least edge joining it under the
/// (length, index) order; loops are dropped. Indices are returned in increasing order.
/// Two stable counting-sort passes (Lemma 3.1), O(n + m).
template <class W>
std::vector<EdgeId> lightest_edge_per_pair(std::size_t n, std::span<const Edge<W>> edges);

/// Reference minimum spanning forest: sort by (length, index), then union-find (Kruskal).
/// Returns the forest's edge indices in increasing order. O(m log m).
template <class W>
std::vector<EdgeId> mst_kruskal_ref(std::size_t n, std::span<const Edge<W>> edges);

}  // namespace sssp::paper
