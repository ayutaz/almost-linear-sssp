#pragma once
// Owner decomposition (Sec 5.1, txt l.960-1018; Definition 5.1, Lemmas 5.2 and 5.3).
//
// Construction on a connected positive-weight graph G with N vertices and separator size b:
//   1. T = MST of G (Lemma 3.5, mst_three_pass). D_T = tree path length, b_T = largest weight
//      on a tree path. The reduced graph omits every nontree edge with w(uv) > D_T(u, v)
//      (l.965-966); D_T comes from root distances and the LCA structure of Lemma 3.9, and the
//      test w <= D_T(u, v) is evaluated as w + 2 rd(lca) <= rd(u) + rd(v) (additions only).
//   2. A0: T is processed in postorder; each vertex sums the residual sizes of its children
//      plus one and is marked when the sum reaches b, returning residual 0 (l.968-971). If no
//      mark was made, one pivot is used (l.971).
//   3. A = A0 plus the branch (degree >= 3) vertices of the minimal subtree of T spanning A0
//      (l.972-975). Every component of T - A then has < b vertices and <= 2 pivot neighbours.
//   4. Corridor cuts (l.983-985): for each two-pivot component, a maximum-weight tree edge on
//      the corridor between its two pivots is deleted; two adjacent pivots lose their common
//      edge. The remaining components are the owner cells, one pivot each.
//   5. p(v) = pivot of v's cell, r_v = b_T(v, p(v)) (0 at pivots), R_v = D_T(v, p(v)) <= h r_v
//      with h = b + 1 (l.986-988).
//   6. The deleted corridor edges form a tree T_A on A; sigma(p, q) = largest weight on
//      T_A[p, q] = b_T(p, q), sigma(p, p) = 0 (l.994-996), answered through the sorted merge
//      hierarchy of T_A (Kruskal tree, O(|A| log(|A| + 1)) to sort) and its LCA (Lemma 3.9).
//
// The local construction (Alg 4 l.6) uses b = L_N^4 (local_separator_size); the outer one
// (Sec 7.1) passes its own b (l.981). Mode B passes a small b to exercise the construction on
// small graphs (Q-03); any b other than the paper's is a deviation of the caller.
//
// [Interpretation, Q-31] The paper fixes neither the root of T, nor the pivot used when no mark
// is made, nor the tie-break among maximum corridor edges. Here T is rooted at vertex 0, the
// unmarked case uses the root as its single pivot (and records it in `marked`), and among
// equal maximum corridor edges the one with the smallest input edge id is cut.
// [Interpretation] Corridors are found by rooting T at the smallest pivot and walking upward
// from every other pivot to its first pivot ancestor; the walked vertices are nonpivot
// degree-two vertices of the minimal subtree, so every vertex is walked at most once (the
// paper's "following degree-two paths in that subtree", l.989-990).
//
// Lengths are opaque reals (Sec 2): only comparisons and additions of W are performed. For
// double weights the sums D_T, R_v are rounded (Q-06: correctness checks use integers).
//
// Range of the arithmetic (overflow rule of include/sssp/paper.hpp). Sec 5.1 has no source: its
// sums are tree distances, auxiliary quantities and not tentative path lengths, so none of them
// raises std::overflow_error.
//   * Root distances D_T(root, v) are formed in TreeLength<W>. Integer W: unsigned __int128,
//     exact (a root path has fewer than 2^32 edges of less than 2^64 each), and both sides of
//     the test w + 2 rd(lca) <= rd(u) + rd(v) are exact 128-bit sums (below 2^98), so the
//     reduced graph is decided exactly whatever the scale of the weights. double: rounded IEEE
//     sums; a side that exceeds the range is +infinity. If only the left side overflows the edge
//     is omitted (it is heavier than its finite tree path); if the right side overflows the
//     edge is kept, also when both do. Keeping an edge is always safe for distances (the
//     reduced graph is then a superset of the paper's); the bound of Lemma 5.2 is not claimed
//     for such an edge.
//     [Deviation candidate: the 128-bit intermediate is a numeric representation the paper's
//     real model does not need; the formula and its decisions are unchanged.]
//   * R_v = D_T(v, p(v)) is returned in W. A value that is not below
//     WeightTraits<W>::infinity() (integers: R_v >= 2^64 - 1; double: beyond DBL_MAX) is
//     returned as infinity(), and so is R of every vertex behind it in the cell (the saturating
//     sum of paper/local/saturating.hpp, Q-38). A user of R_v that needs it as a length must
//     treat infinity() as "does not fit" (Sec 7.1 does: it raises on A_e = R_x + w + R_y). The
//     local construction (Sec 5.2-5.4) reads r_v and sigma only, never R_v.
//     [Deviation candidate: saturation of R_v, the treatment of D-13 applied to Sec 5.1.]

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "paper/model/wide_length.hpp"
#include "paper/primitives/forest_queries.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

inline constexpr std::uint32_t kNoPivot = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::uint32_t kNoCorridor = std::numeric_limits<std::uint32_t>::max();

/// Type of a sum of edge weights along a path of T from the root: unsigned __int128 for integer
/// W (exact for every tree), W itself for floating point (see "Range of the arithmetic" above).
template <class W>
using TreeLength = WideLength<W>;

/// b = L_N^4 for the local construction (Sec 5.1 l.967), saturating at 2^64 - 1.
std::size_t local_separator_size(std::size_t n);

/// One edge of the pivot tree T_A: the corridor between the pivots with indices `lower` and
/// `upper` (indices into OwnerDecomposition::pivots) was cut at input edge `edge` of weight
/// `weight`. `upper` is the first pivot above `lower` when T is rooted at pivots[0].
template <class W>
struct CorridorEdge {
  std::uint32_t lower;
  std::uint32_t upper;
  EdgeId edge;
  W weight;
};

/// Definition 5.1 together with the derived quantities of Sec 5.1. Vertex-indexed vectors have
/// size N, edge-indexed vectors have size m (input edge ids).
template <class W>
struct OwnerDecomposition {
  std::size_t separator = 0;  ///< b
  std::size_t hop_bound = 0;  ///< h = b + 1 (bound on the edges of T[v, p(v)])
  Vertex root = 0;            ///< root of T for the postorder marking

  std::vector<EdgeId> tree_edges;        ///< T, increasing input ids
  std::vector<Vertex> tree_parent;       ///< T rooted at `root` (kNoVertex at the root)
  std::vector<EdgeId> tree_parent_edge;  ///< edge to the parent (kNoEdge at the root)
  /// D_T(root, v). Exact for integer W (it may exceed the range of W); a rounded sum for
  /// double, +infinity beyond the range.
  std::vector<TreeLength<W>> root_distance;
  std::vector<std::uint8_t> reduced;  ///< per edge: 1 iff the edge is in the reduced graph

  std::vector<std::uint8_t> marked;        ///< per vertex: 1 iff in A0
  std::vector<Vertex> pivots;              ///< A, increasing vertex order
  std::vector<std::uint32_t> pivot_index;  ///< vertex -> index into pivots, or kNoPivot
  std::vector<CorridorEdge<W>> corridors;  ///< T_A (|A| - 1 edges)
  std::vector<std::uint8_t> cut;           ///< per edge: 1 iff a deleted corridor edge

  std::vector<std::uint32_t> owner;  ///< vertex -> pivot index of p(v)
  std::vector<W> radius;             ///< r_v = b_T(v, p(v))
  /// R_v = D_T(v, p(v)); WeightTraits<W>::infinity() if the sum does not fit W.
  std::vector<W> owner_distance;
  std::vector<std::uint32_t> owner_hops;  ///< number of edges of T[v, p(v)]

  /// Sorted merge hierarchy of T_A: nodes 0..|A|-1 are the pivots, node |A| + i is the merge
  /// made by corridors[merge_corridor[i]] in increasing (weight, edge id) order.
  std::vector<Vertex> merge_parent;
  std::vector<std::uint32_t> merge_corridor;
  std::optional<ForestQueries> merge_queries;

  [[nodiscard]] std::size_t num_vertices() const noexcept { return owner.size(); }
  [[nodiscard]] bool is_pivot(Vertex v) const { return pivot_index[v] != kNoPivot; }

  /// Index into `corridors` of the largest (weight, edge id) corridor edge on T_A[p, q]; the
  /// "named corridor maximum" of Sec 5.2 (l.1022). kNoCorridor for p = q.
  [[nodiscard]] std::uint32_t corridor_max(std::uint32_t p, std::uint32_t q) const {
    if (p == q) return kNoCorridor;
    const std::optional<Vertex> a = merge_queries ? merge_queries->lca(p, q) : std::nullopt;
    if (!a) throw std::out_of_range("OwnerDecomposition::corridor_max: no merge hierarchy");
    return merge_corridor[*a - pivots.size()];
  }

  /// sigma(p, q) = b_T(p, q) for pivot indices p, q; sigma(p, p) = 0 (l.995-996).
  [[nodiscard]] W sigma(std::uint32_t p, std::uint32_t q) const {
    const std::uint32_t c = corridor_max(p, q);
    return c == kNoCorridor ? W{0} : corridors[c].weight;
  }
};

/// Sec 5.1 on a connected graph with N >= 1 vertices and positive finite weights, separator
/// size b >= 1. Throws std::invalid_argument if N = 0, b = 0, the graph is disconnected, or an
/// edge (including a loop) has weight 0 or (integers) the infinity sentinel. Never throws
/// std::overflow_error: tree distances that leave the range of W are exact (root_distance) or
/// saturated (owner_distance), see "Range of the arithmetic" above.
/// Counters (if `diag` is non-null): owner_decomposition_calls, owner_decomposition_marks,
/// owner_decomposition_pivots, owner_decomposition_branch_pivots,
/// owner_decomposition_unmarked, owner_decomposition_dropped_edges; plus those of the MST.
template <class W>
OwnerDecomposition<W> owner_decomposition(const CsrGraph<W>& g, std::size_t b,
                                          Diagnostics* diag = nullptr);

}  // namespace sssp::paper
