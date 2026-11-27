#pragma once
// Normalisation of an ordered hierarchy (Sec 4.1, l.674-723, Lemma 4.3) and the shared result
// type that Sec 4.2 (child order) and Sec 4.3 (Alg 2 Scan, Alg 3 Evaluate) consume.
//
// What normalize() does, step by step (paper locations in brackets):
//   1. mass   S(C) = h(C) + lengths of C's links + S of its internal children, bottom-up; a
//              leaf has mass 0 [l.686-688]. S(C) bounds the forest path between any two states
//              below C, hence the physical diameter of C, and S(C) >= h(C).
//   2. links  an internal child C is joined to its parent P ("normalisation link") iff
//              4 S(C) >= h(P) [l.688-689]. 4 S(C) is formed by two doublings.
//   3. units  each connected set of normalisation links gets a table hmin * 2^k, k = 0, 1, ...,
//              starting at the least height hmin of the set and reaching the set's largest
//              (= topmost) height; nu(C) is h(C) rounded down on that table, so
//              nu(C) in (h(C)/2, h(C)] [l.689-691]. The table is walked downwards from the top
//              index of the set, as in the proof of Lemma 4.3 [l.704-706].
//              [Interpretation, Q-24] A node with no normalisation link is a one-node set whose
//              table starts at h(C): nu(C) = h(C).
//   4. contract a normalisation link between equal units is contracted: the child's children
//              and links move to the parent, and the parent keeps its own (topmost original)
//              mass [l.707-713]. [Interpretation, Q-23] Only normalisation links are contracted;
//              across an unlinked boundary nu(P) > 2 S(C) >= 2 nu(C), so units cannot be equal
//              there and the result is the same.
//   5. arcs    self-arcs are discarded [l.719]; arcs between different roots go to the root
//              quotient, which is topologically sorted [l.720, l.722]; every other arc is
//              located at its LCA L in the contracted forest and kept as a LIGHT ARC of L iff
//              its length is < nu(L) [l.720-721]. LCA and level ancestors use Lemma 3.9
//              (ForestQueries, O(1) per query). Light arcs are grouped by L with one counting
//              sort (Lemma 3.1).
//
// Lemma 4.3 guarantees, and the tests check: nu(P) >= 2 nu(C) for every surviving child C of
// P; nu(P) is a power-of-two multiple of nu(C) or 2 S(C) < nu(P); sum_C S(C)/nu(C) = O(N);
// light arcs still form an acyclic quotient at every surviving node.
//
// Model of computation (Sec 2, l.292-311): lengths are only compared, added and doubled
// (model::LenOps, counted into Diagnostics::ops). Table positions and node ids are integers.
// [Supplement] A doubling that overflows the weight type is treated as "larger than every
// representable height" (saturation): it only occurs when 4 S(C) or the next table entry
// exceeds the largest finite value, where the comparison it feeds is then decided anyway.
// [Deviation D-23, proposed] For the same reason a mass S(C) that does not fit W saturates at
// WeightTraits<W>::infinity() (sat::add, the D-13 helper): S(C) is only an upper bound, the link
// test 4 S(C) >= h(P) is then true as it is for the true value, and Evaluate uses S(C) only in
// saturating sums (evaluate.hpp). Counter normalize_saturated_masses.
//
// [Interpretation, proposed Q-50] The contracted forest is returned as an OrderedHierarchy whose
// node heights are the UNITS nu(C) and whose connection constant is 2 c0. With these heights,
// "light" in Def 4.1 (length < h(C)) coincides with the paper's light arcs of Sec 4.1-4.2
// (length < nu(C)), heights still increase towards the root (units at least double), and every
// link of an original node D now attached to C satisfies len <= c0 h(D) < 2 c0 nu(D) =
// 2 c0 nu(C). Hence validate_ordered_hierarchy() checks the contracted forest unchanged. The
// original heights are kept in `original_height`. (Keeping h(top) as the height would make the
// Def 4.1 light-arc condition fail on contracted nodes: arcs of length in [nu, h) between
// grandchildren may form cycles. The paper only claims acyclicity below the common unit.)

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/primitives/forest_queries.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// A light arc of an internal node X (Sec 4.1 l.721): an arc of D whose endpoints lie below
/// distinct children of X, whose LCA is X, and whose length is < nu(X). Tail and head states
/// and the length are graph.arcs[arc].
struct LightArc {
  ArcId arc;
  NodeId tail_child;  ///< child of X containing the tail state
  NodeId head_child;  ///< child of X containing the head state
};

/// Result of Sec 4.1 (normalize) and, once order_children() has run, of Sec 4.2.
///
/// All per-node vectors are parallel arrays indexed by the NodeId of `forest` (size
/// forest.nodes.size()). Leaves carry 0 in unit / mass / original_height.
template <class W>
struct NormalizedHierarchy {
  // ---------------------------------------------------------------- Sec 4.1 (normalize)

  /// The contracted forest. Same state graph, projection and leaf-per-state structure as the
  /// input; node ids are renumbered (surviving input nodes in increasing input id). For every
  /// internal node: height == unit[node] and links = the links of all input nodes contracted
  /// into it; c0 == 2 * input c0 (see the Q-50 note in the file comment). It satisfies
  /// Def 4.1 whenever the input does.
  OrderedHierarchy<W> forest;

  /// nu(C): the unit of Lemma 4.3. Internal nodes only. Evaluate's local slots of C have this
  /// length (Sec 4.3 l.789-797); parent units are >= 2x child units.
  std::vector<W> unit;
  /// S(C): the certified mass (topmost original mass after contraction), an upper bound on the
  /// physical distance between any two states below C (Lemma 4.3). Evaluate's expiry is
  /// D0 + S(C) (Sec 4.3 l.791) and supplemental offers use S(i) (l.826).
  /// WeightTraits<W>::infinity() if S(C) does not fit W (saturated, D-23).
  std::vector<W> mass;
  /// h of the surviving input node (the topmost member of a contracted group).
  std::vector<W> original_height;
  /// forest node -> the input node that survived as it.
  std::vector<NodeId> original_node;
  /// input node -> forest node (a contracted input node maps to the node that absorbed it).
  std::vector<NodeId> node_of_input;

  /// Light arcs grouped by their LCA: light(X) = light_arcs[light_begin[X] .. light_begin[X+1]).
  /// Arcs at their LCA with length >= nu(LCA) and self-arcs are not listed (Evaluate relaxes
  /// every arc of forest.graph; this list only constrains the child order).
  std::vector<std::size_t> light_begin;
  std::vector<LightArc> light_arcs;
  /// Arcs of D between different roots (the root quotient, Def 4.1 H4).
  std::vector<ArcId> root_arcs;
  /// The roots of `forest` in a topological order of the root quotient (Alg 3 l.4 processes
  /// roots in this order).
  std::vector<NodeId> root_order;

  /// LCA, level ancestor and depth on `forest` (Lemma 3.9). Engaged after normalize(). Alg 2
  /// l.4 finds "the child of X containing a leaf" as level_ancestor(leaf, depth(X) + 1).
  std::optional<ForestQueries> queries;

  // ---------------------------------------------------------------- Sec 4.2 (order_children)

  /// True once order_children() has filled the fields below. Then, for every internal node X,
  /// forest.nodes[X].children lists the children in RANK ORDER (a topological order of the
  /// light-arc quotient of X: every light arc goes from a lower to a higher rank).
  bool children_ordered = false;
  /// Rank of a node among its parent's children (0-based). For roots: position in root_order.
  std::vector<std::uint32_t> rank;
  /// Representative state x_i of a non-root node i (a state below i) [Sec 4.2 l.729].
  /// kNoState for roots.
  std::vector<StateId> rep;
  /// gap_next[i] = g(i, j) for the child j ranked immediately after i under the same parent
  /// (Lemma 4.4 / 4.5): an upper bound, in either direction, on a physical walk between
  /// pi(rep[i]) and pi(rep[j]). 0 for the last child of a node and for roots.
  /// WeightTraits<W>::infinity() if g(i, j) does not fit W (saturated, D-26).
  /// g(j, i) = g(i, j): the gap to the previous sibling is gap_next of that sibling.
  std::vector<W> gap_next;

  [[nodiscard]] std::span<const LightArc> light(NodeId x) const noexcept {
    return {light_arcs.data() + light_begin[x], light_arcs.data() + light_begin[x + 1]};
  }
  [[nodiscard]] std::size_t num_nodes() const noexcept { return forest.nodes.size(); }
  [[nodiscard]] bool is_leaf(NodeId x) const noexcept { return forest.is_leaf(x); }
};

/// Sec 4.1 / Lemma 4.3 on an ordered hierarchy that satisfies Def 4.1 (the caller's obligation;
/// use validate_ordered_hierarchy to check it). O(M + N) apart from the Lemma 3.9 tables.
///
/// Throws std::invalid_argument if c0 == 0, if a height is not positive (H2), if the forest
/// is malformed (some node unreachable from `roots`), or if the root quotient is cyclic (H4).
///
/// Invalid lengths (policy for every length of Def 4.1, review A-4): normalize() is the entry
/// check of Evaluate (Alg 3 l.2), and it raises std::invalid_argument, never
/// std::overflow_error or std::logic_error, for
///   * an arc of D whose length is negative, NaN, +inf (double) or the integer infinity
///     sentinel WeightTraits<W>::infinity(): every arc, also self-arcs (which step 5 then
///     discards) and arcs between roots (Def 4.1 l.646, "nonnegative arc lengths"; infinity is
///     a separate flag, not a length, Sec 2 l.293);
///   * a link whose length is one of these, or exceeds c0 h(C) (H5, l.656). c0 h(C) is formed
///     by the binary expansion of c0 (doublings and additions, O(log c0) per node, so the
///     O(M + N) bound holds for every accepted c0); if it does not fit W it bounds every finite
///     link. It is exact for integers and for double when c0 = 2^k (every producer has
///     c0 <= 2). For double, c0 = 3 * 2^k rounds once to nearest, which never rejects a link
///     <= c0 h(C); for other c0 the roundings add up, so a double link within a few ulp of
///     c0 h(C) may be judged differently from the exact rule and from
///     validate_ordered_hierarchy (kLinkTooLong, c0 - 1 additions). Without this check S(C)
///     would be unbounded and Evaluate's root loop would not terminate in reasonable time
///     (review F-3);
///   * a height that is NaN, +inf or the sentinel (via model::Len, as for H2).
/// -0.0 is the real 0 and is accepted. std::overflow_error stays reserved for distances that do
/// not fit W (include/sssp/paper.hpp). The other conditions of Def 4.1 / 4.2 (H1, H3, H6, link
/// endpoints, offers) are not checked here; validate_ordered_hierarchy checks them all.
///
/// Counters (if `diag` is non-null): normalize_calls, normalize_input_nodes, normalize_nodes,
/// normalize_normalization_links, normalize_link_sets, normalize_table_steps,
/// normalize_contracted_nodes, normalize_self_arcs, normalize_light_arcs,
/// normalize_root_arcs. Length operations are counted in diag->ops.
template <class W>
NormalizedHierarchy<W> normalize(const OrderedHierarchy<W>& h, Diagnostics* diag = nullptr);

}  // namespace sssp::paper
