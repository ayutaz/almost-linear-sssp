#pragma once
// Ordering the children of a normalised ordered hierarchy (Sec 4.2, l.725-779, Lemma 4.4 and
// Lemma 4.5).
//
// At every internal node X (l.726-730):
//   * each child i gets a provisional representative: the state below i that comes first in
//     the F preorder ([Interpretation, Q-25]; the paper allows any state);
//   * a directed DFS runs on the light child DAG of X (light arcs from normalize(), length
//     < nu(X)), launching at unvisited children in the F preorder of their provisional
//     representatives. A DFS root i is entered at its provisional representative; a DFS tree
//     arc enters its head child at the arc's head state. At child i, the outgoing light arcs
//     are inspected in the CYCLIC F preorder of their tail states, starting at the entry
//     state x_i; ties at one tail keep arc order (l.729-730 allows any);
//   * the final representative rep(i) is the entry state x_i; the rank order is the reverse
//     finishing order, which is topological (Lemma 4.5).
//
// Gaps (proof of Lemma 4.4, l.737-751): the traversal is charged the certified length of every
// F edge it moves along and the length of every DFS tree arc it crosses (each tree arc twice,
// forward and back). Within a child it moves between the selected exit tails along F, starting
// and ending at its entry; between successive DFS roots it moves along F between their
// provisional representatives. The gap g(i, j) of two children consecutive in rank order is the
// total charge between their marks (finishing events). Every charged step certifies a physical
// walk of at most its charge (links by Def 4.1 H5, arcs by H1), so g(i, j) bounds a physical
// walk between pi(rep i) and pi(rep j) in either direction, and
//     sum of the gaps at X <= 2 sum_i S(i) + 2 d nu(X) + 2 S(X)              (Lemma 4.4).
//
// Implementation notes (Lemma 4.5 proof, l.762-777):
//   * F distances are rd(a) + rd(b) - 2 rd(lca_F(a, b)) with root distances rd and the
//     Lemma 3.9 LCA (ForestQueries); no physical path is walked explicitly.
//   * DFS roots: one stable radix sort of the provisional representatives by (parent, F
//     preorder), Lemma 3.1.
//   * Cyclic start in O(1): tails of light arcs, entry positions (heads of light arcs and
//     provisional representatives) are radix-sorted by (child, F preorder, kind); a backward
//     sweep stores, for every entry position, the next tail at or after it, wrapping to the
//     child's first tail (l.773-777). A child with no outgoing light arcs finishes at once.
//     This is the paper's table, not a binary search (no D-4).
//   * [Deviation candidate, proposed D-15] The charge between two marks is kept as ONE
//     nonnegative value instead of the signed pair of Sec 2 (l.294-296): a forest distance adds
//     (rd(a) - rd(l)) + (rd(b) - rd(l)) with l = lca_F(a, b), two subtractions of root distances
//     along one root path (rd is nondecreasing towards the leaves, so both differences are
//     nonnegative; the same form as D-22), an arc adds its length. The gap is that value, so
//     the evaluator can use it as an ordinary length in z + S(i) + g(i, j) (Sec 4.3 l.826).
//
// Range of the arithmetic (overflow rule of include/sssp/paper.hpp). Root distances and charges
// are auxiliary quantities, not path lengths, so they never raise:
//   * integer W: root distances and charges are formed in ForestLength<W> = unsigned __int128
//     and are exact (fewer than 2^32 links of less than 2^64 each on a root path; a charge is
//     saturated at 2^128 - 1, which no input that fits in memory reaches). A gap below the
//     sentinel WeightTraits<W>::infinity() is returned exactly; a larger one is returned as
//     infinity(), which is still an upper bound on the walk (the saturation of D-13 applied to
//     g). [Deviation candidate: the 128-bit intermediate is a numeric representation the
//     paper's real model does not need.]
//   * double: root distances are rounded sums (ForestLength<double> = double). g is exact
//     whenever the root distances it uses are exact (e.g. lengths that are multiples of a common
//     2^-p with sums below 2^(53 - p), at any power-of-two scale); otherwise each root distance
//     carries its rounding and g may differ from the certified charge by up to the spacing of
//     doubles at the root distances involved (g >= 0 always). No comparison of Sec 4.2 depends
//     on a length, so rounding never changes the order, the representatives or the DFS. A
//     root distance beyond the range is +infinity; a forest distance that uses one, and a
//     charge whose rounded sum is not finite, is +infinity (never NaN), also when the true gap
//     would fit.
//
// The link forest F (Sec 4.1, l.676-681) is the union of all links of the normalised forest as
// undirected edges between states. Its trees are rooted at their smallest state and searched
// depth-first, neighbours in the order the links appear (node id, then link order); the
// resulting preorder is "the" fixed preorder of l.726. build_link_forest() exposes it.

#include <cstdint>
#include <type_traits>
#include <vector>

#include "paper/hierarchy/normalize.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

namespace detail {
__extension__ typedef unsigned __int128 UInt128;  // not ISO C++; exact sums of 64-bit lengths
}  // namespace detail

/// Type of a sum of link lengths along a path of F: unsigned __int128 for integer W (exact for
/// every forest), W itself for floating point (see "Range of the arithmetic" above).
template <class W>
using ForestLength = std::conditional_t<std::is_integral_v<W>, detail::UInt128, W>;

/// The link forest F of a normalised hierarchy, rooted and in preorder.
template <class W>
struct LinkForest {
  std::vector<StateId> parent;   ///< parent state in F, kNoState for a tree root
  std::vector<W> parent_length;  ///< certified length of the link to the parent
  /// rd(x): sum of link lengths up to the tree root. Exact for integer W (it may exceed the
  /// range of W); a rounded sum for double, +infinity beyond the range.
  std::vector<ForestLength<W>> root_distance;
  std::vector<std::uint32_t> preorder;  ///< state -> preorder position (over all trees)
  std::vector<StateId> order;           ///< preorder position -> state
};

/// Builds F (one undirected edge per link of nh.forest). Throws std::invalid_argument if the
/// links contain a cycle (then Def 4.1 H5 fails and F is not a forest). Root distances are
/// formed with additions only and never raise on the range of W.
template <class W>
LinkForest<W> build_link_forest(const NormalizedHierarchy<W>& nh, Diagnostics* diag = nullptr);

/// Sec 4.2 on a hierarchy returned by normalize(). Fills nh.rank, nh.rep and nh.gap_next,
/// rewrites every forest.nodes[X].children into rank order and sets nh.children_ordered (see
/// NormalizedHierarchy). Roots get rank = position in nh.root_order, rep = kNoState and
/// gap_next = 0. Calling it again recomputes the same result from the (reordered) children.
///
/// Throws std::invalid_argument if nh.queries is not engaged (not produced by normalize), if
/// the links do not form a forest, or if the light-arc quotient of some node has a cycle (the
/// input violated Def 4.1 H3). Never throws std::overflow_error: a gap that does not fit W is
/// WeightTraits<W>::infinity() in nh.gap_next (see "Range of the arithmetic" above).
///
/// Counters (if `diag` is non-null): child_order_calls, child_order_nodes (internal nodes),
/// child_order_light_arcs, child_order_dfs_roots, child_order_tree_arcs,
/// child_order_arcs_inspected, child_order_forest_distances. Length operations go to
/// diag->ops.
template <class W>
void order_children(NormalizedHierarchy<W>& nh, Diagnostics* diag = nullptr);

}  // namespace sssp::paper
