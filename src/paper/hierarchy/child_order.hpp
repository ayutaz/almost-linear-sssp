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
//   * Charges are accumulated as signed pairs (Sec 2, l.294-296): a forest distance adds
//     (rd(a) + rd(b), 2 rd(l)), an arc adds (len, 0). [Deviation candidate, proposed D-15] Each
//     gap is materialised as ONE real by subtracting the two components of its pair, so that
//     the evaluator can use it as an ordinary length in z + S(i) + g(i, j) (Sec 4.3 l.826).
//     For integer W this is exact. For double it is one rounded subtraction per gap, and the
//     root distances are themselves rounded sums, so g may fall below the certified charge by
//     O(N^2 eps) times the total length (both components are rounded monotonically, hence
//     g >= 0 always). The paper's model would keep g as a pair.
//
// The link forest F (Sec 4.1, l.676-681) is the union of all links of the normalised forest as
// undirected edges between states. Its trees are rooted at their smallest state and searched
// depth-first, neighbours in the order the links appear (node id, then link order); the
// resulting preorder is "the" fixed preorder of l.726. build_link_forest() exposes it.

#include <cstdint>
#include <vector>

#include "paper/hierarchy/normalize.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// The link forest F of a normalised hierarchy, rooted and in preorder.
template <class W>
struct LinkForest {
  std::vector<StateId> parent;          ///< parent state in F, kNoState for a tree root
  std::vector<W> parent_length;         ///< certified length of the link to the parent
  std::vector<W> root_distance;         ///< rd(x): sum of link lengths up to the tree root
  std::vector<std::uint32_t> preorder;  ///< state -> preorder position (over all trees)
  std::vector<StateId> order;           ///< preorder position -> state
};

/// Builds F (one undirected edge per link of nh.forest). Throws std::invalid_argument if the
/// links contain a cycle (then Def 4.1 H5 fails and F is not a forest). Root distances are
/// formed with additions only.
template <class W>
LinkForest<W> build_link_forest(const NormalizedHierarchy<W>& nh, Diagnostics* diag = nullptr);

/// Sec 4.2 on a hierarchy returned by normalize(). Fills nh.rank, nh.rep and nh.gap_next,
/// rewrites every forest.nodes[X].children into rank order and sets nh.children_ordered (see
/// NormalizedHierarchy). Roots get rank = position in nh.root_order, rep = kNoState and
/// gap_next = 0. Calling it again recomputes the same result from the (reordered) children.
///
/// Throws std::invalid_argument if nh.queries is not engaged (not produced by normalize), if
/// the links do not form a forest, or if the light-arc quotient of some node has a cycle (the
/// input violated Def 4.1 H3).
///
/// Counters (if `diag` is non-null): child_order_calls, child_order_nodes (internal nodes),
/// child_order_light_arcs, child_order_dfs_roots, child_order_tree_arcs,
/// child_order_arcs_inspected, child_order_forest_distances. Length operations go to
/// diag->ops.
template <class W>
void order_children(NormalizedHierarchy<W>& nh, Diagnostics* diag = nullptr);

}  // namespace sssp::paper
