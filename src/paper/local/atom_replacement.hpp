#pragma once
// Replacing connected subgraphs by local hierarchies (Sec 5.5, txt l.1302-1364): Definition 5.11
// (replaceable atom), Definition 5.12 (partition hierarchy), the replacement procedure
// (l.1316-1319) and Lemmas 5.13-5.15.
//
// cap_local_hierarchy (Lemma 5.13) runs Algorithm 4 on the atom's internal graph H once and caps
// the result at the cutoff c:
//   * it retains the maximal subtrees of height strictly below c, including the leaves not
//     covered by any such subtree (l.1316-1317); the removed nodes have height >= c;
//   * if several subtrees remain, it gives them a common root of height c (l.1317). Its links
//     follow the proof of Lemma 5.13, item 2 (l.1331-1334): zero-length certificates between
//     copies of the same vertex of H, then one certificate per edge of the given spanning tree
//     between the native copies of its endpoints, taken greedily as a spanning tree of the
//     quotient on the retained subtrees (union-find, Lemma 3.6). Every link is < 2c because
//     every tree edge is (Def 5.11); the root's links are certificate-only (kNoArc, H6), and the
//     result declares c0 = max(c0 of Alg 4, 2).
//   States, arcs, projections (onto the vertices of H) and native entries are those of Alg 4.
//
// substitute_atoms (Lemma 5.15) applies cap_local_hierarchy to every atom of a partition
// hierarchy (one Alg 4 call per atom, Lemma 5.15 proof) and rebuilds the ordered hierarchy on
// the underlying state graph:
//   * every copy keeps the physical image of its old state (l.1319, Def 5.12 l.1312);
//   * an old arc entering an atom, or a seed at an atom state, goes to the endpoint's native
//     inward copy (l.1318; redirect_offers for seeds);
//   * every copy of an old tail gets one copy of each of the tail's old outside arcs (l.1318-
//     1319, Lemma 5.14): the arc keeps its length and its head (the head's entry);
//   * the forest of the partition hierarchy is kept node for node; an atom leaf is replaced by
//     the capped root, a singleton leaf by the leaf of its state; links keep their lengths and
//     name the entries of their endpoint states, and a backing arc becomes its copy at the
//     tail's entry (which joins the same two entries).
//
// Interpretations and supplements (docs/research/open-questions.md, proposed numbers):
//   [Interpretation, proposed Q-53] D (PartitionHierarchy::graph) holds the outside arcs only:
//     an arc with both endpoints in one atom is internal and belongs to the atom's internal
//     graph H (Sec 7.4 l.1816-1818 collects "every internal edge" of an atom separately from the
//     interatom exits). substitute_atoms rejects such arcs instead of guessing a meaning.
//   [Interpretation] The underlying state graph of Def 5.12 is D plus both orientations of every
//     internal edge; Lemma 5.14's witnesses refer to it.
//   [Supplement] The partition forest is a vector of HierarchyNode plus a parallel array
//     leaf_atom (the design rule of ordered_hierarchy.hpp: extra per-node data lives in parallel
//     arrays). The checks below (parts, atom leaves, Def 5.11 cutoff and tree, Def 5.12 cutoff
//     <= parent height, forest structure) are linear and throw std::invalid_argument; Def 4.1 on
//     the parts is the caller's obligation and is checked on the output by
//     validate_ordered_hierarchy.
//   [Supplement] A one-state atom (empty H) is accepted: Alg 4 returns a single leaf, which is
//     its own capped root.
//
// Lengths are only compared, copied and added (Sec 2): 2c is formed by one saturating addition.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/base_case.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

inline constexpr std::uint32_t kNoAtom = std::numeric_limits<std::uint32_t>::max();

/// Definition 5.11 (l.1303-1307). Local vertex i of the internal graph H is states[i]. Every
/// internal edge has a positive length, is usable in both directions and projects to a physical
/// walk of that length between the images of its endpoints (caller's obligation). The given
/// spanning tree of H has every edge shorter than 2 * cutoff, and cutoff > 0.
template <class W>
struct ReplaceableAtom {
  std::vector<StateId> states;        ///< underlying states, distinct
  std::vector<Edge<W>> edges;         ///< internal graph H on 0 .. states.size()-1
  std::vector<EdgeId> spanning_tree;  ///< ids into edges
  W cutoff{};                         ///< c
};

/// Definition 5.12 (l.1308-1315). `graph` is the directed state graph of the OUTSIDE arcs (no
/// arc has both endpoints in one atom), `projection` the physical image of every state. `nodes`
/// is the rooted forest on the parts: a singleton leaf carries its state (leaf_atom = kNoAtom);
/// an atom leaf has no children, state = kNoState and leaf_atom = its index into `atoms`.
/// Internal nodes carry heights and links naming actual states (Def 4.1 on the parts, connection
/// constant c0). An atom with a parent has cutoff <= the parent's height.
template <class W>
struct PartitionHierarchy {
  StateGraph<W> graph;
  std::vector<Vertex> projection;
  std::vector<ReplaceableAtom<W>> atoms;
  std::vector<HierarchyNode<W>> nodes;
  std::vector<std::uint32_t> leaf_atom;  ///< per node
  std::uint32_t c0 = 1;
};

/// Lemma 5.13: the capped local hierarchy of one atom.
template <class W>
struct CappedLocalHierarchy {
  /// Alg 4's states, arcs, native entries and arc -> edge map with the capped forest (exactly
  /// one root, every height <= c). Projections are vertices of H.
  LocalHierarchyResult<W> local;
  NodeId root = kNoNode;
  bool added_root = false;        ///< true iff a common root of height c was created
  std::size_t removed_nodes = 0;  ///< Alg 4 nodes of height >= c dropped
};

/// Lemma 5.13 on the atom with internal graph `internal` (vertices 0..k-1), spanning tree
/// `spanning_tree` (edge ids of `internal`) and cutoff `cutoff`; Alg 4 runs with `p` and pushes
/// its own FallbackRecord when it takes the base case. Throws std::invalid_argument if k = 0,
/// the cutoff is not a positive finite length, or `spanning_tree` is not a spanning tree of
/// `internal` with every edge < 2c (Def 5.11); plus whatever local_hierarchy throws. Throws
/// std::logic_error if Alg 4's result has a height below c above a height of at least c (H2).
/// Counters: atom_cap_calls, atom_cap_added_roots, atom_cap_removed_nodes.
template <class W>
CappedLocalHierarchy<W> cap_local_hierarchy(const CsrGraph<W>& internal,
                                            std::span<const EdgeId> spanning_tree, W cutoff,
                                            Diagnostics* diag = nullptr,
                                            const LocalHierarchyParams& p = kPaperLocalParams);

/// Lemma 5.15: the ordered hierarchy after substituting every atom.
template <class W>
struct AtomSubstitution {
  OrderedHierarchy<W> hierarchy;       ///< on the new states; c0 = max of the inputs' c0
  std::vector<StateId> old_state;      ///< new state -> underlying state
  std::vector<StateId> entry;          ///< underlying state -> its native (inward) copy
  std::vector<ArcId> old_arc;          ///< new arc -> copied arc of D; kNoArc for local arcs
  std::vector<NodeId> atom_root;       ///< per atom: its capped root in `hierarchy`
  std::vector<NodeId> partition_node;  ///< per node of the partition forest: its new node
};

/// Lemma 5.15 on `ph` (Alg 4 parameters `p` for every atom). Throws std::invalid_argument if the
/// parts do not partition the states (each state in exactly one singleton leaf or atom, each
/// atom in exactly one leaf), if an arc of D has both endpoints in one atom, if an atom with a
/// parent has a cutoff above the parent's height, or if the forest is malformed; and
/// std::out_of_range for a state or node id out of range; plus what cap_local_hierarchy throws.
/// Counters: atom_replace_calls, atom_replace_atoms, atom_replace_states,
/// atom_replace_exit_copies, plus the atom_cap_* and local_hierarchy_* counters.
template <class W>
AtomSubstitution<W> substitute_atoms(const PartitionHierarchy<W>& ph, Diagnostics* diag = nullptr,
                                     const LocalHierarchyParams& p = kPaperLocalParams);

/// Replacement procedure l.1318: every seed moves to the native copy of its state. Throws
/// std::out_of_range for an offer at a state that is not an underlying state.
template <class W>
std::vector<StateOffer<W>> redirect_offers(const AtomSubstitution<W>& s,
                                           std::span<const StateOffer<W>> offers);

}  // namespace sssp::paper
