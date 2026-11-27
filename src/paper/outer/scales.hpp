#pragma once
// Numerical scales and atoms of the outer construction: arXiv 2611.34911v1, Sec 8.2 (txt
// l.2076-2104) and Lemma 8.4 (l.2105-2142), Mode A. This is the part of Alg 8 l.13 between the
// outer state graph of Sec 8.1 (outer_states.hpp) and the atom substitution of Lemma 5.15
// (local/atom_replacement.hpp). The numeric decisions are those of scales_numeric.hpp; this
// file applies them to the blocks, joins and arcs of the outer states.
//
// Stages (docs/research/04_outer_scales_pipeline.md Sec 2):
//   S0  New names (l.2077-2079): the name of a macro edge becomes the least physical weight of
//       any macro edge carrying the same route name (outer_name_weights). [Q-65] "Physical
//       weight" is the input weight w_e (Lemma 7.1 bounds w), the minimum is over macro edges.
//   S1  F (l.2080-2081) is chosen once per graph from integer parameters and passed in as an
//       exact exponent. [Q-71] It must come from the exact h = b + 1: for the paper's
//       b = L_N^64, which no word holds, paper_outer_scale_factor forms b in multiword integers.
//   S2  Block scales (Sec 8.1 l.1999, Sec 8.2 l.2083-2086): a singleton block (every pivot copy,
//       every special state, any one-state ordinary block) keeps only its finite theta join.
//       For a nonsingleton positive ordinary block U: s_U = min r_v over its states, with the
//       original radius r_v = b_T(v, p(v)) [Q-35]; z_U minimises theta over U and a_U is the new
//       name of the macro attaining B0(z_U). [Supplement] Ties go to the least id.
//   S3  Alignment (l.2086, l.2130-2131): U is aligned iff F^2 s_U >= a_U. The tests are grouped
//       by name and decided by median pruning (Sec 6.1), one threshold test per median.
//   S4  Supports and the global grid (l.2081-2091, Lemma 8.4 item 2): one key per route name
//       that names a finite join and one per corridor of T_A [Q-65: a corridor join uses the
//       support of its own weight]; build_outer_global_grid.
//   S5  Cutoffs (l.2091-2093, l.2100-2104): an aligned block takes the floor of s_U on its
//       name's global grid. Every other nonsingleton block [Q-66: unaligned or without a finite
//       attachment] takes it on the separate grid of the vertex v attaining s_U: [Q-67] the
//       graph-piece grid of Sec 7.3 (Calendars), read in input units as origin * 2^j and
//       extended once per piece by kSeparateGridExtraOctaves doublings at the top; the grid
//       pointer of v (Calendars::radius_index) gives the floor of r_v = s_U directly.
//       Both give s_U / 2 < t_U <= s_U [Q-36].
//   S6  Atoms (Sec 8.1 l.2000-2001, Lemma 8.4 item 3, Def 5.11): the components of the block
//       MST edges lighter than t_U. A component with at least two states is a replaceable atom
//       with cutoff t_U; [Q-74] its internal graph holds every within-window edge with both
//       endpoints in it (heavier chords included), its spanning tree the MST edges below t_U;
//       one-state components are singleton parts.
//   S7  Classes (l.2091-2094): a block MST edge below t_U is internal to an atom; one at or
//       above t_U is classified on the global grid of a_U (aligned) or on the separate grid of
//       its block; a state-to-owner join at theta(z) on the grid of the name of z; a corridor
//       join on the grid of its own key. Every class is [g, 2g) (the upper end may saturate).
//   S8  Partition hierarchy (l.2094-2099, Def 5.12): build_outer_merge_plan (parts, outside
//       arcs, classes in processing order) and build_outer_partition (simultaneous merges,
//       links, child and root order through the shared engine local/assembly_core.hpp).
//   S9  Substitution (l.2099, Lemma 5.15) is the caller's one call
//       substitute_atoms(build_outer_partition(..).hierarchy, diag, local_params); its result
//       is the ordered hierarchy of Sec 4 with c0 <= 2.
//
// Relations the construction relies on are checked where they are used, and a violation throws
// std::logic_error: a nonsingleton block is positive and ordinary; theta >= s_U on it (Lemma
// 8.2); a finite theta has a name with a macro edge; every classified value lies in the index
// range of its key; a significant weight of a separate-grid block lies below the top of the
// extended table (w <= 4 K^2 R_v); an unaligned block has every MST edge below a_U / 2 and
// completes below its first global class (2 p_top <= g, l.2093-2094). check_outer_scales
// re-derives all invariants of note 04 Sec 2.1-2.8 on an output.
//
// Model of computation (Sec 2): lengths are only compared, copied, added and doubled; names,
// keys, grid indices and class indices are integers.
// Numerics (Q-06, D-13): products by powers of F saturate at WeightTraits<W>::infinity() and
// are only compared with finite lengths, which keeps every decision exact (scales_numeric.hpp).
//
// Fallback policy: Sec 8.2 never substitutes a search; it pushes no FallbackRecord. The atoms
// record their Alg 4 base cases when they are substituted (local_hierarchy).

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/atom_replacement.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/outer/outer_states.hpp"
#include "paper/outer/scales_numeric.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/calendars.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::outer {

inline constexpr std::uint32_t kNoScaleKey = std::numeric_limits<std::uint32_t>::max();

/// [Q-67] Doublings added once at the top of a piece's table: the Sec 7.3 grid reaches the
/// first point >= 4 max R; a block edge is at most 4 K^2 R_v = 2^10 R_v (l.2101-2102), i.e. 8
/// more octaves, plus one point as the search sentinel. The paper says "a fixed number".
inline constexpr unsigned kSeparateGridExtraOctaves = 9;

/// F from the exact h = b + 1 for b = base^exponent held in multiword integers (S1, Q-71):
/// log2 F = 7 + p + 4k + bitlen(h^2 (2h + 1)), as outer_scale_factor. Integer word work on
/// sizes (Sec 2 l.302-304). Throws std::invalid_argument if base = 0 or b would need more than
/// 2^20 bits.
OuterScaleFactor outer_scale_factor_power_separator(unsigned log2_p, std::uint64_t base,
                                                    unsigned exponent, unsigned log2_k);

/// F of Mode A: b = L_N^64 (Sec 7.1 l.1644), K = 16. Never derived from the saturated word
/// that OwnerDecomposition::hop_bound holds.
OuterScaleFactor paper_outer_scale_factor(unsigned log2_p, std::uint64_t n);

/// S0 (l.2077): per dictionary index of `cls`, the least input weight of a macro edge carrying
/// that route name; infinity() if no macro edge does. Throws std::invalid_argument if `cls`
/// does not belong to `g` or a macro edge has a name outside the dictionary.
template <class W>
std::vector<W> outer_name_weights(const CsrGraph<W>& g,
                                  const owner::OuterEdgeClassification<W>& cls);

/// S3 for one name group (l.2130-2131): flags[i] = 1 iff F^2 s[i] >= a, decided by median
/// pruning over the scales (pivots::median_prune; the test is monotone in s): at most
/// floor(log2 n) + 1 threshold tests, added to *tests if given. Throws std::invalid_argument
/// unless a and every tested scale are positive and finite.
template <class W>
std::vector<std::uint8_t> outer_alignment_by_median_pruning(const OuterScaleFactor& f, W a,
                                                            std::span<const W> s,
                                                            std::uint64_t* tests = nullptr);

enum class OuterBlockKind : std::uint8_t {
  kSingleton,   ///< one state: only its finite theta join (l.2083-2085)
  kAligned,     ///< F^2 s_U >= a_U: global grid of a_U
  kUnaligned,   ///< finite attachment, F^2 s_U < a_U: separate grid
  kUnattached,  ///< theta infinite on U: separate grid (Q-66)
};

/// S2, S3, S5 for one block of the outer states. Only `kind` and `size` are meaningful for a
/// singleton block.
template <class W>
struct OuterBlockScale {
  OuterBlockKind kind = OuterBlockKind::kSingleton;
  std::size_t size = 0;                   ///< states of the block
  W s{};                                  ///< s_U
  Vertex s_vertex = kNoVertex;            ///< least vertex with r_v = s_U
  StateId z = kNoState;                   ///< z_U, least state minimising theta; kNoState if none
  std::uint32_t name = owner::kNoName;    ///< route name (dictionary index) of B0(z_U)
  std::uint32_t key = kNoScaleKey;        ///< grid key of `name`
  W t{};                                  ///< t_U
  std::uint32_t piece = owner::kNoPiece;  ///< separate grid: graph piece of s_vertex
  std::size_t table_index = 0;            ///< separate grid: t_U = piece_table[piece][table_index]
};

enum class OuterEventGrid : std::uint8_t {
  kAtom,      ///< block MST edge below t_U: inside an atom (index = atom id)
  kSeparate,  ///< class [t_U 2^j, t_U 2^(j+1)) of the block's separate grid (index = j)
  kGlobal,    ///< class [points[i], points[i+1]) of the global grid (index = i)
};

/// S7: the class of one join record. For kAtom, lower = 0 and upper = t_U.
template <class W>
struct OuterEvent {
  OuterEventGrid grid = OuterEventGrid::kGlobal;
  std::size_t index = 0;
  W lower{};
  W upper{};
};

/// Output of S0-S7.
template <class W>
struct OuterScales {
  OuterScaleFactor factor;                 ///< F as given
  std::vector<W> name_weight;              ///< S0, per dictionary index
  std::vector<std::uint32_t> name_key;     ///< per dictionary index: grid key, or kNoScaleKey
  std::size_t num_name_keys = 0;           ///< corridor c of T_A has key num_name_keys + c
  std::vector<W> key_weight;               ///< per key: the name weight or the corridor weight
  std::vector<OuterBlockScale<W>> blocks;  ///< per block of the outer states
  OuterGlobalGrid<W> grid;                 ///< S4; supports are indexed by key
  /// Per graph piece of the calendars: origin * 2^j in input units, j = 0 .. (extended top);
  /// empty for a piece no separate-grid block uses.
  std::vector<std::vector<W>> piece_table;
  std::vector<std::uint32_t> state_atom;  ///< per state: atom id or kNoAtom
  std::vector<ReplaceableAtom<W>> atoms;  ///< S6, numbered by their least state
  std::vector<std::uint32_t> atom_block;  ///< per atom
  std::vector<OuterEvent<W>> events;      ///< S7, per join record of OuterStates::joins
  std::uint64_t alignment_tests = 0;      ///< threshold tests made by S3
};

/// S0-S7 of Sec 8.2 on the outer states `os` of (g, od, cls, aux) with the calendars `cal` of
/// Alg 7 (OwnerDistancesResult::calendars) and the scale factor `f`.
/// Throws std::invalid_argument if the inputs do not belong together (sizes), and
/// std::logic_error if a relation listed in the file comment fails.
/// Counters (SSSP_DIAGNOSTICS) "outer_scales.*": calls, scale_factor_log2, name_keys,
/// blocks_singleton, blocks_aligned, blocks_unaligned, blocks_unattached, alignment_tests,
/// grid_points, grid_components, atoms, atom_states, atom_edges, events_atom, events_separate,
/// events_global.
template <class W>
OuterScales<W> build_outer_scales(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                  const owner::OuterEdgeClassification<W>& cls,
                                  const owner::AuxGraph<W>& aux, const owner::Calendars<W>& cal,
                                  const OuterStates<W>& os, const OuterScaleFactor& f,
                                  Diagnostics* diag = nullptr);

/// Checks an output of build_outer_scales against Sec 8.2 and Lemma 8.4 (note 04 Sec 2.1-2.8):
/// new names (a <= w < 2P(2h+1) a, attained; a <= B0 <= theta < F a on finite joins), block
/// kinds and scales (s_U = least original radius, z_U, theta >= s_U), alignment (F^2 s >= a;
/// aligned a / F^2 <= s <= F a; unaligned: MST edges below a / 2 and one theta), the grid
/// (increasing, doubling inside a component) and the supports (theta and corridor weights
/// inside their own), the cutoffs (s / 2 < t <= s; a grid point of the name; origin * 2^j of
/// the piece), the atoms (components of the MST edges below t_U; fewer than b states; one
/// block, one owner, one radius >= t_U; injective projection; internal edges = the
/// within-window edges inside; a spanning tree of edges below t_U; finite theta >= t_U), and
/// the classes (each brackets its weight; unaligned 2 p_top <= g). Returns one message per
/// violation.
template <class W>
std::vector<std::string> check_outer_scales(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                            const owner::OuterEdgeClassification<W>& cls,
                                            const owner::AuxGraph<W>& aux,
                                            const owner::Calendars<W>& cal,
                                            const OuterStates<W>& os, const OuterScales<W>& sc);

// ---------------------------------------------------------------- S8, first half: the merge plan

inline constexpr std::uint32_t kNoBlock = std::numeric_limits<std::uint32_t>::max();

/// One merge class (l.2094-2099): the join records merged simultaneously into nodes of height
/// `height` (the class lower boundary, Q-68); arcs shorter than `upper` connect the children.
template <class W>
struct OuterMergeClass {
  OuterEventGrid grid = OuterEventGrid::kGlobal;  ///< kSeparate or kGlobal
  std::uint32_t block = kNoBlock;                 ///< separate grid: its block
  std::size_t index = 0;                          ///< class index as in OuterEvent
  W height{};
  W upper{};
  std::vector<std::uint32_t> joins;  ///< indices into OuterStates::joins.joins, increasing
};

/// What the class-wise merge of Sec 8.2 (l.2094-2099) starts from:
///   * the parts of Def 5.12: one per state outside the atoms and one per atom (Q-74), numbered
///     by their least state;
///   * D, the state graph of the outside arcs: every arc of the outer states except those with
///     both endpoints in one atom (they are the atom's internal graph, Q-53 / Q-74);
///   * the classes in processing order: [Q-68] the separate grids first, block by block in
///     increasing class, then the global classes in increasing grid index. Any order consistent
///     with the heights along each root path gives the same forest: an unaligned block is
///     complete below its first global class and an unattached block meets no other join.
template <class W>
struct OuterMergePlan {
  std::vector<std::uint32_t> part;  ///< per state
  std::size_t num_parts = 0;
  std::vector<StateId> part_state;       ///< per part: its state, kNoState for an atom
  std::vector<std::uint32_t> part_atom;  ///< per part: its atom, or kNoAtom
  StateGraph<W> graph;                   ///< D
  std::vector<ArcId> arc_of;             ///< arc of D -> arc of OuterStates::states.graph
  std::vector<OuterMergeClass<W>> classes;
};

/// The merge plan of `sc` = build_outer_scales(.., os, ..). Counting sorts only (Lemma 3.1).
/// Throws std::invalid_argument if `sc` does not belong to `os`.
template <class W>
OuterMergePlan<W> build_outer_merge_plan(const OuterStates<W>& os, const OuterScales<W>& sc);

// ------------------------------------------------------- S8, second half: the partition hierarchy

/// The partition hierarchy of Sec 8.2 (l.2094-2099; Def 5.12), the input of Lemma 5.15.
/// Per-node vectors are indexed by the nodes of `hierarchy` (leaves included, with unspecified
/// values).
template <class W>
struct OuterPartition {
  /// graph = D (outside arcs), projection = that of the outer states, atoms = those of S6,
  /// leaves = the parts (node id = part id of the merge plan), c0 = 2.
  PartitionHierarchy<W> hierarchy;
  std::vector<ArcId> arc_of;          ///< arc of D -> arc of OuterStates::states.graph
  std::vector<NodeId> leaf_of_state;  ///< state -> the leaf of its part
  std::vector<W> node_upper;          ///< internal nodes: the class upper end (links are below it)
  std::vector<OuterEventGrid> node_grid;  ///< internal nodes: kSeparate or kGlobal
  std::vector<NodeId> roots;              ///< parentless nodes in condensation order (H4)
};

/// S8 (l.2094-2099): processes the classes of the merge plan with the simultaneous merges of
/// Sec 3.3 on the parts and one auxiliary vertex per pivot (assembly::ClassMerger: auxiliary
/// nodes without states are suppressed, unary nodes compressed), then extracts the connecting
/// trees by the child-quotient searches of Sec 5.4 and orders the children by the arcs lighter
/// than each node (assembly::link_and_order; every arc is assigned by the parts containing its
/// endpoints and keeps its actual endpoint states). [Q-68] A node of the class [a, 2a) has
/// height a and links strictly below 2a, hence c0 = 2 as in Lemma 5.10.
/// Throws std::invalid_argument if `sc` does not belong to `os`, and std::logic_error if a
/// height decreases towards the root, an atom's cutoff exceeds its parent's height (Def 5.12),
/// the children of a node are not connected by arcs below its class upper end, or the light
/// arcs between children (or the arcs between roots) form a cycle (Lemma 8.3 would fail).
/// Pushes no FallbackRecord. Counters "outer_scales.*": partition_calls, partition_parts,
/// partition_nodes (internal), partition_links, partition_roots, separate_classes,
/// global_classes; plus the dsu.* counters of the merges.
template <class W>
OuterPartition<W> build_outer_partition(const OuterStates<W>& os, const OuterScales<W>& sc,
                                        Diagnostics* diag = nullptr);

}  // namespace sssp::paper::outer
