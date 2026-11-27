#pragma once
// The outer state graph and its threshold components: Sec 8.1 of arXiv 2611.34911v1
// (txt l.1966-2001; Lemma 8.1 l.2002, Lemma 8.2 l.2027, Lemma 8.3 l.2053), i.e. the part of
// Alg 8 l.12 that precedes the numerical scales of Sec 8.2.
//
// Input: the connected positive graph G with its outer owner decomposition (Sec 5.1 with the
// outer b, Sec 7.1), the edge classification of Sec 7.1 (radial / macro, d_e, names), the
// auxiliary graph G' with T' (Sec 7.2), and J from OwnerDistances (Alg 7) in the calendar units
// of Lemma 7.3, J' = 2^7 J (OwnerDistancesResult::scaled_distance, indexed by aux vertex).
//
// Steps:
//   1. Deletion (l.1967-1968): a radial edge with w > J_x + J_y + d_e is deleted (it becomes
//      kOmitted). [Lead decision, D-21] The test is made in calendar units, exactly:
//      2^7 w > J'_x + J'_y + 2^7 d_e (outer_radial_deleted). The unscaled J is never formed
//      inside Sec 8. Owner tree edges survive by J_v >= r_v >= w_e (l.1969-1970); this is not
//      enforced but checked by the tests, as is w <= 2 (J_x + J_y) for survivors (Lemma 8.1).
//   2. Two-shell construction (l.1973-1983): Def 5.4 (build_radial_states) on G with the
//      coordinate J' (shells and the near test K^2 are invariant under the common factor 2^7),
//      K = 16 (kOuterShellRatioLog2 = 4), and the outer classification after deletion; macro
//      keys are the route names (dictionary indices of OuterEdgeClassification::name), the keys
//      Sec 8.2 replaces by least macro weights (l.2077).
//   3. Directed owner-tree states (l.1984-1997): every nonpivot v gets v_in (inward) and v_out
//      (outward), appended after the ordinary states, each a singleton block (l.1998). With
//      parent(v) the next vertex on the owner tree T[v, p(v)] (the T' parent, read through the
//      physical owner edge) and w its edge:
//        kTreeEntry   every ordinary inward copy of v -> v_in, length 0
//        kInwardTree  v_in -> parent_in, or -> the ordinary inward copy of the pivot, length w
//        kOutwardTree parent_out, or the ordinary outward copy of the pivot, -> v_out, length w
//        kTreeExit    v_out -> native outward copy of v, length 0
//      These arcs carry RadialArcKind::kOwnerTree in `states.arc_kind`; macro arcs keep
//      ordinary tails and heads (l.1990).
//   4. Joins (l.1998-2000): compute_joins (Sec 5.3) on the augmented graph and the outer owner
//      decomposition: I0, O0, I*, O*, theta in macro-free block order; block MST joins,
//      state-to-owner joins at theta, corridor joins of T_A.
// s_U = min_{v in U} r_v and the cutoffs t_U (l.2000-2001) belong to Sec 8.2 (P3-OUTER-SCALES).
//
// Lemma 8.1: survivors have w <= 2 (J_x + J_y); globally shortest radial runs swing by at most
// 16; the ordinary graph realizes every geodesic; the augmented macro-free quotient is acyclic.
// Lemma 8.2: B0 and theta are finite together, r <= theta and B0 <= theta <= 4 K^2 h B0.
// Lemma 8.3: the joins of weight <= W (< W) give the SCCs of the arcs of length <= W (< W).
//
// [Interpretation] "Ordinary inward copy of v may enter v_in" (l.1985) is read as: from every
// ordinary inward copy of v (native and lower), as for the far and pivot rules (Q-32).
// [Interpretation] The owner tree of the special states is T' (aux.tree_parent), i.e. the
// uncut MST edges, with the pivots as roots; for a vertex whose T' parent is the root a the
// physical parent is the pivot at the other end of its owner edge.
// [Supplement] Special states get window 0; their phase is inward for v_in and outward for v_out.
//
// Numerics (Q-06, D-21): for std::uint64_t the deletion test is evaluated in unsigned 128-bit
// arithmetic (2^7 w by a shift, sums of at most three 71-bit values), so no edge is dropped or
// rejected for not fitting 64 bits (Q-58). For double the doublings by 2^7 are exact and the
// sums round; exactness is established for integer and dyadic weights (Q-55 plan A). A
// non-finite 2^7 w raises std::overflow_error.
//
// Fallback policy: Sec 8.1 never substitutes a search; it pushes no FallbackRecord. (The atoms of
// Alg 7 record their Alg 4 base cases in OwnerDistances.)

#include <cstddef>
#include <cstdint>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "paper/owner/aux_metric.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::outer {

/// K = 2^4 = 16 for the outer two-shell construction (Sec 8.1 l.1977). Lemma 8.1's swing factor
/// needs exactly this value, so it is not a parameter (Mode B too).
inline constexpr unsigned kOuterShellRatioLog2 = 4;

/// Finer kind of every arc of OuterStates::states.
enum class OuterArcKind : std::uint8_t {
  kOrdinary,     ///< an arc of Def 5.4 (its kind is in states.arc_kind)
  kTreeEntry,    ///< ordinary inward copy of v -> v_in, length 0
  kInwardTree,   ///< v_in -> parent_in or the pivot's ordinary inward copy, physical owner edge
  kOutwardTree,  ///< parent_out or the pivot's ordinary outward copy -> v_out, physical owner edge
  kTreeExit,     ///< v_out -> native outward copy of v, length 0
};

/// Step 1 (l.1967) in calendar units (D-21): true iff 2^7 w > J'_x + J'_y + 2^7 d_e, where
/// J' = 2^7 J. Exact for std::uint64_t (128-bit); for double the doublings are exact and the
/// sums round (Q-06). Throws std::overflow_error for double if 2^7 w or 2^7 d_e is not finite.
template <class W>
bool outer_radial_deleted(W w, W detour, W scaled_jx, W scaled_jy);

template <class W>
struct OuterStates {
  std::size_t hop_bound = 0;          ///< h of the outer owner decomposition
  std::vector<std::uint8_t> deleted;  ///< per input edge: 1 iff a radial edge deleted in step 1
  std::size_t num_deleted = 0;

  /// Def 5.4 on (G, J', classification after deletion, K = 16), followed by the special states
  /// and their arcs. states.coordinate[v] = J'_v (0 at pivots); states.classification is the
  /// classification after deletion with the route names as macro keys.
  RadialStateGraph<W> states;
  std::size_t num_ordinary_states = 0;  ///< states [0, this) are ordinary, the rest special
  std::size_t num_ordinary_arcs = 0;    ///< arcs [0, this) are ordinary
  std::size_t num_ordinary_blocks = 0;  ///< blocks [0, this) are ordinary
  std::vector<OuterArcKind> arc_kind;   ///< per arc of `states`

  // Per input vertex; kNoState / kNoVertex / kNoEdge at pivots.
  std::vector<StateId> special_inward;    ///< v_in
  std::vector<StateId> special_outward;   ///< v_out
  std::vector<Vertex> owner_parent;       ///< next vertex towards p(v) on T[v, p(v)]
  std::vector<EdgeId> owner_parent_edge;  ///< the physical owner edge to owner_parent

  /// Step 4: Sec 5.3 on the augmented graph (num_pivots = |A| auxiliary pivot vertices).
  ThresholdJoins<W> joins;

  [[nodiscard]] bool is_special(StateId x) const noexcept { return x >= num_ordinary_states; }
};

/// Sec 8.1 (Alg 8 l.12, first part). `scaled_j` = OwnerDistancesResult::scaled_distance (2^7 J
/// per aux vertex). Throws std::invalid_argument if the inputs do not belong together (sizes,
/// aux of a different graph) or J' is not a valid coordinate (J'_a != 0, or a nonpivot with
/// J' = 0 or infinite: Lemma 7.2 gives J >= r > 0); std::overflow_error as outer_radial_deleted;
/// std::logic_error if the macro-free block quotient has a cycle (compute_joins; impossible by
/// Lemma 8.1).
/// Counters (SSSP_DIAGNOSTICS) "outer_states.*": calls, edges_scanned, deleted_radial_edges,
/// special_states, special_arcs, states, arcs; plus those of build_radial_states and
/// compute_joins.
template <class W>
OuterStates<W> build_outer_states(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                  const owner::OuterEdgeClassification<W>& cls,
                                  const owner::AuxGraph<W>& aux, const std::vector<W>& scaled_j,
                                  Diagnostics* diag = nullptr);

}  // namespace sssp::paper::outer
