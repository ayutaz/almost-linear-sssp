#pragma once
// Sec 7.4 "Local Searches, Clipping, and External Arrivals" and Algorithm 7 OwnerDistances
// (arXiv 2611.34911v1, txt l.1780-1963; Alg 7 l.1843-1867, Lemma 7.7 l.1871, Lemma 7.8 l.1919,
// Lemma 7.9 l.1944), i.e. Alg 8 l.11, Mode A (Mode B passes other Alg 4 parameters for the atoms
// and reaches nontrivial auxiliary graphs through other separator sizes b; Q-04).
//
// Input: the auxiliary graph G' and its owner-spoke MST T' rooted at the common pivot a
// (AuxGraph, Sec 7.2). Output: J_v = d_{w'}(v, a) for every auxiliary vertex (J_a = 0).
//
// Steps (Alg 7):
//   l.2  build_calendars (Sec 7.3, Lemma 7.6).
//   l.3  find_atoms (Sec 7.3 l.1705-1709, Lemma 7.5).
//   l.4  build_families: an atom belongs to the family (P, q) of its representative's tree piece
//        and width (l.1781-1782).
//   l.5  Within a family the atoms are ordered by the T' preorder of their representatives
//        (l.1782-1783; a counting pass over one DFS preorder, i.e. integer sorting). The gap
//        g(i, j) of consecutive representatives is their T' distance R_i + R_j - 2 R_lca with
//        the LCA of Lemma 3.9 (l.1785). The family horizon is max_Y (R_rep(Y) + L_Y)
//        (l.1786-1787). Activation lists (l.1796-1800): the family is appended to the starting
//        marker of every cell of its (C, q) calendar from cell 0 through the last cell whose
//        start is <= the horizon, or to the calendar's end.
//   l.6  One edge scan over G' (l.1816-1818) sorts every edge into: a root edge (endpoint a), an
//        internal edge of one atom (both ends in Y, heavier chords included), or an interatom
//        exit listed at its tail. For an edge between graph pieces only the direction from the
//        smaller range to the larger (the tail of its interpiece marker arc) is kept
//        (l.1823-1824). LocalHierarchy (Alg 4) runs on every G'[Y] and is capped at the width q
//        (cap_local_hierarchy, Lemma 5.13; the discovery tree has every edge < 2q, so Y is a
//        replaceable atom at cutoff q, l.1829-1831).
//   l.7  normalize (Lemma 4.3) and order_children (Lemma 4.5) on every capped hierarchy.
//   l.8  One split-findmin layout (Lemma 3.11): family order, then representative preorder
//        within the family, then each atom's local leaf layout (PortionedEvaluator with the
//        shared structure), all keys infinite, split at every family boundary.
//   l.9  initial_offer R_v at the native inward copy of every nonroot v, and the direct offers
//        w'(a v) of every root edge at the native copy of v (Q-47).
//   l.10-13  For every marker in marker_order (a topological order of the marker DAG) and every
//        family record (F, cell [a0, b0)) of its activation list:
//        * l.12 (l.1801-1815): every live run of unregistered atoms of F is checked by one
//          interval-minimum query; if its minimum z is due (a0 <= z < b0) the atom i holding it
//          is registered (two splits isolate it), and each of the at most two new nonempty
//          subruns gets the proposal u = max{b0, z + L_i + g(i, j)} at the native copy of the
//          representative v of its boundary atom j, discarded if u > R_v (and ignored by the
//          evaluator unless it lowers the pending key). New subruns are checked again in the
//          same cell. The registered atom gets the inclusive deadline z + L_Y (l.1825).
//        * l.13 (l.1835-1841): every registered atom of F that has not finished is visited by
//          PortionedEvaluator::visit_cell([a0, b0)) (Lemma 7.7). Every emission at a copy of u
//          relaxes the kept exits of u: suppressed if the emitted value exceeds R_u, discarded at
//          the head v if value + w' exceeds R_v (l.1821-1823, Lemma 7.8); otherwise an
//          external_offer at the native copy of v.
//   l.14 J_v = the least emitted value over the copies of v (J_a = 0).
//
// Units (Lemma 7.3 l.1723-1724, "after scaling all quantities by a fixed power of two"): the
// calendars store every value as 2^7 times the length it stands for (calendars.hpp). Alg 7 runs
// entirely in those units: every auxiliary length (w', R, L_Y, gaps) is scaled by calendar_scale
// before use, so the cells [i q, (i + 1) q) are exact for integer lengths. The scaled distances
// 2^7 J are returned in scaled_distance; distance = J is obtained by seven exact halvings.
// [D-21] The halvings are a deviation: halving is not one of the Sec 2 operations. Every emitted
// value is a sum of scaled lengths, so for integers the low seven bits are zero (checked) and for
// double the halvings are exact for normal values. Sec 8 consumes scaled_distance, taken before
// the halvings.
//
// Interpretations (open questions) and deviations (PAPER_MAPPING §3):
//   [Q-47] "direct offers from the root" (l.1860): the w' of every G' edge incident to a, posted
//     at the native copy of its other endpoint (initial offers, before the first marker).
//   [Q-27, Q-28] slot resumption and padding are those of PortionedEvaluator.
//   [Supplement] the T' preorder visits children in the BFS order aux.tree_order (the paper
//     fixes no child order; any DFS preorder gives the tour bound of l.1783-1785).
//   [D-22] g(i, j) is materialised as one real (R_i - R_lca) + (R_j - R_lca) by two
//     subtractions of nonnegative differences, so that the proposal z + L_i + g(i, j) is an
//     ordinary length (exact for integers; for double rounded, Q-06 and Q-55 option A; the
//     paper's "signed pairs of sums", Sec 2 l.299, would avoid the subtraction).
//   [Interpretation] The inherited physical deadline of l.1833 ("the smaller of its own and
//     inherited physical deadlines") is read as the evaluator's internal cap min(c, D0 + S) of
//     the opened root, c being the atom's deadline (PortionedEvaluator, l.1837-1838); the atom's
//     own deadline c = z + L_Y (l.1825) is set once, at registration, and never changed
//     afterwards.
//   [Q-55, option A] For double, R_v (summed down T') and R_rep + L_Y (another summation order)
//     can differ by rounding (about one ulp). The place where this rounding bites is the atom
//     deadline z + L_Y (l.1825), which can close the atom one ulp before v is due; the family
//     horizon (l.1786-1787) is not the cause. Decision (lead, 2026-11-27): option A, i.e.
//     correctness is established for integer and dyadic weights only (Q-06). With non-dyadic
//     doubles a vertex that is never emitted raises std::logic_error instead of returning a
//     wrong J.
//
// Fallback policy: every atom goes through Alg 4 (local_hierarchy via cap_local_hierarchy),
// which pushes FallbackRecord{"local: k <= 2*L_k^4 (paper base case)", k, m} when it takes its
// base case (with Mode A parameters every atom does, Q-04). OwnerDistances itself never
// substitutes another algorithm and pushes no record of its own.
//
// Invariant checks (std::logic_error, every build): a family run whose minimum lies below the
// current cell (Lemma 7.8 violated), and a nonroot vertex without emission (Lemma 7.9 violated).
// PortionedEvaluator throws std::invalid_argument for an offer below a visited cell.

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "paper/local/local_hierarchy.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/calendars.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::owner {

inline constexpr std::uint32_t kNoFamily = 0xffffffffu;

/// Alg 7 l.4-5: families (P, q), representative preorder, consecutive gaps and horizons.
/// Lengths are scaled (2^7 x, calendars.hpp).
template <class W>
struct Families {
  std::vector<std::uint32_t> preorder;    ///< per aux vertex: T' preorder index (root 0)
  std::vector<std::uint32_t> family;      ///< per atom
  std::vector<std::uint32_t> position;    ///< per atom: index into `order`
  std::vector<std::uint32_t> offsets;     ///< family f: order[offsets[f] .. offsets[f + 1])
  std::vector<std::uint32_t> order;       ///< atoms, family by family, reps in T' preorder
  std::vector<std::uint32_t> tree_piece;  ///< per family: P
  std::vector<std::uint32_t> calendar;    ///< per family: the (C, q) calendar
  std::vector<W> horizon;                 ///< per family: max_Y (R_rep + L_Y), scaled
  /// per position k of `order`: g(order[k], order[k + 1]) if both are in one family, else 0.
  std::vector<W> gap;

  [[nodiscard]] std::size_t num_families() const noexcept { return tree_piece.size(); }
};

/// Alg 7 l.4-5 (Sec 7.4 l.1781-1787). The representatives are taken in T' preorder and then
/// stably counting-sorted by width calendar and by tree piece, so families are numbered by
/// (tree piece, calendar) and keep the preorder inside. Throws std::invalid_argument if `cal` /
/// `atoms` do not belong to `aux`, and std::overflow_error if a scaled length does not fit W.
/// Counter "owner_distances.families".
template <class W>
Families<W> build_families(const AuxGraph<W>& aux, const Calendars<W>& cal, const Atoms<W>& atoms,
                           Diagnostics* diag = nullptr);

struct OwnerDistancesOptions {
  /// Alg 4 parameters for every atom (Mode A: the paper's; Mode B: e.g. {1, 2, b}).
  LocalHierarchyParams local = kPaperLocalParams;
  /// Fill OwnerDistancesResult::arrivals (Lemma 7.8 arrival log).
  bool record_arrivals = false;
  /// Fill OwnerDistancesResult::emissions.
  bool record_emissions = false;
  /// Count Scan calls per local nominal slot (Lemma 7.7: at most twice) through EvalHooks.
  bool count_slot_scans = false;
};

enum class ArrivalKind : std::uint8_t {
  kInteratom,     ///< an exit relaxed after an emission (Lemma 7.8)
  kSupplemental,  ///< a family proposal u = max{b0, z + L_i + g(i, j)} (l.1803-1805)
};

/// One offer that lowered a pending key of an atom state (an "effective decrease", Lemma 7.8).
/// Values are scaled.
template <class W>
struct ArrivalRecord {
  ArrivalKind kind = ArrivalKind::kInteratom;
  Vertex tail = kNoVertex;  ///< emitting vertex (kNoVertex for a family proposal)
  Vertex head = kNoVertex;
  W value{};
  std::uint32_t marker = kNoMarker;       ///< marker being processed when the offer arrived
  std::uint32_t destination = kNoMarker;  ///< starting marker of the head's cell of `value`
  bool on_time = false;                   ///< `destination` had not been processed yet
};

template <class W>
struct OwnerDistancesResult {
  std::vector<W> distance;                 ///< J_v per aux vertex (input units), J_a = 0
  std::vector<W> scaled_distance;          ///< 2^7 J_v
  Calendars<W> calendars;                  ///< Alg 7 l.2
  Atoms<W> atoms;                          ///< Alg 7 l.3
  Families<W> families;                    ///< Alg 7 l.4-5
  std::vector<W> registration;             ///< per atom: z (scaled); infinity() if never registered
  std::vector<std::uint32_t> atom_cells;   ///< per atom: cells visited from its registration on
  std::vector<std::uint8_t> radial_atom;   ///< per atom: 1 iff Alg 4 took the radial branch
  std::size_t activation_records = 0;      ///< family records over all activation lists
  std::vector<ArrivalRecord<W>> arrivals;  ///< if record_arrivals
  std::vector<std::pair<Vertex, W>> emissions;  ///< (vertex, scaled value) if record_emissions
  std::uint32_t max_slot_scans = 0;             ///< if count_slot_scans
};

/// Alg 7 on (T', G') = aux. Throws std::invalid_argument for an inconsistent aux (via
/// build_calendars), std::overflow_error if 2^7 w' of an edge of G' (or 2^7 R) does not fit W,
/// std::logic_error if an invariant of Lemmas 7.8-7.9 fails (see the file comment), plus what
/// cap_local_hierarchy / local_hierarchy throw for unsupported Alg 4 parameters.
/// Counters "owner_distances.*" (SSSP_DIAGNOSTICS): calls, families, atoms, atom_states,
/// internal_edges, root_edges, exits, exits_wrong_direction, initial_offers, root_offers,
/// activation_records, markers, run_checks, run_checks_failed, registrations,
/// supplemental_proposals, supplemental_capped, supplemental_decreases, atom_cell_visits,
/// emissions, exit_relaxations, exits_suppressed, exits_capped, interatom_decreases; plus the
/// counters of build_calendars, find_atoms, cap_local_hierarchy, local_hierarchy, normalize,
/// order_children and PortionedEvaluator, and Alg 4's fallback records.
template <class W>
OwnerDistancesResult<W> owner_distances(const AuxGraph<W>& aux, Diagnostics* diag = nullptr,
                                        const OwnerDistancesOptions& options = {});

}  // namespace sssp::paper::owner
