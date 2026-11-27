#pragma once
// The outer construction: Algorithm 8 l.8-15 of arXiv 2611.34911v1 (Sec 8.3, txt l.2167-2196;
// proof of Theorem 1 l.2158-2217), on a connected graph with positive weights (the contracted
// source component of Alg 8 l.2). It chains the components of Sec 5-8 and adds no step of its
// own beyond the per-vertex minimum of l.15:
//
//   l.8   MST, owner cells, deletion of tree-dominated nontree edges   owner_decomposition
//   l.9   NameRoutes on the projected interowner edges (randomised)    classify_interowner_edges
//   l.10  radial / macro, G' and T'                                    build_aux_graph
//   l.11  OwnerDistances (Alg 7): J                                    owner_distances
//   l.12  deletion of radial edges with w > J_x + J_y + d_e            build_outer_states
//   l.13  outer states and joins (Sec 8.1)                             build_outer_states
//         scales, grids, atoms, partition hierarchy (Sec 8.2)          build_outer_scales,
//                                                                      build_outer_partition
//         capped local hierarchies (Lemma 5.15)                        substitute_atoms
//   l.14  seed the native inward copy of s (inside its atom: the       redirect_offers
//         native local subcopy)
//   l.15  Evaluate; minimum emitted value per vertex                   evaluate
//
// Alg 8 l.16-17 (restoration, Lemma 8.5) is restore_predecessor_tree, called by the entry points.
//
// [Q-73] l.15: the minimum runs over every state of the substituted hierarchy that projects to
// the vertex (ordinary copies, the special states v_in / v_out, local subcopies). No state
// underestimates (Lemma 4.6) and one copy is protected. A vertex without a finite emission is a
// violated invariant: std::logic_error, never infinity (the input is connected).
// [Q-71] F is formed from the exact b = base^exponent in multiword integers
// (outer_scale_factor_power_separator), never from the saturated word held by the decomposition.
// [Q-72] One pivot (always the case with the paper's b for N < 2^64) is the paper's
// construction, not a fallback: no macro edge, every nonsingleton block unattached, J = the
// distances from the pivot by Alg 7.
//
// Parameters (Q-04: one implementation). Mode A passes kPaperOuterParams (b = L_N^64, the
// NameRoutes constants of Alg 5, Alg 4 with 2 L_k^4 and b = L_k^4) and never anything else.
// Mode B (src/experimental/paper_experimental/sssp.hpp) passes other time-only constants;
// K = 16, the supports and F are not parameters.
//
// Fallback policy: this function never replaces a stage by a direct search and pushes no
// FallbackRecord of its own; an unsupported input throws. The records that appear in `diag` are
// those of its components: Alg 4's base case on every Sec 7.3 atom (Alg 7 l.6) and on every
// Sec 8.2 atom (Lemma 5.15), and NameRoutes' retry cap (D-10) if a Mode B caller sets one.
// [Q-69] One vertex can lie in an Alg 7 atom and, through its copies, in several Sec 8.2 atoms,
// so the vertex counts of these records can sum to more than N; fallback_shares() reports them
// per reason and never adds across reasons.
//
// Numerics (Q-06, Q-55 plan A): exactness is established for std::uint64_t and for dyadic
// double inputs whose sums are exact. Components throw (std::overflow_error, std::domain_error,
// std::logic_error) when a quantity they need is not representable; nothing is retried by
// another method. Observed in the tests of note 04 PL-7, PL-10 and PL-11 (Mode B forcing):
//   * std::uint64_t: the Sec 7.3 calendars reject R >= 2^54 (D-20); Sec 4.2 order_children
//     and the Sec 5.1 root distances throw on some inputs with weights near 2^62 whose
//     distances all fit (open fixes P3-CHILDORDER-OVERFLOW, P3-OWNERDECOMP-OVERFLOW);
//   * non-dyadic double: most runs are rejected, mainly by the tight-arc test of Lemma 8.5
//     (restore_predecessor_tree, after this function); every result that was returned passed
//     the tree verification within the documented tolerance.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/atom_replacement.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/outer/outer_states.hpp"
#include "paper/outer/scales.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/owner_distances.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::outer {

/// The separator size b = base^exponent of the outer owner decomposition, kept as a power so
/// that h = b + 1 and F can be formed exactly when b does not fit a word (Q-71).
struct OuterSeparator {
  std::uint64_t base = 1;
  unsigned exponent = 1;
  friend constexpr bool operator==(const OuterSeparator&, const OuterSeparator&) = default;
};

/// Sec 7.1 l.1644: b = L_N^64.
OuterSeparator paper_outer_separator(std::uint64_t n);

/// The constants of Alg 8 l.8-15 that only affect time. Defaults are the paper's.
struct OuterConstructionParams {
  /// nullopt: paper_outer_separator(N). A value: Mode B (deviation D-7).
  std::optional<OuterSeparator> separator;
  /// Alg 5 (l.9). Defaults: 16t, 32t, the paper's P and H, no retry cap.
  pivots::NameRoutesParams name_routes{};
  /// Alg 4 on every atom of Sec 7.3 (Alg 7 l.6) and of Sec 8.2 (Lemma 5.15).
  LocalHierarchyParams local = kPaperLocalParams;
};

/// The paper's values; Mode A never uses anything else.
inline const OuterConstructionParams kPaperOuterParams{};

/// Every intermediate object of one run, in the order of Alg 8, for the stage checkers.
template <class W>
struct OuterConstruction {
  OuterSeparator separator;               ///< b as used
  OwnerDecomposition<W> od;               ///< l.8
  owner::OuterEdgeClassification<W> cls;  ///< l.9-10
  owner::AuxGraph<W> aux;                 ///< l.10
  owner::OwnerDistancesResult<W> jr;      ///< l.11
  OuterStates<W> os;                      ///< l.12-13 (Sec 8.1)
  OuterScaleFactor factor;                ///< F (Sec 8.2 l.2080)
  OuterScales<W> sc;                      ///< l.13 (Sec 8.2, S0-S7)
  OuterPartition<W> part;                 ///< l.13 (Sec 8.2, S8)
  AtomSubstitution<W> sub;                ///< l.13 (Lemma 5.15)
  std::vector<StateOffer<W>> offers;      ///< l.14, after redirection
  std::vector<W> state_value;             ///< l.15: emitted value per state of sub.hierarchy
  std::vector<W> dist;                    ///< l.15: minimum per vertex
};

/// Alg 8 l.8-15 on (g, source) keeping every stage. `g` must be connected with positive finite
/// weights and no loops (Alg 8 l.2); parallel edges are allowed.
/// Throws std::out_of_range if source is not a vertex; std::invalid_argument for an invalid
/// graph or parameters (as owner_decomposition, name_routes, local_hierarchy; separator base 0
/// or exponent 0); std::logic_error if a vertex receives no finite value (Q-73) or a component
/// finds one of its invariants violated; std::overflow_error / std::domain_error as the
/// components do when a length they need is not representable.
/// `diag` (may be null) receives the counters and FallbackRecords of every component, plus
/// (SSSP_DIAGNOSTICS) "outer.constructions", "outer.vertices", "outer.edges", "outer.pivots",
/// "outer.macro_edges", "outer.radial_edges_deleted", "outer.states", "outer.new_states",
/// "outer.atoms", "outer.scale_factor_log2".
template <class W>
OuterConstruction<W> run_outer_construction(
    const CsrGraph<W>& g, Vertex source, model::WordRng& rng,
    const OuterConstructionParams& params = kPaperOuterParams, Diagnostics* diag = nullptr);

/// Alg 8 l.8-15: the exact distances from `source` (run_outer_construction(..).dist).
template <class W>
std::vector<W> outer_construction(const CsrGraph<W>& g, Vertex source, model::WordRng& rng,
                                  const OuterConstructionParams& params = kPaperOuterParams,
                                  Diagnostics* diag = nullptr);

/// [Q-69] Fallback share of one reason: the sums of FallbackRecord::vertices and ::edges over
/// the records with that reason, and their count. The sums are to be read against N and m of
/// the graph the records refer to; `vertices` can exceed N. Shares of different reasons are
/// never added.
struct FallbackShare {
  std::string reason;
  std::uint64_t records = 0;
  std::uint64_t vertices = 0;
  std::uint64_t edges = 0;
};

/// Per distinct reason of diag.fallbacks, in order of first appearance.
std::vector<FallbackShare> fallback_shares(const Diagnostics& diag);

/// Appends the diagnostics of one stage: operation counts and counters are summed, fallbacks
/// kept in order.
void absorb_diagnostics(Diagnostics& into, Diagnostics&& from);

extern template OuterConstruction<std::uint64_t> run_outer_construction(
    const CsrGraph<std::uint64_t>&, Vertex, model::WordRng&, const OuterConstructionParams&,
    Diagnostics*);
extern template OuterConstruction<double> run_outer_construction(const CsrGraph<double>&, Vertex,
                                                                 model::WordRng&,
                                                                 const OuterConstructionParams&,
                                                                 Diagnostics*);
extern template std::vector<std::uint64_t> outer_construction(const CsrGraph<std::uint64_t>&,
                                                              Vertex, model::WordRng&,
                                                              const OuterConstructionParams&,
                                                              Diagnostics*);
extern template std::vector<double> outer_construction(const CsrGraph<double>&, Vertex,
                                                       model::WordRng&,
                                                       const OuterConstructionParams&,
                                                       Diagnostics*);

}  // namespace sssp::paper::outer
