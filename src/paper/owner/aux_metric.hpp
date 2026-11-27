#pragma once
// Sec 7.1 "Edge Classification and the Auxiliary Metric" and Sec 7.2 "The Owner Coordinate"
// (arXiv 2611.34911v1, txt l.1642-1690; Lemma 7.1 l.1653, Lemma 7.2 l.1668), i.e. Alg 8 l.9-10,
// Mode A (Mode B passes other separator sizes and NameRoutes parameters).
//
// Input: a connected positive graph G with N vertices and its owner decomposition (Sec 5.1,
// Def 5.1) with the outer separator size b = L_N^64 (Sec 7.1 l.1643-1644, outer_separator_size);
// r_v = b_T(v, p(v)), R_v = D_T(v, p(v)), h = b + 1. Only edges of the reduced graph are used
// (Alg 8 l.8 deletes the strictly tree-dominated nontree edges).
//
// Steps:
//   1. projected_interowner_edges (l.1646-1647): every reduced edge xy with p(x) != p(y) becomes
//      the pivot-routing edge p(x)p(y) of length A_e = R_x + w + R_y (Definition 6.1; the tree
//      walks x -> p(x) and y -> p(y) plus the edge make it a physical walk).
//   2. NameRoutes (Alg 5, Lemma 6.8) on that instance gives names lambda_e <= A_e and the power of
//      two P (l.1647-1649).
//   3. classify_outer_edges (l.1649-1652): d_e = P lambda_e, formed once per dictionary entry by
//      log2 P doublings (l.1649). An interowner edge is radial iff 2 d_e <= w, macro otherwise;
//      intraowner edges are radial with d_e = 0.
//   4. build_aux_graph (l.1663-1667): identify all pivots into the root a, keep the radial edges
//      with w' = w + d_e, discard the positive loops at a (radial edges between two pivots), keep
//      the owner spokes (the uncut MST edges) as T', rooted at a. J_v = d_{w'}(v, a).
//
// Lemma 7.1: every cut corridor edge is macro; a macro of weight w has lambda/(2h+1) <= w < 2P
// lambda. Lemma 7.2: r_v <= J_v <= R_v <= h r_v, |J_x - J_y| <= w' <= 3w/2 on radial edges, T' is
// an MST of G'. check_aux_metric checks these statements (and the construction) on an output.
//
// J itself is computed by OwnerDistances (Alg 7, P3-OWNER-DISTANCES). owner_coordinates_ref is a
// reference (Dijkstra, Lemma 3.2, on G') for tests and checkers only; it is NOT a step of the
// paper and must not replace Alg 7 in Mode A or Mode B.
//
// Model of computation (Sec 2): lengths are only compared, added and doubled. d_e = P lambda_e
// uses doublings; a doubling that does not fit the type saturates at WeightTraits<W>::infinity()
// (Q-38, paper/local/saturating.hpp). A saturated d_e exceeds every finite w, so the radial test
// 2 d_e <= w stays exact and the edge is macro; saturated values never enter G' (only radial edges
// do, and their d_e <= w/2). The integer h is used only in the saturating products of the
// checker. For double weights the sums A_e and w + d_e are rounded (Q-06: exactness is checked
// with integer and dyadic weights).
//
// Vertex numbering of G': the root a is 0 (kAuxRoot); the nonpivots follow in increasing input
// vertex order. Edges of G' follow the increasing input edge order and keep the input orientation.
// [Interpretation] The paper fixes neither numbering; this one makes G' and T' deterministic.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::owner {

inline constexpr Vertex kAuxRoot = 0;
inline constexpr std::uint32_t kNoInstanceEdge = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::uint32_t kNoName = std::numeric_limits<std::uint32_t>::max();

/// b = L_N^64 for the outer construction (Sec 7.1 l.1643-1644), saturating at 2^64 - 1 (for
/// every N >= 1 the value saturates, so the decomposition has a single pivot; Q-04, Q-48).
std::size_t outer_separator_size(std::size_t n);

/// Classification of an input edge (Sec 7.1 l.1650-1652). kOmitted: loops and edges outside
/// the reduced graph (Alg 8 l.8).
enum class OuterEdgeClass : std::uint8_t { kOmitted, kRadial, kMacro };

/// The pivot-routing instance of Sec 7.1 l.1646 (Definition 6.1): t = |A| pivots (pivot
/// indices of the owner decomposition), one edge (p(x), p(y), A_e) per reduced interowner edge
/// in increasing input edge order.
template <class W>
struct ProjectedEdges {
  pivots::PivotRoutingInstance<W> instance;
  std::vector<EdgeId> input_edge;            ///< instance edge -> input edge
  std::vector<std::uint32_t> instance_edge;  ///< input edge -> instance edge, or kNoInstanceEdge
};

/// Step 1. Throws std::invalid_argument if `od` does not belong to `g` (sizes differ) and
/// std::overflow_error if a sum A_e overflows.
template <class W>
ProjectedEdges<W> projected_interowner_edges(const CsrGraph<W>& g, const OwnerDecomposition<W>& od);

/// Steps 2-3 output. Edge-indexed vectors have size m (input edge ids).
template <class W>
struct OuterEdgeClassification {
  std::size_t hop_bound = 0;  ///< h = b + 1 (saturated as in OwnerDecomposition)
  unsigned p_log2 = 0;        ///< P = 2^p_log2
  std::vector<W> dictionary;  ///< names lambda (NameRoutes dictionary, strictly increasing)
  /// P lambda per dictionary entry (l.1649), infinity() if a doubling does not fit.
  std::vector<W> dictionary_detour;
  std::vector<OuterEdgeClass> edge_class;  ///< per input edge
  std::vector<std::uint8_t> interowner;    ///< per input edge: 1 iff reduced and p(x) != p(y)
  std::vector<std::uint32_t> name;         ///< per input edge: dictionary index, or kNoName
  std::vector<W> detour;                   ///< per input edge: d_e (0 if not interowner)
};

/// Step 3 from given route names (one name per edge of `projected.instance`). Throws
/// std::invalid_argument if the names do not match the instance (count, index outside the
/// dictionary, or lambda_e > A_e, which Definition 6.1 forbids). The walk bound
/// d_{G_P}(x, y) <= P lambda_e is the caller's responsibility (pivots::check_route_names).
/// Counters "aux_metric.*": interowner_edges, radial_interowner_edges, macro_edges,
/// intraowner_edges, omitted_edges, saturated_detours.
template <class W>
OuterEdgeClassification<W> classify_outer_edges(const CsrGraph<W>& g,
                                                const OwnerDecomposition<W>& od,
                                                const ProjectedEdges<W>& projected,
                                                const pivots::RouteNames<W>& names,
                                                Diagnostics* diag = nullptr);

/// Steps 1-3 (Alg 8 l.9 and the classification of l.10): NameRoutes on the projected
/// interowner edges with the caller's generator and parameters (defaults: the paper's), then
/// classify_outer_edges. NameRoutes' counters and fallback records (D-10 retry cap, Mode B only)
/// are added to `diag`.
template <class W>
OuterEdgeClassification<W> classify_interowner_edges(const CsrGraph<W>& g,
                                                     const OwnerDecomposition<W>& od,
                                                     model::WordRng& rng,
                                                     const pivots::NameRoutesParams& params = {},
                                                     Diagnostics* diag = nullptr);

/// Step 4 output: G' and T' (Sec 7.2). Aux-vertex-indexed vectors have num_vertices() entries.
template <class W>
struct AuxGraph {
  std::size_t hop_bound = 0;         ///< h
  CsrGraph<W> graph;                 ///< G': radial edges with w' = w + d_e, loops at a removed
  std::vector<EdgeId> input_edge;    ///< aux edge -> input edge
  std::vector<EdgeId> aux_edge;      ///< input edge -> aux edge, or kNoEdge
  std::vector<Vertex> aux_vertex;    ///< input vertex -> aux vertex (every pivot -> kAuxRoot)
  std::vector<Vertex> input_vertex;  ///< aux vertex -> input vertex (kNoVertex at the root)

  std::vector<EdgeId> tree_edges;        ///< T' (owner spokes) as aux edges, increasing
  std::vector<Vertex> tree_parent;       ///< parent in T' (kNoVertex at the root)
  std::vector<EdgeId> tree_parent_edge;  ///< aux edge to the parent (kNoEdge at the root)
  std::vector<Vertex> tree_order;        ///< BFS order of T' from the root (parents first)
  std::vector<W> radius;                 ///< r_v = root bottleneck in T' (0 at the root)
  std::vector<W> owner_distance;         ///< R_v = root distance in T' (0 at the root)

  [[nodiscard]] std::size_t num_vertices() const noexcept { return input_vertex.size(); }
};

/// Step 4. Throws std::invalid_argument if `cls` does not belong to `g` / `od`, std::logic_error
/// if an uncut tree edge is not radial or T' does not span G', and std::overflow_error if a sum
/// w + d_e overflows. Counters "aux_metric.*": aux_vertices, aux_edges, root_loops.
template <class W>
AuxGraph<W> build_aux_graph(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                            const OuterEdgeClassification<W>& cls, Diagnostics* diag = nullptr);

/// Reference J_v = d_{w'}(v, a) by Dijkstra (Lemma 3.2) on G' from the root. For tests and
/// checkers only; not a step of the paper (Alg 7 computes J).
template <class W>
std::vector<W> owner_coordinates_ref(const AuxGraph<W>& aux);

/// Checks the construction and Lemmas 7.1-7.2 on an output, with any coordinate vector J
/// (indexed by aux vertex): classification as in Sec 7.1, cut corridor edges macro, the macro
/// bounds, J_a = 0, r <= J <= R <= h r, J_x <= J_y + w' and 2 w' <= 3 w on every edge of G', T'
/// the owner spokes rooted at a with depth <= h and the root bottleneck / distance r, R, and the
/// MST cycle property of T' in G'. Cost O(m h + N). Returns one message per violation.
template <class W>
std::vector<std::string> check_aux_metric(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                          const OuterEdgeClassification<W>& cls,
                                          const AuxGraph<W>& aux, const std::vector<W>& J);

}  // namespace sssp::paper::owner
