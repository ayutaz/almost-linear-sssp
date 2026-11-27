#pragma once
// The local state graph of Sec 5.2 (txt l.1020-1107): Definition 5.4 (radial state
// representation), Lemmas 5.5-5.7.
//
// Def 5.4 has three parameters (l.1057-1060): a nonnegative vertex coordinate, the
// adjacent-shell ratio K, and the classification of edges into radial and macro. The generic
// builder build_radial_states takes exactly these, so Sec 8.1 (coordinate J_v, K = 16) can
// reuse it; the local instance (coordinate r_v, K = least power of two >= 8h, macro iff an
// intercell edge has w <= N^3 sigma(p, q)) is classify_local_edges + local_radial_states.
//
// Construction (l.1029-1055):
//   * Near components: radial edges with both coordinates positive and ratio <= K^2 (an
//     isolated positive vertex is its own component). With a = the least coordinate of the
//     component, the native shell i(v) satisfies a K^i <= r_v < a K^(i+1); it is found by a
//     search over near edges from the minimum, adjacent shells differing by at most two, with
//     cached powers a K^j (doublings only).
//   * Copies: inward window j >= 0 holds shells j, j+1; outward window j >= 0 holds shells
//     j-1, j. A positive vertex has the copies in_i (native inward), in_{i-1} (if i >= 1), out_i
//     (native outward) and out_{i+1}. A zero-coordinate vertex (a pivot) has in and out.
//   * Arcs, with the original physical length except switches (length 0):
//       within window  near edge with both endpoints in the window, both orientations
//       near transfer  inward copy at j -> native inward of a neighbour with shell < j;
//                      outward copy at j -> native outward of a neighbour with shell > j
//       far transfer   every inward copy of the higher end -> native inward of the lower end;
//                      every outward copy of the lower end -> native outward of the higher end
//       switch         every inward copy -> native outward (pivot: in -> out)
//       pivot          radial edge p-v with coordinate(p) = 0: every inward copy of v -> in(p),
//                      out(p) -> native outward of v
//       macro          each direction x -> y: every copy of x -> native inward of y
//   * Blocks: connected components of the within-window arcs.
//
// [Interpretation, Q-32] "point into its inward copy from positive inward states" (l.1052) is
// read as: from every inward copy of the positive endpoint, as in the far rule.
// [Interpretation, Q-33] Outward window 0 would hold shells -1 and 0; shell -1 does not exist,
// so it holds shell 0 only (the out_{i+1} copy of shell i is window i + 1 >= 1).
// [Supplement] A radial edge between two zero-coordinate vertices is rejected: the paper shows
// no positive radial edge joins two pivots (l.1025) and gives no arc rule for one.
//
// Numerics: K^2 r and a K^j are formed by doublings that saturate at WeightTraits<W>::infinity()
// (Q-38, paper/local/saturating.hpp); N^3 sigma by three integer multiplications through the
// binary expansion of N, saturating likewise.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

inline constexpr std::uint32_t kNoNearComponent = std::numeric_limits<std::uint32_t>::max();

/// Edge classification of Def 5.4. kOmitted edges (loops, edges outside the reduced graph)
/// produce no arcs.
enum class EdgeClass : std::uint8_t { kOmitted, kRadial, kMacro };

/// Per input edge: its class and, for macro edges, an opaque named key (the local construction
/// stores the index into OwnerDecomposition::corridors of the named corridor maximum, l.1022);
/// kNoCorridor otherwise.
struct EdgeClassification {
  std::vector<EdgeClass> edge_class;
  std::vector<std::uint32_t> macro_key;
};

enum class Phase : std::uint8_t { kInward, kOutward };

enum class RadialArcKind : std::uint8_t {
  kWithinWindow,
  kNearTransfer,
  kFarTransfer,
  kSwitch,
  kPivot,
  kMacro,
};

/// Definition 5.4. States are indexed 0..num_states-1; `graph` holds every arc.
template <class W>
struct RadialStateGraph {
  StateGraph<W> graph;
  std::vector<Vertex> projection;     ///< state -> physical vertex
  std::vector<Phase> phase;           ///< state -> inward / outward
  std::vector<std::uint32_t> window;  ///< state -> window index j (0 for zero coordinates)
  std::vector<std::uint32_t> block;   ///< state -> block id (0..num_blocks-1)
  std::size_t num_blocks = 0;
  std::vector<RadialArcKind> arc_kind;  ///< arc -> kind
  std::vector<EdgeId> arc_edge;         ///< arc -> input edge (kNoEdge for switches)

  // Per vertex.
  std::vector<StateId> native_inward;         ///< in_{i(v)}; the native entry (Def 5.4)
  std::vector<StateId> lower_inward;          ///< in_{i(v)-1}, kNoState if i(v) = 0 or coordinate 0
  std::vector<StateId> native_outward;        ///< out_{i(v)}
  std::vector<StateId> upper_outward;         ///< out_{i(v)+1}, kNoState for coordinate 0
  std::vector<std::uint32_t> shell;           ///< native shell i(v); 0 for coordinate 0
  std::vector<std::uint32_t> near_component;  ///< kNoNearComponent for coordinate 0
  std::vector<W> component_anchor;            ///< least coordinate a of each near component

  std::vector<W> coordinate;          ///< the vertex coordinate (r_v locally)
  EdgeClassification classification;  ///< the edge classes used
  unsigned shell_ratio_log2 = 0;      ///< K = 2^shell_ratio_log2

  [[nodiscard]] std::size_t num_states() const noexcept { return graph.num_states; }
  [[nodiscard]] std::size_t num_vertices() const noexcept { return native_inward.size(); }
};

/// Least k with 2^k >= 8h: K for the local construction (Sec 5.2 l.1029). h >= 1.
unsigned local_shell_ratio_log2(std::size_t hop_bound);

/// Def 5.4 on `g` with the given coordinate (size N, nonnegative), classification (size m) and
/// K = 2^shell_ratio_log2 (>= 2). Throws std::invalid_argument on size mismatches, a zero
/// exponent, a negative/NaN coordinate, or a radial edge whose endpoints both have coordinate
/// 0. Counters (if `diag` is non-null): radial_states_calls, radial_states_states,
/// radial_states_arcs, radial_states_blocks, radial_states_near_components.
template <class W>
RadialStateGraph<W> build_radial_states(const CsrGraph<W>& g, const std::vector<W>& coordinate,
                                        const EdgeClassification& classification,
                                        unsigned shell_ratio_log2, Diagnostics* diag = nullptr);

/// The local classification of Sec 5.2 (l.1021-1028): an edge of the reduced graph that is not
/// a loop is macro iff its endpoints lie in different owner cells p, q and
/// w <= N^macro_exponent * sigma(p, q) (paper: exponent 3); its key is the named corridor
/// maximum od.corridor_max(p, q). Every other reduced nonloop edge is radial; loops and edges
/// outside the reduced graph are omitted. N^e sigma is computed once per corridor edge.
template <class W>
EdgeClassification classify_local_edges(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                        unsigned macro_exponent = 3);

/// Alg 4 l.7-8: classify_local_edges(g, od) and build_radial_states with coordinate r_v and
/// K = 2^local_shell_ratio_log2(od.hop_bound).
template <class W>
RadialStateGraph<W> local_radial_states(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                        Diagnostics* diag = nullptr);

}  // namespace sssp::paper
