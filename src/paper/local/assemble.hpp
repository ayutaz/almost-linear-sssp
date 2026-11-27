#pragma once
// Assembling the local hierarchy (Sec 5.4, txt l.1182-1243; Alg 4 l.10-12) and the
// post-processing in the proof of Lemma 5.10 (l.1265-1300).
//
// Input: the owner decomposition (Sec 5.1), the local radial state representation (Sec 5.2) and
// its threshold joins (Sec 5.3). Output: the ordered hierarchy on the states with physical
// connectors, connection constant c0 = 2 (Lemma 5.10).
//
// Steps (paper order):
//   1. Block scales (l.1183-1189). A singleton block keeps only its finite attachment join. For
//      a block U with >= 2 states, s_U is the lower boundary a K^i of the lowest native shell i
//      present (a = least coordinate of the near component; outward window 0 holds shell 0 only,
//      Q-33, so its s_U is a). Then s_U <= r_x < K^2 s_U on U (Q-34, local definition).
//   2. Attachment (l.1201-1203): z_U minimises theta over U; sigma_U is its named corridor key
//      (ThresholdJoins::b0_key). Infinite minimum: U has no finite attachment.
//   3. Alignment (l.1206-1214), c = 10: U is aligned iff s_U >= sigma_U / N^c, decided per
//      named key by successive medians (BFPRT, Lemma 3.1), one test s * N^c >= sigma per median.
//      Unaligned blocks and blocks without attachment use a separate grid from s_U.
//   4. Supports (l.1214-1224, Q-37): one interval per corridor of T_A, from
//      min(sigma, aligned s_U assigned to it) to N^16 sigma; every finite theta(z) lies in the
//      interval of its own key (checked). Intervals are sorted, overlaps merged, and a doubling
//      grid starts at the minimum of each merged component; a class that would cross the next
//      component start stops there (l.1279). Each key records the largest grid point <= its
//      lower end and the first grid point above its upper end (search sentinel).
//   5. t_U (l.1232-1236): aligned: the grid lower boundary of the class of s_U; separate grid:
//      t_U = s_U. Hence t_U in [s_U/2, s_U] (Q-36, local closed interval).
//   6. Y components (l.1189-1200): the components of the block-MST edges lighter than t_U. Each
//      gets its exact Kruskal merge tree (its < b edges sorted, equal weights merged
//      simultaneously; node height = weight, links of weight = height).
//   7. Separate grids (l.1209-1213, l.1236): classes [s_U 2^j, s_U 2^(j+1)) for the remaining
//      block-MST edges of an unaligned or unattached block.
//   8. Global classes (l.1237-1238, l.1244-1250): all other joins (owner joins, corridor joins,
//      block-MST edges of aligned blocks at least t_U) classified on the grid of their named key
//      by binary search; one simultaneous merge per class [g, g') at height g. Merges use the
//      disjoint sets of Lemma 3.6 on the join vertices (states and auxiliary pivot vertices).
//      Empty auxiliary nodes are discarded and unary nodes compressed on the fly (l.1281).
//   9. Links (l.1282-1292): every state arc is assigned to the pair of children it joins at the
//      lowest common ancestor of its endpoints (Lemma 3.9 LCA and level ancestor); at each node
//      the assigned arcs shorter than the class upper bound (at most the height for an exact
//      merge) form the child adjacency, and a DFS takes one discovery arc per child as link.
//  10. Order (l.1292-1293): children are topologically sorted by the arcs shorter than the node
//      height; roots (the SCCs of the whole state graph) in condensation order (l.1280).
//
// Correctness invariant (proof of Lemma 5.10, l.1271-1275): the states below every node form
// one SCC of the state arcs shorter than its class upper bound (at most the height for an exact
// merge), and its children are the SCCs of the arcs shorter than its height. This relies on
// Lemma 5.9 and on the scale relations of l.1189-1213; where the latter are needed (theta >= s_U
// on a nonsingleton block, values inside their support, separate-grid classes ending below the
// first attachment class) they are checked and a violation throws std::logic_error.
//
// Numerics (Q-38): N^c * s, N^16 * sigma and the grid doublings saturate at
// WeightTraits<W>::infinity() (paper/local/saturating.hpp). A saturated support end merges the
// support with every later one [deviation candidate D-13: exact reals might keep them apart; this
// moves grid points only, every class still has upper <= 2 * lower]. [Supplement] Grid points
// above the largest finite join weight are not generated: no value is classified there, so no
// decision changes.
//
// Lengths are only copied, compared and added (Sec 2); the per-block power a K^i is cached per
// near component and formed by doublings.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/base_case.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// c of the alignment test s_U >= sigma_U / N^c (Sec 5.4 l.1206).
inline constexpr unsigned kLocalAlignExponent = 10;
/// Upper end N^16 sigma of a named key's support (Sec 5.4 l.1216).
inline constexpr unsigned kLocalSupportExponent = 16;

inline constexpr std::uint32_t kNoComponent = std::numeric_limits<std::uint32_t>::max();

/// Sec 5.4 l.1206: in the nontrivial range N > 2 L_N^4 we have N >= 32 and h = b + 1 <= N; the
/// constants c = 10 and N^16 depend on both. A caller with other parameters (Mode B, tests)
/// must meet them (Q-03, conservative variant k >= max(32, b + 1)).
constexpr bool local_assembly_applies(std::size_t n, std::size_t b) {
  return n >= 32 && b < n;
}

enum class BlockScaleKind : std::uint8_t {
  kSingleton,   ///< one state: only its attachment join (l.1183)
  kAligned,     ///< s_U >= sigma_U / N^c: global grid
  kUnaligned,   ///< finite attachment, s_U < sigma_U / N^c: separate grid from s_U
  kUnattached,  ///< theta infinite on U: separate grid from s_U (l.1213)
};

/// Which merge created a hierarchy node.
enum class NodeOrigin : std::uint8_t {
  kExactLocal,    ///< exact Kruskal tree of a Y component (height = weight, upper inclusive)
  kSeparateGrid,  ///< local doubling class [s_U 2^j, s_U 2^(j+1)) of one block
  kGlobalGrid,    ///< global doubling class [g, g')
};

template <class W>
struct BlockScale {
  BlockScaleKind kind = BlockScaleKind::kSingleton;
  std::size_t size = 0;                   ///< states of the block
  W s{};                                  ///< s_U (nonsingleton blocks)
  W t{};                                  ///< t_U (nonsingleton blocks)
  StateId z = kNoState;                   ///< z_U, the theta minimiser (finite attachment only)
  std::uint32_t sigma_key = kNoCorridor;  ///< named corridor key of B0(z_U)
};

/// Support [lower, upper] of one named corridor key and its grid indices: lower_index is the
/// largest grid point <= lower, upper_index the first grid point > upper (grid.size() if none).
template <class W>
struct GridSupport {
  W lower{};
  W upper{};
  std::size_t lower_index = 0;
  std::size_t upper_index = 0;
};

/// Output of the assembly: the Alg 4 result plus the scales and classes it used (for tests,
/// later stages and diagnostics). Per-node vectors are indexed by hierarchy NodeId (leaves
/// included, with unspecified values).
template <class W>
struct LocalAssembly {
  LocalHierarchyResult<W> result;
  std::vector<BlockScale<W>> blocks;       ///< per block of the state representation
  std::vector<std::uint32_t> y_component;  ///< per state: Y id, kNoComponent in singleton blocks
  std::vector<W> node_upper;               ///< class upper bound of each internal node
  std::vector<std::uint8_t> node_upper_inclusive;  ///< 1 iff the bound is inclusive (exact)
  std::vector<NodeOrigin> node_origin;
  std::vector<W> grid;                   ///< global class lower boundaries, increasing
  std::vector<W> grid_upper;             ///< matching class upper boundaries (<= 2 * lower)
  std::vector<GridSupport<W>> supports;  ///< per corridor (named key, Q-37)
};

/// Sec 5.4 + Lemma 5.10 post-processing. Throws std::invalid_argument if the three inputs do
/// not belong together or if !local_assembly_applies(N, b), and std::logic_error if a scale
/// relation the construction relies on fails (see the header comment).
/// Counters (if `diag` is non-null): local_assemble_calls, local_assemble_blocks,
/// local_assemble_aligned_blocks, local_assemble_unaligned_blocks,
/// local_assemble_unattached_blocks, local_assemble_y_components (>= 2 states),
/// local_assemble_grid_points, local_assemble_global_classes, local_assemble_nodes,
/// local_assemble_links, plus the dsu.* counters of the merges.
template <class W>
LocalAssembly<W> assemble_local_hierarchy(const OwnerDecomposition<W>& od,
                                          const RadialStateGraph<W>& rs,
                                          const ThresholdJoins<W>& joins,
                                          Diagnostics* diag = nullptr);

}  // namespace sssp::paper
