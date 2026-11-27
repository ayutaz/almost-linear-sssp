#pragma once
// Threshold components and joins (Sec 5.3, txt l.1109-1180; Lemmas 5.8 and 5.9).
//
// Admissible paths use radial arcs and switches only (every arc kind except macro). For a
// state z (l.1110-1117):
//   I0(z)  least weight of a macro edge whose head reset can reach z,
//   O0(z)  least weight of a macro edge traversable after a path from z to one of its tail
//          copies, B0 = max(I0, O0);
//   I*(z), O*(z)  the same with the whole access path measured by its bottleneck (largest
//          physical arc weight, the macro edge included), theta = max(I*, O*).
// Empty minima are +infinity (WeightTraits<W>::infinity()).
//
// Computation (l.1112-1124): I0 and O0 by min propagation over the blocks in topological order
// of the macro-free block quotient (Lemma 5.7) and its reverse; within a block every state
// reaches every other, so both are constant on a block. I* in block order: the posted costs of a
// block (head-reset seeds and relaxed transfers) are spread by the multisource bottleneck on the
// block MST, min_x max(c_x, b_S(x, v)), with one postorder and one preorder sweep; each outgoing
// transfer is relaxed by max(value, length). O* is the same with all directions reversed. Block
// MSTs come from Lemma 3.5 (mst_three_pass, one call per block with >= 2 states).
//
// Join records (l.1124-1128): one auxiliary vertex per pivot; every block MST edge joins its
// endpoints at its physical weight; every state with finite theta joins its owner's auxiliary
// vertex at theta; auxiliary vertices are joined across the corridors (T_A edges) at the
// corridor weights. Lemma 5.9: for every W > 0, the components of the joins of weight <= W,
// restricted to states, are the SCCs of the state arcs of length <= W (and likewise with <).
//
// [Interpretation] Attaining macro edges are kept for the (weight, edge id)-least edge. B0's
// attaining edge is that of the larger of I0, O0, and of I0 when they are equal; its named key
// (b0_key) is the classification's macro key, i.e. the named corridor maximum (l.1022, used by
// Sec 5.4 as sigma_U).
//
// Only comparisons of lengths are performed (no additions): the bottleneck quantities are maxima
// and minima of arc lengths.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

enum class JoinKind : std::uint8_t { kBlockMst, kOwner, kCorridor };

/// An undirected weighted join between join vertices a and b. Join vertices 0..S-1 are the
/// states, S + p the auxiliary vertex of pivot p. `ref` is the arc id of the within-window arc
/// for kBlockMst, the pivot index for kOwner, and the corridor index for kCorridor.
template <class W>
struct JoinRecord {
  std::uint32_t a;
  std::uint32_t b;
  W weight;
  JoinKind kind;
  std::uint32_t ref;
};

template <class W>
struct ThresholdJoins {
  std::size_t num_states = 0;
  std::size_t num_pivots = 0;
  // Per state; +infinity is WeightTraits<W>::infinity().
  std::vector<W> i0, o0, b0;
  std::vector<W> i_star, o_star, theta;
  std::vector<EdgeId> i0_edge, o0_edge, b0_edge;  ///< attaining macro edge, kNoEdge if infinite
  std::vector<std::uint32_t> b0_key;              ///< named key of b0_edge, kNoCorridor if none
  /// Blocks in a topological order of the macro-free block quotient.
  std::vector<std::uint32_t> block_order;
  /// Join records: block MST edges, then owner joins, then corridor joins.
  std::vector<JoinRecord<W>> joins;

  [[nodiscard]] std::size_t num_join_vertices() const noexcept { return num_states + num_pivots; }
  [[nodiscard]] std::uint32_t aux_vertex(std::uint32_t pivot) const noexcept {
    return static_cast<std::uint32_t>(num_states + pivot);
  }
};

/// Sec 5.3 for the local radial state representation `rs` of the owner decomposition `od`
/// (Alg 4 l.9-10, first half). Throws std::invalid_argument if rs was built for a different
/// vertex set, and std::logic_error if the macro-free block quotient has a cycle (it is a DAG by
/// Lemma 5.7). Counters (if `diag` is non-null): threshold_joins_calls, threshold_joins_records,
/// threshold_joins_finite_theta; plus those of the block MSTs.
template <class W>
ThresholdJoins<W> compute_joins(const OwnerDecomposition<W>& od, const RadialStateGraph<W>& rs,
                                Diagnostics* diag = nullptr);

}  // namespace sssp::paper
