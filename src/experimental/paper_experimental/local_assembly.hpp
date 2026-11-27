#pragma once
// Mode B (paper-experimental) assembly of the local hierarchy, task P5-EXP-LOCAL
// (arXiv 2611.34911v1, Sec 5.4 l.1182-1243, Alg 4 l.10-12, proof of Lemma 5.10 l.1265-1300).
//
// Inputs are the shared Mode A stages, called unchanged: owner_decomposition (Sec 5.1),
// local_radial_states (Sec 5.2) and compute_joins (Sec 5.3). Only the assembly differs, in two
// independent choices (params.hpp):
//
//   LocalAlignment (Sec 5.4 l.1206-1221)
//     kPaper   aligned iff s_U * N^10 >= sigma_U, support [., N^16 sigma]. Requires N >= 32 and
//              h = b + 1 <= N (l.1206), exactly as paper::assemble_local_hierarchy.
//     kLambda  deviation D-8 (research note 02 §4.3, Q-03 view (a)): aligned iff
//              s_U * Lambda >= sigma_U with Lambda = 2N * 3hK^2, support end
//              2 (3hK^2)^2 N^3 sigma. Requires N >= 2 only (any b >= 1). Why it is enough:
//                * an unaligned block has every within-window edge <= 3hK^2 s_U (Q-34) and at
//                  most N - 1 of them, so its block-MST length is < N 3hK^2 sigma / Lambda =
//                  sigma / 2 (the replaced l.1208). Each such edge is even < sigma / (2N), so
//                  the top separate-grid class [p, 2p) has p <= that edge and ends below
//                  sigma / N <= sigma / 2, while the first attachment class starts above
//                  sigma / 2 (theta(z_U) >= sigma_U by Lemma 5.2, and the grid of sigma's
//                  support doubles from a point <= sigma): the check of Lemma 5.10
//                  l.1276-1278 holds;
//                * theta <= 3hK^2 B0 (Lemma 5.8) and B0 <= N^3 sigma (macro test, Sec 5.2), so
//                  theta <= 3hK^2 N^3 sigma; an aligned block's edges are <= 3hK^2 s_U <=
//                  3hK^2 theta <= (3hK^2)^2 N^3 sigma; both lie in the support (the replaced
//                  l.1216-1221).
//              Lemmas 5.2-5.9 use N >= 2 and K >= 8h only (research note 02 §4.2). Only the
//              constants of the two tests change; every step and every check of the paper's
//              assembly is kept (a violated relation still throws std::logic_error).
//
//   LocalAssembly
//     kGrid         the paper's doubling grid with simultaneous merges per class (c0 = 2).
//     kExactHeight  deviation D-11 (research note 02 §3.4 [Supplement], Q-39): no grid; all
//                   join records sorted by (weight, index) and the records of one weight merged
//                   simultaneously at height = that weight. By Lemma 5.9 every node is an SCC of
//                   the arcs <= its height and its children are the SCCs of the arcs < it, so
//                   light arcs (< h) are acyclic between children and the arcs <= h connect
//                   them: links are at most the height (c0 = 1). The alignment is not used.
//                   Requires N >= 2. Extra cost: one O(J log J) comparison sort of the J join
//                   records (the grid's point is to avoid it, Lemma 5.10 time bound).
//
// The post-processing of Lemma 5.10 (links from the arcs at the lowest common ancestor, child
// order by the light arcs, roots in condensation order) is the same for every variant.
//
// Relation to Mode A: src/paper/local/assemble.cpp fixes c = 10 and N^16 and offers no way to
// pass other constants, so the engine below repeats its steps with the two tests made
// parameters (the Mode A file is not edited by this task). With (kPaper, kGrid) the engine
// reproduces paper::assemble_local_hierarchy field for field, counters included; tests check
// this on every family, which ties the D-8 path to Mode A. The Mode B entry point
// (local_hierarchy.hpp) does not use the engine for (kPaper, kGrid): it calls Mode A directly.
//
// Numerics (Q-38, D-13): s * Lambda and the support end are formed by doublings and integer
// multiplications that saturate at WeightTraits<W>::infinity() (paper/local/saturating.hpp).
//
// Output fields (paper::LocalAssembly): for kGrid as in Mode A. For kExactHeight `blocks`,
// `y_component`, `grid`, `grid_upper` and `supports` are empty, every internal node has
// node_upper = height, node_upper_inclusive = 1 and origin kExactLocal.

#include <cstddef>
#include <cstdint>

#include "experimental/paper_experimental/params.hpp"
#include "paper/local/assemble.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "sssp/result.hpp"

namespace sssp::paper_experimental {

/// D-8: s * Lambda, Lambda = 2N * 3hK^2 with h = b + 1 and K = 2^k_log2, saturating at the
/// type's infinity (also when b + 1 overflows).
template <class W>
W lambda_scaled(W s, std::uint64_t n, std::uint64_t b, unsigned k_log2);

/// D-8: the support end 2 (3hK^2)^2 N^3 sigma, saturating like lambda_scaled.
template <class W>
W lambda_support_end(W sigma, std::uint64_t n, std::uint64_t b, unsigned k_log2);

/// Whether the engine accepts N vertices with separator b under `alignment` and `assembly`:
/// (kPaper, kGrid) needs paper::local_assembly_applies(N, b); every other choice needs N >= 2.
bool local_assembly_variant_applies(std::size_t n, std::size_t b, LocalAlignment alignment,
                                    LocalAssembly assembly);

/// Sec 5.4 + Lemma 5.10 post-processing with the chosen alignment and assembly (see the file
/// comment). Throws std::invalid_argument if the inputs do not belong together or if
/// !local_assembly_variant_applies(N, b, alignment, assembly), and std::logic_error if a scale
/// relation of Sec 5.4 fails (grid only).
/// Counters (if `diag` is non-null): kGrid writes those of paper::assemble_local_hierarchy
/// (same names and meaning). kExactHeight writes local_assemble_calls, local_assemble_blocks,
/// local_assemble_nodes, local_assemble_links, local_assemble_exact_classes (distinct join
/// weights merged) and the dsu.* counters of the merges.
template <class W>
paper::LocalAssembly<W> assemble_local_hierarchy(const paper::OwnerDecomposition<W>& od,
                                                 const paper::RadialStateGraph<W>& rs,
                                                 const paper::ThresholdJoins<W>& joins,
                                                 LocalAlignment alignment, LocalAssembly assembly,
                                                 Diagnostics* diag = nullptr);

}  // namespace sssp::paper_experimental
