#pragma once
// Mode B (paper-experimental) entry to Algorithm 8, ShortestPaths(G, s) (arXiv 2611.34911v1,
// Sec 8.3, Alg 8 l.2-17, Theorem 1, Lemma 8.5).
//
// The steps are those of the Mode A entry (include/sssp/paper.hpp) and the code is shared
// (ADR 0003, Q-04: one implementation, parameterised): zero-edge contraction (l.2), the size
// test of l.4, the base case (l.5) or the outer construction (l.8-15,
// paper::outer::outer_construction) and the restoration of Lemma 8.5 (l.6, l.16-17). Only the
// time-only constants come from PaperParams, so that the outer branch runs on feasible sizes:
//   D-7   outer.threshold (N > c L_N^e) and outer.separator (b), local.threshold, local.separator
//   D-17  outer.sample_factor, residual_factor, p_log2 (Alg 5)
//   D-10  outer.retry_cap (Alg 5; not in the paper)
// The minimal forcing set (docs/research/04_outer_scales_pipeline.md Sec 3.4) is
// outer.threshold = 0 * L^0 with outer.separator = constant(b).
//
// Failure rule (CLAUDE.md, Mode B): an input whose correctness cannot be established is a
// failure. Nothing here substitutes a reference search for the construction: a parameter set
// the shared components cannot run is rejected with std::invalid_argument, and every exception
// of the construction or of the restoration (std::logic_error, std::overflow_error,
// std::domain_error; std::invalid_argument when a label is not attained by a tight original
// arc, Lemma 8.5) propagates to the caller. Not supported on the outer branch:
//   * local.alignment = kLambda (D-8) and local.assembly = kExactHeight (D-11): the atoms of
//     Alg 7 and of Sec 8.2 run the shared Alg 4 (paper::local_hierarchy) only;
//   * local.separator = L_k^e with e != 4: the shared Alg 4 takes one b for all atoms;
//   * outer.h_rule = kProofMinimum (D-9): H depends on the pivot count t, which is known only
//     after Alg 8 l.8.
//
// Fallback records: the base-case branch records exactly one FallbackRecord
// {paper::outer::kOuterBaseCaseScanReason, N, m'} as Mode A does. The outer branch records what
// its components record (Alg 4's base case per atom; NameRoutes' D-10 record) and nothing
// else. The parameter log is written by PaperParams::record ("params.*" counters), also when
// SSSP_DIAGNOSTICS = 0. Report the fallback share per reason
// (paper::outer::fallback_shares, Q-69) next to every timing.
//
// Lowered local thresholds (P5-EXP-THRESHOLDS, Q-03): the outer branch runs Alg 4 on atoms of
// sizes that depend on the input, so a local setting under which some size k selects the
// radial branch while the Sec 5.4 l.1206 preconditions fail (least_unsupported_local_size,
// thresholds.hpp) is rejected with std::invalid_argument before the construction starts, on
// every input. conservative_local_params(b) and q03_threshold_settings() give settings that
// pass this check.
//
// Branch record (written like the parameter log, also with SSSP_DIAGNOSTICS = 0):
//   mode_b.branch.outer_base_case = 1      Alg 8 l.5 was taken
//   mode_b.branch.outer_construction = 1   Alg 8 l.8-15 was taken; then also
//   mode_b.local.calls / .radial / .base_case
//       the Alg 4 calls on the atoms of Alg 7 (radial flag reported by owner_distances) and of
//       Sec 8.2 (Lemma 5.15; the same size test on the atom's state count), and how many took
//       each branch. base_case must equal the number of Alg 4 FallbackRecords; a mismatch
//       throws std::logic_error.

#include <cstdint>
#include <string_view>

#include "experimental/paper_experimental/params.hpp"
#include "paper/outer/shortest_paths.hpp"
#include "sssp/graph.hpp"
#include "sssp/paper.hpp"
#include "sssp/result.hpp"

namespace sssp::paper_experimental {

/// Algorithm name for benchmark runners and reports.
inline constexpr std::string_view kPaperExperimentalName = "paper-experimental";

/// The parameters of the shared outer construction for a contracted graph of n vertices.
/// Throws std::invalid_argument for a parameter set listed as unsupported in the file comment.
paper::outer::OuterConstructionParams outer_construction_params(const PaperParams& params,
                                                                std::uint64_t n);

/// Algorithm 8 on (g, s) under `params`: exact distances (infinity if unreachable) and a
/// shortest-path tree of original edges. `seed` initialises the generator of NameRoutes (Alg 8
/// l.9). result.diag holds the parameter log, the contraction and restoration counters, the
/// counters of the branch taken and its FallbackRecords.
/// Throws std::out_of_range if s is not a vertex, std::invalid_argument for an unsupported
/// parameter set when the outer branch is taken (including a local setting that is unsupported
/// at some size, see the file comment), and whatever the construction throws; it never
/// returns a result obtained by another algorithm in place of the outer construction.
template <class W>
SsspResult<W> paper_sssp(const CsrGraph<W>& g, Vertex s, const PaperParams& params,
                         std::uint64_t seed = paper::kDefaultPaperSeed);

extern template SsspResult<std::uint64_t> paper_sssp(const CsrGraph<std::uint64_t>&, Vertex,
                                                     const PaperParams&, std::uint64_t);
extern template SsspResult<double> paper_sssp(const CsrGraph<double>&, Vertex, const PaperParams&,
                                              std::uint64_t);

}  // namespace sssp::paper_experimental
