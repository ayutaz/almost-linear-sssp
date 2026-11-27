#pragma once
// Mode B (paper-experimental) entry point of Algorithm 4, LocalHierarchy(H), task P5-EXP-LOCAL
// (arXiv 2611.34911v1, Sec 5.4 l.1244-1261, Lemma 5.10).
//
// Branch (Alg 4 l.2) decided by PaperParams::local_branch(k):
//   kBaseCase     Alg 4 l.3-5, the paper's base case (paper::local_hierarchy_base_case). Pushes
//                 FallbackRecord{paper::local_base_case_reason(params.local_size_params(k)),
//                 k, m} (CLAUDE.md fallback policy, Q-10).
//   kRadial       Alg 4 l.6-12: the shared Mode A stages owner_decomposition (Sec 5.1, b from
//                 the parameters, D-7), local_radial_states (Sec 5.2) and compute_joins
//                 (Sec 5.3), then the assembly chosen by params.local():
//                   (paper alignment, grid)  paper::local_hierarchy itself, so the output is
//                                            Mode A's (hierarchy, maps, fallbacks, counters);
//                   (Lambda, grid)           D-8, paper_experimental::assemble_local_hierarchy;
//                   (any, exact height)      D-11, the same engine (alignment unused, so
//                                            local_branch accepts every k >= 2 under either
//                                            alignment, Q-60).
//                 No FallbackRecord: every vertex is handled by the construction.
//   kUnsupported  the size test selects the radial branch but the variant's preconditions do
//                 not hold (Sec 5.4 l.1206 for the paper's formulas with the grid, N >= 2 for
//                 D-8 and for D-11). This is a failure: std::invalid_argument is thrown; no
//                 other construction is substituted (CLAUDE.md, Mode B).
//
// The parameter log (params.record) is written into the result's Diagnostics before anything
// else, so every hierarchy carries the parameter set it was built with.

#include <cstdint>

#include "experimental/paper_experimental/params.hpp"
#include "paper/local/base_case.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper_experimental {

/// Output of the Mode B Alg 4: the hierarchy with native entries and arc edges, the branch
/// taken and the diagnostics (parameter log, fallbacks, counters).
template <class W>
struct LocalHierarchyRun {
  paper::LocalHierarchyResult<W> local;
  LocalBranch branch = LocalBranch::kBaseCase;
  Diagnostics diag;
};

/// Algorithm 4 on `g` (k = g.num_vertices(), m = g.num_edges()) under `params`.
/// Throws std::invalid_argument for kUnsupported, for a disconnected input on the radial
/// branch (Sec 5.1 l.961) and for a weight that is not positive and finite (Sec 5.1 l.962); and
/// std::logic_error if a Sec 5.4 scale relation fails (it would be a bug or a counterexample to
/// D-8, never silently repaired).
/// Counters: those of the Mode A stages and of the assembly, the local_hierarchy_* counters of
/// paper::local_hierarchy, and paper_experimental_local_{calls,base_case,radial,paper_grid,
/// lambda_grid,exact_height}.
template <class W>
LocalHierarchyRun<W> local_hierarchy(const CsrGraph<W>& g, const PaperParams& params);

}  // namespace sssp::paper_experimental
