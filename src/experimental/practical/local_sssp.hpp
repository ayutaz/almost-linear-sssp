#pragma once
// Mode C (practical-optimized) single-source shortest paths through ONE local hierarchy, task
// P5-PRACTICAL-LOCAL. Benchmark / report name: "practical-local", mode "practical-optimized".
//
// This is NOT the paper's algorithm (arXiv 2611.34911v1). Algorithm 8 builds the outer
// construction of Sec 6-8 and uses local hierarchies (Alg 4) only inside owner cells; here the
// whole contracted source component is handed to one Alg 4-style local hierarchy and evaluated
// once (Alg 3). The O(m + n log log log n) bound of Theorem 1 is not claimed: the radial state
// graph of Sec 5.2 on the whole component has no such size bound. The components are the
// paper's, called unchanged from src/paper/ and src/experimental/paper_experimental/.
//
// Steps of local_sssp(g, s, options):
//   1. Alg 8 l.2: contract zero-edge components, drop loops, keep the source component
//      (paper::outer::contract_zero_edges). N = contracted vertices, m' = contracted edges.
//   2. Mode C: not in the paper. Size dispatch on N (never on lengths):
//        N = 1                         nothing to construct; FallbackRecord
//                                      {kTrivialComponentReason, 1, 0}.
//        N <= direct_search_max_vertices   baseline binary-heap Dijkstra on the contracted
//                                      graph; FallbackRecord{direct_search_reason(t), N, m'}.
//        otherwise                     steps 3-5, no FallbackRecord.
//      The paper's thresholds (Sec 2: local k > 2 L_k^4, outer N > L_N^64) are not used: the
//      hierarchy is built for every N >= 2 above the direct-search limit.
//   3. Alg 4 l.6-10 with the separator b from the options (paper: b = L_k^4, deviation D-7):
//      paper::owner_decomposition (Sec 5.1), paper::local_radial_states (Sec 5.2),
//      paper::compute_joins (Sec 5.3), then paper_experimental::assemble_local_hierarchy with
//        Assembly::kExactHeight  the exact-height assembly, deviation D-11 (default), or
//        Assembly::kLambdaGrid   the Sec 5.4 grid with the Lambda alignment, deviation D-8.
//      Mode C: not in the paper. The connectivity search of the Alg 4 entry points is skipped
//      (step 1 returns a connected graph by construction; owner_decomposition still rejects a
//      disconnected input), and the stage outputs are released as soon as the next stage has
//      consumed them.
//   4. Sec 4.1-4.3: paper::normalize and paper::order_children run ONCE and
//      paper::evaluate_normalized runs on that result with the single initial offer
//      (native entry of the source, 0). Mode C: not in the paper as a separate step, it only
//      avoids the copy of the hierarchy into a SeededOrderedHierarchy that paper::evaluate
//      takes; Def 4.1 / 4.2 are not re-validated (validate_seeded_ordered_hierarchy is a test
//      tool). The label of a contracted vertex is the least emission over its states
//      (Lemma 5.10 with Lemma 4.6). A vertex without a finite label raises std::logic_error:
//      it is a failure, never repaired by another search.
//      supplemental_offers = false is deviation D-14 (Alg 2 l.6 skipped; Lemma 4.6 does not
//      use them, the visit bound of Lemma 4.7 is lost).
//   5. Alg 8 l.6, Lemma 8.5: paper::outer::restore_predecessor_tree gives the original
//      distances and the first-discovery tree of tight original arcs. It throws
//      std::invalid_argument if a finite label is not attained by a tight arc, so labels that
//      are not exact distances are reported, not returned.
//
// Option log. Every choice above is listed by option_log() and written by record_options() into
// Diagnostics.counters as "practical.<key>=<value>" = 1 (the convention of Mode B's
// "params.<key>=<value>"). Like FallbackRecords it is part of the result and is written also
// when SSSP_DIAGNOSTICS = 0, so a timing from the release build carries the choices it was
// measured with. An option is not a fallback: record_options pushes no FallbackRecord.
//   key                          values
//   assembly                     exact-height(D-11) | lambda-grid(D-8)
//   separator                    L_N | <b>
//   direct_search_max_vertices   <t>
//   supplemental_offers          on | off(D-14)
//   thresholds                   none(N>=2)            fixed
//   normalization                once-reused           fixed
//   validation                   skipped               fixed
//   connectivity_check           skipped               fixed
//   intermediates                released-per-stage    fixed
//   tree                         lemma-8.5-restore     fixed
// Not done by this version (so not logged): a contiguous layout of the hierarchy and of
// Evaluate's work arrays, and fewer allocations inside the stages; both need a modified copy
// of the components.
//
// Numerics. Integer weights: exact. double: exact whenever every sum the components form is
// exactly representable (e.g. dyadic weights of bounded size). For general doubles Evaluate's
// gaps and supplemental offers are only certified in exact arithmetic (evaluate.hpp, D-19 /
// D-15), and step 5 demands exactly tight arcs: an input whose rounding breaks either makes
// local_sssp throw (std::invalid_argument or std::logic_error); it never returns approximate
// labels. Sums of weights are also formed for masses and scaled thresholds, so for integer
// weights std::overflow_error can be raised on inputs where a Dijkstra search would not
// overflow (weights near 2^63; P3-EVAL-ARITH-FIX).
//
// Fallback share: FallbackRecord.vertices / edges count CONTRACTED vertices and edges (N, m'),
// as in Mode A (ADR 0008).

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::practical {

/// Name used by benchmark runners and reports.
inline constexpr std::string_view kLocalSsspName = "practical-local";
/// Mode of the benchmark schema (bench.schema.MODES).
inline constexpr std::string_view kLocalSsspMode = "practical-optimized";

/// How the hierarchy is assembled from the join records (step 3).
enum class Assembly : std::uint8_t {
  kExactHeight,  ///< deviation D-11: one node per distinct join weight, no grid (c0 = 1)
  kLambdaGrid,   ///< deviation D-8: Sec 5.4 doubling grid with the Lambda alignment (c0 = 2)
};

struct PracticalOptions {
  Assembly assembly = Assembly::kExactHeight;
  /// Separator size b of the owner cells (Def 5.1). 0 means L_N = ceil(log2(N + 2)) of the
  /// contracted size N; any other value is used as a constant b. (Paper: b = L_k^4, D-7.)
  std::uint64_t separator = 0;
  /// Mode C: not in the paper. A contracted source component with 2 <= N <= this many vertices
  /// is solved by binary-heap Dijkstra and recorded as a fallback; 0 disables the direct search.
  /// The default is a placeholder that has not been tuned by measurement.
  std::uint64_t direct_search_max_vertices = 16;
  /// Alg 2 l.6. false is deviation D-14.
  bool supplemental_offers = true;
};

/// One logged choice. `note` says what it is and whether the paper has it.
struct OptionEntry {
  std::string key;
  std::string value;
  std::string note;
};

/// Every choice of local_sssp under `options`, tunable and fixed (see the file comment).
std::vector<OptionEntry> option_log(const PracticalOptions& options);

/// Writes option_log(options) into diag.counters as "practical.<key>=<value>" = 1. Idempotent.
void record_options(const PracticalOptions& options, Diagnostics& diag);

/// FallbackRecord reason for a one-vertex contracted source component.
inline constexpr std::string_view kTrivialComponentReason =
    "practical-local: N = 1 (single contracted vertex, nothing to construct)";

/// FallbackRecord reason for the direct search with limit t:
/// "practical-local: N <= <t> (direct search, binary-heap Dijkstra)".
std::string direct_search_reason(std::uint64_t direct_search_max_vertices);

/// Distances from s (infinity if unreachable) and a shortest-path tree of original edges, by
/// the steps of the file comment. result.diag holds the option log, at most one FallbackRecord
/// (none iff the hierarchy handled the contracted component), and, when SSSP_DIAGNOSTICS = 1,
/// the counters practical_local_calls, practical_local_{trivial,direct_searches,hierarchies},
/// practical_local_states / _arcs / _nodes (size of the hierarchy before normalisation) and
/// those of every stage called (contract_*, owner / radial / join / local_assemble_*,
/// normalize_*, child_order_*, eval_*, restore_*).
/// Throws std::out_of_range if s is not a vertex; std::overflow_error, std::invalid_argument or
/// std::logic_error as described under "Numerics".
template <class W>
SsspResult<W> local_sssp(const CsrGraph<W>& g, Vertex s, const PracticalOptions& options = {});

extern template SsspResult<std::uint64_t> local_sssp(const CsrGraph<std::uint64_t>&, Vertex,
                                                     const PracticalOptions&);
extern template SsspResult<double> local_sssp(const CsrGraph<double>&, Vertex,
                                              const PracticalOptions&);

}  // namespace sssp::practical
