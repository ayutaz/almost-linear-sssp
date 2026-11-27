// Mode B entry point of Algorithm 4 (Sec 5.4 l.1244-1261), task P5-EXP-LOCAL. See
// local_hierarchy.hpp.

#include "experimental/paper_experimental/local_hierarchy.hpp"

#include <stdexcept>
#include <string>
#include <utility>

#include "experimental/paper_experimental/local_assembly.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "paper/primitives/search.hpp"

namespace sssp::paper_experimental {

namespace {

/// Size counters of Lemma 5.10, under the names paper::local_hierarchy uses.
template <class W>
void count_sizes(const paper::LocalHierarchyResult<W>& r, Diagnostics& diag) {
  SSSP_COUNT(diag.counters["local_hierarchy_states"] += r.hierarchy.num_states());
  SSSP_COUNT(diag.counters["local_hierarchy_arcs"] += r.hierarchy.graph.arcs.size());
  SSSP_COUNT(diag.counters["local_hierarchy_nodes"] += r.hierarchy.nodes.size());
  static_cast<void>(r);
  static_cast<void>(diag);
}

[[maybe_unused]] const char* variant_counter(const LocalParams& l) {
  if (l.assembly == LocalAssembly::kExactHeight) return "paper_experimental_local_exact_height";
  return l.alignment == LocalAlignment::kPaper ? "paper_experimental_local_paper_grid"
                                               : "paper_experimental_local_lambda_grid";
}

}  // namespace

template <class W>
LocalHierarchyRun<W> local_hierarchy(const CsrGraph<W>& g, const PaperParams& params) {
  LocalHierarchyRun<W> run;
  params.record(run.diag);  // the parameter set travels with the result
  const std::size_t k = g.num_vertices();
  const LocalParams& l = params.local();
  run.branch = params.local_branch(k);  // Alg 4 l.2 and the variant's preconditions
  SSSP_COUNT(++run.diag.counters["paper_experimental_local_calls"]);

  if (run.branch == LocalBranch::kUnsupported) {
    const char* why = "Lambda alignment (D-8) needs N >= 2";
    if (l.assembly == LocalAssembly::kExactHeight) {
      why = "exact-height assembly (D-11) needs N >= 2";  // the alignment is unused (Q-60)
    } else if (l.alignment == LocalAlignment::kPaper) {
      why = "Sec 5.4 formulas need N >= 32 and b + 1 <= N (l.1206)";
    }
    throw std::invalid_argument("paper_experimental::local_hierarchy: k = " + std::to_string(k) +
                                " passes the size test but the " + why +
                                "; Mode B treats this instance as a failure");
  }
  const paper::LocalHierarchyParams size = params.local_size_params(k);

  if (l.alignment == LocalAlignment::kPaper && l.assembly == LocalAssembly::kGrid) {
    // The paper's Sec 5.4: Mode A's entry point, base case (with its FallbackRecord) or radial.
    SSSP_COUNT(++run.diag.counters[run.branch == LocalBranch::kBaseCase
                                       ? "paper_experimental_local_base_case"
                                       : "paper_experimental_local_radial"]);
    if (run.branch == LocalBranch::kRadial) SSSP_COUNT(++run.diag.counters[variant_counter(l)]);
    run.local = paper::local_hierarchy(g, &run.diag, size);
    return run;
  }

  SSSP_COUNT(++run.diag.counters["local_hierarchy_calls"]);
  if (run.branch == LocalBranch::kBaseCase) {
    // Alg 4 l.3-5: the paper's base case, recorded as a fallback (Q-10).
    SSSP_COUNT(++run.diag.counters["paper_experimental_local_base_case"]);
    SSSP_COUNT(++run.diag.counters["local_hierarchy_base_case"]);
    run.local = paper::local_hierarchy_base_case(g, &run.diag);
    run.diag.fallbacks.push_back({paper::local_base_case_reason(size), k, g.num_edges()});
    count_sizes(run.local, run.diag);
    return run;
  }

  // Alg 4 l.6-12 with the D-8 or D-11 assembly.
  SSSP_COUNT(++run.diag.counters["paper_experimental_local_radial"]);
  SSSP_COUNT(++run.diag.counters[variant_counter(l)]);
  SSSP_COUNT(++run.diag.counters["local_hierarchy_radial"]);
  if (paper::connected_components(g).count != 1) {
    throw std::invalid_argument(
        "paper_experimental::local_hierarchy: Alg 4 requires a connected input (Sec 5.1 "
        "l.961); call it per connected component");
  }
  const std::size_t b = size.separator == 0 ? paper::local_separator_size(k) : size.separator;
  // l.6: MST, tree-dominated edges removed, owner cells (validates positive finite weights).
  const auto od = paper::owner_decomposition(g, b, &run.diag);
  // l.7-8: classification, windows and macro resets.
  const auto rs = paper::local_radial_states(g, od, &run.diag);
  // l.9-10 (first half): bottleneck minima and join records.
  const auto joins = paper::compute_joins(od, rs, &run.diag);
  // l.10 (second half) - l.12 with the Mode B alignment / assembly.
  run.local =
      std::move(assemble_local_hierarchy(od, rs, joins, l.alignment, l.assembly, &run.diag).result);
  count_sizes(run.local, run.diag);
  return run;
}

template LocalHierarchyRun<std::uint64_t> local_hierarchy<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, const PaperParams&);
template LocalHierarchyRun<double> local_hierarchy<double>(const CsrGraph<double>&,
                                                           const PaperParams&);

}  // namespace sssp::paper_experimental
