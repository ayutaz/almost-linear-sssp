// Entry point of Algorithm 4 (Sec 5.4 l.1244-1261). See local_hierarchy.hpp.

#include "paper/local/local_hierarchy.hpp"

#include <stdexcept>
#include <string>
#include <utility>

#include "paper/local/assemble.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "paper/primitives/search.hpp"

namespace sssp::paper {

std::string local_base_case_reason(const LocalHierarchyParams& p) {
  return "local: k <= " + std::to_string(p.threshold_coeff) + "*L_k^" +
         std::to_string(p.threshold_exponent) + " (paper base case)";
}

namespace {

/// Size counters of Lemma 5.10 (O(N) nodes, O(m + N) arcs).
template <class W>
void count_sizes(const LocalHierarchyResult<W>& r, Diagnostics* diag) {
  if (diag == nullptr) return;
  SSSP_COUNT(diag->counters["local_hierarchy_states"] += r.hierarchy.num_states());
  SSSP_COUNT(diag->counters["local_hierarchy_arcs"] += r.hierarchy.graph.arcs.size());
  SSSP_COUNT(diag->counters["local_hierarchy_nodes"] += r.hierarchy.nodes.size());
  static_cast<void>(r);
}

}  // namespace

template <class W>
LocalHierarchyResult<W> local_hierarchy(const CsrGraph<W>& g, Diagnostics* diag,
                                        const LocalHierarchyParams& p) {
  const std::size_t k = g.num_vertices();
  if (diag != nullptr) SSSP_COUNT(++diag->counters["local_hierarchy_calls"]);

  // Alg 4 l.2: the size test (integer arithmetic on k only).
  if (!local_construction_applies(k, p)) {
    // Alg 4 l.3-5: the paper's base case. Recorded as a fallback (CLAUDE.md policy, Q-10).
    auto out = local_hierarchy_base_case(g, diag);
    if (diag != nullptr) {
      SSSP_COUNT(++diag->counters["local_hierarchy_base_case"]);
      diag->fallbacks.push_back({local_base_case_reason(p), k, g.num_edges()});
    }
    count_sizes(out, diag);
    return out;
  }

  // Alg 4 l.6-12: the radial construction on a connected positive-weight input.
  if (diag != nullptr) SSSP_COUNT(++diag->counters["local_hierarchy_radial"]);
  if (connected_components(g).count != 1) {
    throw std::invalid_argument(
        "local_hierarchy: Alg 4 requires a connected input (Sec 5.1 l.961); call it per "
        "connected component");
  }
  const std::size_t b = p.separator == 0 ? local_separator_size(k) : p.separator;
  if (!local_assembly_applies(k, b)) {
    throw std::invalid_argument(
        "local_hierarchy: the Sec 5.4 assembly requires N >= 32 and h = b + 1 <= N (l.1206); got "
        "N = " +
        std::to_string(k) + ", b = " + std::to_string(b));
  }
  // l.6: MST (Lemma 3.5), tree-dominated edges removed, owner cells at separator size b.
  const auto od = owner_decomposition(g, b, diag);
  // l.7-8: radial / macro classification, windows and macro resets.
  const auto rs = local_radial_states(g, od, diag);
  // l.9-10 (first half): adjacency and bottleneck minima, join records.
  const auto joins = compute_joins(od, rs, diag);
  // l.10 (second half) - l.12: grids, merges, local trees, links and child order.
  auto out = std::move(assemble_local_hierarchy(od, rs, joins, diag).result);
  count_sizes(out, diag);
  return out;
}

template LocalHierarchyResult<std::uint64_t> local_hierarchy<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, Diagnostics*, const LocalHierarchyParams&);
template LocalHierarchyResult<double> local_hierarchy<double>(const CsrGraph<double>&, Diagnostics*,
                                                              const LocalHierarchyParams&);

}  // namespace sssp::paper
