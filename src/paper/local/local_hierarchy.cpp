// Entry point of Algorithm 4 (Sec 5.4 l.1244-1261). See local_hierarchy.hpp.

#include "paper/local/local_hierarchy.hpp"

#include <stdexcept>
#include <string>

#include "paper/primitives/search.hpp"

namespace sssp::paper {

std::string local_base_case_reason(const LocalHierarchyParams& p) {
  return "local: k <= " + std::to_string(p.threshold_coeff) + "*L_k^" +
         std::to_string(p.threshold_exponent) + " (paper base case)";
}

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
    return out;
  }

  // Alg 4 l.6-12: the radial construction on a connected positive-weight input.
  if (diag != nullptr) SSSP_COUNT(++diag->counters["local_hierarchy_radial"]);
  if (connected_components(g).count != 1) {
    throw std::invalid_argument(
        "local_hierarchy: Alg 4 requires a connected input (Sec 5.1 l.961); call it per "
        "connected component");
  }
  throw std::logic_error(
      "local_hierarchy: the radial construction of Alg 4 l.6-12 needs the Sec 5.4 assembly, "
      "which is not implemented yet (P3-LOCAL-ASSEMBLE)");
}

template LocalHierarchyResult<std::uint64_t> local_hierarchy<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, Diagnostics*, const LocalHierarchyParams&);
template LocalHierarchyResult<double> local_hierarchy<double>(const CsrGraph<double>&, Diagnostics*,
                                                              const LocalHierarchyParams&);

}  // namespace sssp::paper
