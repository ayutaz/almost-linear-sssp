// Assembling the local hierarchy, Sec 5.4 and the proof of Lemma 5.10. See assemble.hpp.
//
// The steps live in the shared core (assembly_core.hpp, task P5-LOCAL-UNIFY); this file checks
// the precondition of l.1206 and passes the paper's constants: c = 10 (kLocalAlignExponent)
// and N^16 (kLocalSupportExponent).

#include "paper/local/assemble.hpp"

#include <stdexcept>
#include <string>

#include "paper/local/assembly_core.hpp"

namespace sssp::paper {

template <class W>
LocalAssembly<W> assemble_local_hierarchy(const OwnerDecomposition<W>& od,
                                          const RadialStateGraph<W>& rs,
                                          const ThresholdJoins<W>& joins, Diagnostics* diag) {
  const std::size_t n = od.num_vertices();
  const assembly::Rule rule = assembly::Rule::paper(n);  // c = 10 and N^16, with the grid
  assembly::check_inputs(od, rs, joins, rule.context);
  if (!local_assembly_applies(n, od.separator)) {
    throw std::invalid_argument(
        "assemble_local_hierarchy: Sec 5.4 l.1206 requires N >= 32 and h = b + 1 <= N (got N = " +
        std::to_string(n) + ", b = " + std::to_string(od.separator) + ")");
  }
  return assembly::assemble(od, rs, joins, rule, diag);
}

template LocalAssembly<std::uint64_t> assemble_local_hierarchy<std::uint64_t>(
    const OwnerDecomposition<std::uint64_t>&, const RadialStateGraph<std::uint64_t>&,
    const ThresholdJoins<std::uint64_t>&, Diagnostics*);
template LocalAssembly<double> assemble_local_hierarchy<double>(const OwnerDecomposition<double>&,
                                                                const RadialStateGraph<double>&,
                                                                const ThresholdJoins<double>&,
                                                                Diagnostics*);

}  // namespace sssp::paper
