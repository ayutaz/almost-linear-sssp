#pragma once
// Internal steps of Algorithm 8 (arXiv 2611.34911v1, Sec 8.3) behind the public entry point
// include/sssp/paper.hpp. Exposed for tests only.

#include <cstdint>
#include <vector>

#include "paper/model/params.hpp"
#include "paper/outer/contract_restore.hpp"
#include "sssp/paper.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// Alg 8 l.4 negated: the outer construction runs iff N > L_N^64 (Sec 2, txt l.286), evaluated on
/// the fly with capped integer powers (no precomputed cut-off: the satisfying set is not
/// upward-closed, ADR 0008). False for every N < 2^64.
constexpr bool takes_outer_construction(std::uint64_t big_n) {
  return model::outer_nontrivial(big_n);
}

/// Alg 8 l.8-15 on the contracted graph with the paper's constants (b = L_N^64, the NameRoutes
/// constants of Alg 5, Alg 4 with 2 L_k^4): the exact contracted distances by the outer
/// construction, outer::outer_construction(c.graph, c.source, WordRng(seed),
/// outer::kPaperOuterParams, &diag). Never replaced by a direct search; it throws as that
/// function does. In Mode A it is unreachable through paper_sssp because
/// takes_outer_construction is false for every N < 2^64 (Q-04); called directly, b exceeds N and
/// the decomposition has one pivot (Q-72).
template <class W>
std::vector<W> outer_construction(const outer::ZeroContraction<W>& c, std::uint64_t seed,
                                  Diagnostics& diag);

extern template std::vector<std::uint64_t> outer_construction(
    const outer::ZeroContraction<std::uint64_t>&, std::uint64_t, Diagnostics&);
extern template std::vector<double> outer_construction(const outer::ZeroContraction<double>&,
                                                       std::uint64_t, Diagnostics&);

}  // namespace sssp::paper
