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

/// Alg 8 l.8-15 on the contracted graph: exact contracted distances by the outer construction.
/// NOT IMPLEMENTED (task P3-OUTER-PIPELINE): always throws std::logic_error("outer construction
/// not implemented"). It is never replaced by a direct search; in Mode A it is unreachable
/// because takes_outer_construction is false for every N < 2^64.
template <class W>
std::vector<W> outer_construction(const outer::ZeroContraction<W>& c, Diagnostics& diag);

extern template std::vector<std::uint64_t> outer_construction(
    const outer::ZeroContraction<std::uint64_t>&, Diagnostics&);
extern template std::vector<double> outer_construction(const outer::ZeroContraction<double>&,
                                                       Diagnostics&);

}  // namespace sssp::paper
