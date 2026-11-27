#pragma once
// Entry point of Algorithm 4, LocalHierarchy(H) (Sec 5.4, txt l.1244-1261; Lemma 5.10).
//
//   l.2     if k <= 2 L_k^4: base case (Sec 2 l.286-290), see paper/local/base_case.hpp
//   l.6-12  otherwise: MST, owner cells with b = L_k^4 (Sec 5.1), radial states (Sec 5.2),
//           threshold joins (Sec 5.3) and the assembly of Sec 5.4.
//
// The size test is integer word arithmetic on k (never on lengths): L_k = ceil(log2(k + 2)) and
// c * L_k^e are evaluated with products capped at k (Sec 2 l.304), so no product wraps
// (paper/model/params.hpp). Mode A uses kPaperLocalParams (c = 2, e = 4). Following open
// question Q-04 there is one implementation: Mode B passes its own LocalHierarchyParams, and any
// value other than the paper's is a deviation to be listed in docs/PAPER_MAPPING.md.
//
// Fallback policy (CLAUDE.md): whenever the base case is taken, a FallbackRecord
// {"local: k <= 2*L_k^4 (paper base case)", k, m} is pushed (Q-10; with Mode B parameters the
// reason names the coefficient and exponent actually used). With Mode A constants every input
// with k <= 209952 takes this branch.
//
// Status: the radial branch (l.6-12) needs the Sec 5.4 assembly, which is not implemented yet.
// The branch throws std::logic_error instead of substituting the base case or any other
// construction, so a run that needs it can never be reported as the paper's algorithm.

#include <cstdint>
#include <string>

#include "paper/local/base_case.hpp"
#include "paper/model/params.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// Size test of Alg 4 l.2: the radial construction applies iff k > threshold_coeff *
/// L_k^threshold_exponent.
struct LocalHierarchyParams {
  std::uint64_t threshold_coeff = model::kLocalCoeff;   ///< paper: 2 (Sec 2 l.286)
  unsigned threshold_exponent = model::kLocalExponent;  ///< paper: 4 (Sec 2 l.286)

  friend constexpr bool operator==(const LocalHierarchyParams&,
                                   const LocalHierarchyParams&) = default;
};

/// The paper's values; Mode A never uses anything else.
inline constexpr LocalHierarchyParams kPaperLocalParams{};

/// k > c * L_k^e, evaluated exactly with capped products (no wrap-around for any k, c, e).
constexpr bool local_construction_applies(std::uint64_t k, const LocalHierarchyParams& p) {
  return model::exceeds_scaled_pow(k, p.threshold_coeff, model::L(k), p.threshold_exponent);
}

/// The FallbackRecord reason for the base-case branch under `p`:
/// "local: k <= <c>*L_k^<e> (paper base case)"; for the paper's values this is
/// "local: k <= 2*L_k^4 (paper base case)".
std::string local_base_case_reason(const LocalHierarchyParams& p);

/// Algorithm 4 on `g` (k = g.num_vertices(), m = g.num_edges()).
/// * k <= c L_k^e: returns local_hierarchy_base_case(g) and pushes
///   FallbackRecord{local_base_case_reason(p), k, m} into diag->fallbacks.
/// * otherwise: the radial construction. Requires a connected input (Sec 5.1 l.961; throws
///   std::invalid_argument otherwise). Not implemented yet: throws std::logic_error.
/// Throws std::invalid_argument for a zero (or, for integers, infinity-sentinel) weight: Alg 4
/// requires positive lengths and zero components are contracted by the caller (Sec 5.1 l.962).
/// Counters (if `diag` is non-null): local_hierarchy_calls, local_hierarchy_base_case,
/// local_hierarchy_radial, plus those of the branch taken.
template <class W>
LocalHierarchyResult<W> local_hierarchy(const CsrGraph<W>& g, Diagnostics* diag = nullptr,
                                        const LocalHierarchyParams& p = kPaperLocalParams);

}  // namespace sssp::paper
