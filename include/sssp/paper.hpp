#pragma once
// Public entry point of the Mode A (paper-faithful) single-source algorithm: Algorithm 8,
// ShortestPaths(G, s), of "Linear Single-Source Shortest Path, Almost" (arXiv 2611.34911v1,
// Sec 8.3, Theorem 1).
//
// What runs (docs/adr/0008-mode-a-outer-base-case-variants.md): zero-edge contraction (Alg 8
// l.2), the outer threshold N > L_N^64 on the contracted size N (l.3-4), the base case (l.5) and
// the predecessor-tree restoration of Lemma 8.5 (l.6). For every N < 2^64 the threshold is false,
// so every feasible input runs the base case; each run records exactly one FallbackRecord with
// the contracted size, also when SSSP_DIAGNOSTICS = 0. Timings of these variants are timings of
// the paper's preprocessing plus a Dijkstra search, not of the paper's outer construction.
//
// Variants (benchmark / CLI names):
//   "paper-faithful"     base case by the scan Dijkstra of Lemma 3.2, as Alg 8 l.5 states.
//   "paper-faithful-d1"  base case by the binary-heap Dijkstra of Lemma 3.2: deviation D-1.

#include <array>
#include <cstdint>
#include <optional>
#include <string_view>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

enum class PaperVariant : std::uint8_t {
  kFaithful,    ///< "paper-faithful": scan Dijkstra base case (Alg 8 l.5, Lemma 3.2)
  kFaithfulD1,  ///< "paper-faithful-d1": binary-heap Dijkstra base case (deviation D-1)
};

struct PaperOptions {
  PaperVariant variant = PaperVariant::kFaithful;
};

/// Every variant, in registration order.
inline constexpr std::array<PaperVariant, 2> kPaperVariants{PaperVariant::kFaithful,
                                                            PaperVariant::kFaithfulD1};

/// Algorithm name used by benchmark runners and reports.
constexpr std::string_view algorithm_name(PaperVariant v) noexcept {
  return v == PaperVariant::kFaithful ? "paper-faithful" : "paper-faithful-d1";
}

/// Inverse of algorithm_name; std::nullopt for any other string.
constexpr std::optional<PaperVariant> parse_algorithm_name(std::string_view name) noexcept {
  for (const PaperVariant v : kPaperVariants) {
    if (algorithm_name(v) == name) return v;
  }
  return std::nullopt;
}

/// Algorithm 8 on (g, s): exact distances (infinity if unreachable) and a shortest-path tree of
/// original edges. result.diag holds the contraction counters (contract_*), outer_base_cases,
/// outer_base_pair_edges_{kept,removed}, restore_* counters, operation counts, and one
/// FallbackRecord{reason of the variant (ADR 0008), N, m'} where N and m' are the vertex and edge
/// counts of the contracted source component. Throws std::out_of_range if s is not a vertex and,
/// for integer weights, std::overflow_error if a distance would reach the infinity sentinel.
///
/// Overflow rule (lead's reading, 2026-11-27), for this function and every component on its path:
///   * MAY throw std::overflow_error: a tentative path length d[u] + w of a reachable u reaches
///     the infinity sentinel, even if every true distance is representable (as the baselines
///     Bellman-Ford and Dijkstra do; e.g. the triangle 0-1:1, 1-2:1, 0-2:2^64 - 1).
///   * MUST NOT throw: only an auxiliary quantity that is not a path length reaches the sentinel,
///     such as the mass S(C) of Sec 4.1, upper bounds like the expiry D0 + S(C), slot ends or
///     other grid values. These saturate at infinity (deviation candidate D-23,
///     src/paper/hierarchy/evaluate.hpp); e.g. one edge of weight 2^63, where S(C) = 2^64.
/// Evaluate (Sec 4.3) is stricter than the first point: it skips a relaxation value beyond the
/// range and throws only if a state that received nothing else never emits. For double weights
/// Evaluate rejects a slot that rounds back to its start with std::domain_error (D-24). Neither
/// is on this function's path yet (base case only). Sec 4.2 order_children can still overflow on
/// sums of forest path lengths, an auxiliary quantity (open finding of P3-EVAL-ARITH-FIX).
template <class W>
SsspResult<W> paper_sssp(const CsrGraph<W>& g, Vertex s, PaperOptions options = {});

extern template SsspResult<std::uint64_t> paper_sssp(const CsrGraph<std::uint64_t>&, Vertex,
                                                     PaperOptions);
extern template SsspResult<double> paper_sssp(const CsrGraph<double>&, Vertex, PaperOptions);

}  // namespace sssp::paper
