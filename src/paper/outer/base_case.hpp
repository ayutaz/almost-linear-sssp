#pragma once
// Outer base case of Algorithm 8 (arXiv 2611.34911v1, Sec 8.3, Alg 8 l.4-5; Sec 2, txt
// l.286-288; Lemma 3.2), Mode A.
//
// Alg 8 takes this branch when the contracted source component has N <= L_N^64 vertices, which
// holds for every N < 2^64 (src/paper/model/params.hpp). l.5: "retain a lightest edge per
// vertex pair and compute exact source distances by the scan-based Dijkstra search (Lemma 3.2)".
//
// Two variants (docs/adr/0008-mode-a-outer-base-case-variants.md):
//   * outer_base_case     - the scan version of Lemma 3.2, O(N^2 + m): paper-faithful.
//   * outer_base_case_d1  - the binary-heap version of Lemma 3.2, O((N + m) log N): deviation
//                           D-1 (same exact output, different running time only).
// Both always push one FallbackRecord (reason below, vertices = N, edges = m of the graph passed
// in, i.e. before the l.5 filter). It is recorded also when SSSP_DIAGNOSTICS = 0: a fallback is
// part of the result, not a diagnostic counter (CLAUDE.md, fallback policy).
//
// Lengths are only copied, compared and added (Sec 2 model); vertex ids and edge ids are
// integer indices.

#include <cstdint>
#include <string_view>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::outer {

/// FallbackRecord reason of the faithful base case (ADR 0008).
inline constexpr std::string_view kOuterBaseCaseScanReason =
    "outer: N <= L_N^64 (paper base case, scan Dijkstra)";
/// FallbackRecord reason of the binary-heap base case, deviation D-1 (ADR 0008).
inline constexpr std::string_view kOuterBaseCaseD1Reason =
    "outer: N <= L_N^64 (paper base case, binary-heap Dijkstra, deviation D-1)";

/// Alg 8 l.5 with the scan-based Dijkstra search of Lemma 3.2 on (g, s). g may have zero
/// weights, parallel edges and self-loops (Lemma 3.2 only needs non-negative lengths). Returns
/// exact distances (infinity if unreachable) and a shortest-path tree whose parent_edge ids are
/// edge ids of g. diag holds the operation counts, the counters outer_base_cases (= 1),
/// outer_base_pair_edges_kept and outer_base_pair_edges_removed (loops and non-lightest
/// parallel edges), and exactly one FallbackRecord{kOuterBaseCaseScanReason, N, m}.
/// Throws std::out_of_range if s is not a vertex and, for integer weights, std::overflow_error if
/// a distance would reach the infinity sentinel.
template <class W>
SsspResult<W> outer_base_case(const CsrGraph<W>& g, Vertex s);

/// As outer_base_case, with the binary-heap search of Lemma 3.2 (deviation D-1) and the
/// FallbackRecord reason kOuterBaseCaseD1Reason.
template <class W>
SsspResult<W> outer_base_case_d1(const CsrGraph<W>& g, Vertex s);

extern template SsspResult<std::uint64_t> outer_base_case(const CsrGraph<std::uint64_t>&, Vertex);
extern template SsspResult<double> outer_base_case(const CsrGraph<double>&, Vertex);
extern template SsspResult<std::uint64_t> outer_base_case_d1(const CsrGraph<std::uint64_t>&,
                                                             Vertex);
extern template SsspResult<double> outer_base_case_d1(const CsrGraph<double>&, Vertex);

}  // namespace sssp::paper::outer
