#pragma once
// Test support for the numerical scales of the outer construction (Sec 8.2 of arXiv
// 2611.34911v1, P3-OUTER-SCALES): runs Alg 8 l.8-12 up to the outer states (through
// outer_states_test_support.hpp) and Sec 8.2 on top, and provides the hand inputs of
// docs/research/04_outer_scales_pipeline.md Sec 3.5 (levels L3 and L4) and their variants.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "paper/outer/outer_states.hpp"
#include "paper/outer/scales.hpp"
#include "paper/outer/scales_numeric.hpp"
#include "paper/owner/aux_metric.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "unit/paper/outer/outer_states_test_support.hpp"

namespace sssp::test::outer {

/// F for a test case: from the exact h = b + 1 when b + 1 fits a word, otherwise b is the
/// saturated paper separator L_N^64 and F comes from the exact multiword b (Q-71).
template <class W>
paper::outer::OuterScaleFactor scale_factor_of(const OuterCase<W>& c) {
  const std::size_t b = c.od.separator;
  if (b == std::numeric_limits<std::size_t>::max()) {
    return paper::outer::paper_outer_scale_factor(c.cls.p_log2, c.n);
  }
  return paper::outer::outer_scale_factor(c.cls.p_log2, std::uint64_t{b} + 1,
                                          paper::outer::kOuterShellRatioLog2);
}

/// The inputs of Sec 8.2 and its output.
template <class W>
struct ScalesCase {
  OuterCase<W> c;
  paper::outer::OuterStates<W> os;
  paper::outer::OuterScaleFactor f;
  paper::outer::OuterScales<W> sc;
  Diagnostics diag;  ///< diagnostics of Sec 8.2 only
};

template <class W>
ScalesCase<W> run_scales(OuterCase<W> c) {
  ScalesCase<W> r;
  r.c = std::move(c);
  r.os = paper::outer::build_outer_states(r.c.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.scaled_distance);
  r.f = scale_factor_of(r.c);
  r.sc = paper::outer::build_outer_scales(r.c.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.calendars, r.os,
                                          r.f, &r.diag);
  return r;
}

/// An OuterCase from hand-given route names (one dictionary index per projected interowner
/// edge, in increasing input edge order), for weights whose sums do not fit the reference
/// all-pairs oracle of exact_names. P = 1.
template <class W>
OuterCase<W> make_outer_case_named(std::size_t n, const std::vector<Edge<W>>& edges, std::size_t b,
                                   std::vector<W> dictionary, std::vector<std::uint32_t> name) {
  OuterCase<W> c;
  c.n = n;
  c.edges = edges;
  c.g = CsrGraph<W>::from_edges(n, edges);
  c.od = paper::owner_decomposition(c.g, b);
  const auto pe = paper::owner::projected_interowner_edges(c.g, c.od);
  paper::pivots::RouteNames<W> names;
  names.p_log2 = 0;
  names.dictionary = std::move(dictionary);
  names.name = std::move(name);
  c.cls = paper::owner::classify_outer_edges(c.g, c.od, pe, names);
  c.aux = paper::owner::build_aux_graph(c.g, c.od, c.cls);
  c.jr = paper::owner::owner_distances(c.aux, &c.diag, {});
  return c;
}

/// Level L3 of note 04 Sec 3.5: the path 0 - 1 - 2 - 3 with weights (corridor, 3, 3) plus the
/// chord 1 - 3 of weight 5; b = 2. Edge ids: 0 = 0-1, 1 = 1-2, 2 = 2-3, 3 = 1-3.
/// By hand: T is the path (the chord closes a cycle of lighter edges); A = {0, 2}; the corridor
/// 0 .. 2 is cut at its heaviest edge 0-1, so p(1) = p(3) = 2 and r_1 = r_3 = 3; the chord is
/// intraowner (radial, d = 0) and survives (5 <= J_1 + J_3 = 6); 0-1 is interowner and macro.
template <class W>
std::vector<Edge<W>> level3_edges(W corridor) {
  return {{0, 1, corridor}, {1, 2, W{3}}, {2, 3, W{3}}, {1, 3, W{5}}};
}

/// Level L4 of note 04 Sec 3.5: a tree; b = 3. Edge ids: 0 = 0-1, 1 = 1-2, 2 = 0-3, 3 = 3-4,
/// 4 = 4-5. By hand: A = {0, 3}; p(1) = p(2) = 0 with r = 4, 4; p(4) = p(5) = 3 with r = 8, 8
/// (b_T(5, 3) = max(1, 8)); J = 4, 8 and 8, 9. The only interowner edge 0-3 joins two pivots.
template <class W>
std::vector<Edge<W>> level4_edges() {
  return {{0, 1, W{4}}, {1, 2, W{4}}, {0, 3, W{100}}, {3, 4, W{8}}, {4, 5, W{1}}};
}

}  // namespace sssp::test::outer
