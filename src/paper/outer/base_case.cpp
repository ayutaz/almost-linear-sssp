#include "paper/outer/base_case.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/primitives/dijkstra.hpp"
#include "paper/primitives/edge_pairs.hpp"

namespace sssp::paper::outer {

namespace {

/// Alg 8 l.5, first half: one edge per unordered vertex pair, the (weight, id)-least of its
/// parallel edges (the Sec 3.2 edge order "compare edges by the pair (length, index)"). Loops
/// carry no path and are dropped. Edge i of `graph` is edge `input_edge[i]` of the input; kept
/// edges stay in increasing id order.
template <class W>
struct PairFilter {
  CsrGraph<W> graph;
  std::vector<EdgeId> input_edge;
};

/// Sec 2 (txt l.287): "a scan retaining a lightest edge per vertex pair", by the shared
/// primitive (edge_pairs.hpp: two stable counting-sort passes, Lemma 3.1, and one counted length
/// comparison per further edge of a pair). O(m + N) word operations plus at most m comparisons.
template <class W>
PairFilter<W> retain_lightest_pairs(const CsrGraph<W>& g, OpCounters& ops) {
  PairFilter<W> f;
  f.input_edge = lightest_edge_per_pair<W>(g.num_vertices(), g.edges(), &ops);
  std::vector<Edge<W>> edges;
  edges.reserve(f.input_edge.size());
  for (const EdgeId id : f.input_edge) edges.push_back(g.edge(id));
  f.graph = CsrGraph<W>::from_edges(g.num_vertices(), std::move(edges));
  return f;
}

enum class Search : std::uint8_t { kScan, kBinaryHeap };

template <class W>
SsspResult<W> base_case(const CsrGraph<W>& g, Vertex s, Search search) {
  const std::size_t n = g.num_vertices();
  if (s >= n) {
    throw std::out_of_range("outer base case: source " + std::to_string(s) + " out of range");
  }
  SsspResult<W> r;
  // Alg 8 l.4: the base branch is recorded before any work, with the size of the graph this
  // branch solves (CLAUDE.md fallback policy; ADR 0008).
  r.diag.fallbacks.push_back(FallbackRecord{
      std::string(search == Search::kScan ? kOuterBaseCaseScanReason : kOuterBaseCaseD1Reason), n,
      g.num_edges()});
  r.diag.counters["outer_base_cases"] += 1;

  // Alg 8 l.5: retain a lightest edge per vertex pair.
  const PairFilter<W> f = retain_lightest_pairs(g, r.diag.ops);
  r.diag.counters["outer_base_pair_edges_kept"] += f.graph.num_edges();
  r.diag.counters["outer_base_pair_edges_removed"] += g.num_edges() - f.graph.num_edges();

  // Alg 8 l.5: exact source distances by Dijkstra's search (Lemma 3.2) without a radius.
  SearchArena<W> arena;
  const W unbounded = WeightTraits<W>::infinity();
  if (search == Search::kScan) {
    dijkstra_scan_search(f.graph, s, unbounded, arena, &r.diag.ops);
  } else {
    dijkstra_search(f.graph, s, unbounded, arena, &r.diag.ops);  // deviation D-1
  }
  r = detail::export_result(n, arena, std::move(r));
  for (EdgeId& e : r.parent_edge) {
    if (e != kNoEdge) e = f.input_edge[e];  // back to the edge ids of g
  }
  return r;
}

}  // namespace

template <class W>
SsspResult<W> outer_base_case(const CsrGraph<W>& g, Vertex s) {
  return base_case(g, s, Search::kScan);
}

template <class W>
SsspResult<W> outer_base_case_d1(const CsrGraph<W>& g, Vertex s) {
  return base_case(g, s, Search::kBinaryHeap);
}

template SsspResult<std::uint64_t> outer_base_case(const CsrGraph<std::uint64_t>&, Vertex);
template SsspResult<double> outer_base_case(const CsrGraph<double>&, Vertex);
template SsspResult<std::uint64_t> outer_base_case_d1(const CsrGraph<std::uint64_t>&, Vertex);
template SsspResult<double> outer_base_case_d1(const CsrGraph<double>&, Vertex);

}  // namespace sssp::paper::outer
