#include "paper/outer/base_case.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/primitives/dijkstra.hpp"

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

/// Stable counting sort of `ids` by key(id) in [0, n), one pass per count, prefix sum and output
/// write (Lemma 3.1 proof). Integer work on vertex ids only.
template <class Key>
std::vector<EdgeId> counting_sort(const std::vector<EdgeId>& ids, std::size_t n, Key key) {
  std::vector<std::size_t> start(n + 1, 0);
  for (const EdgeId id : ids) ++start[key(id) + 1];
  for (std::size_t i = 0; i < n; ++i) start[i + 1] += start[i];
  std::vector<EdgeId> out(ids.size());
  for (const EdgeId id : ids) out[start[key(id)]++] = id;
  return out;
}

/// Sec 2 (txt l.287): "a scan retaining a lightest edge per vertex pair". Two stable
/// counting-sort passes (larger endpoint, then smaller endpoint) group the parallel edges of each
/// pair in increasing id order; one comparison per further group member keeps the least weight,
/// the earliest id on ties. O(m + N) word operations plus m comparisons of lengths.
template <class W>
PairFilter<W> lightest_edge_per_pair(const CsrGraph<W>& g,
                                     [[maybe_unused]] OpCounters& ops) {  // unused if counts off
  const std::size_t n = g.num_vertices();
  const std::size_t m = g.num_edges();
  auto lo = [&g](EdgeId id) {
    const auto& e = g.edge(id);
    return e.u < e.v ? e.u : e.v;
  };
  auto hi = [&g](EdgeId id) {
    const auto& e = g.edge(id);
    return e.u < e.v ? e.v : e.u;
  };
  std::vector<EdgeId> ids;
  ids.reserve(m);
  for (std::size_t i = 0; i < m; ++i) {
    const auto id = static_cast<EdgeId>(i);
    if (g.edge(id).u != g.edge(id).v) ids.push_back(id);  // loops are dropped
  }
  ids = counting_sort(counting_sort(ids, n, hi), n, lo);

  std::vector<char> keep(m, 0);
  for (std::size_t i = 0; i < ids.size();) {
    EdgeId best = ids[i];
    std::size_t j = i + 1;
    for (; j < ids.size() && lo(ids[j]) == lo(best) && hi(ids[j]) == hi(best); ++j) {
      SSSP_COUNT(++ops.comparisons);
      // Strict: an equal weight keeps the earlier (smaller) id.
      if (g.edge(ids[j]).weight < g.edge(best).weight) best = ids[j];
    }
    keep[best] = 1;
    i = j;
  }

  PairFilter<W> f;
  std::vector<Edge<W>> edges;
  for (std::size_t i = 0; i < m; ++i) {
    if (keep[i] == 0) continue;
    edges.push_back(g.edge(static_cast<EdgeId>(i)));
    f.input_edge.push_back(static_cast<EdgeId>(i));
  }
  f.graph = CsrGraph<W>::from_edges(n, std::move(edges));
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
  const PairFilter<W> f = lightest_edge_per_pair(g, r.diag.ops);
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
