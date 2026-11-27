#pragma once
// Shared pieces of the Dijkstra baselines. Every variant runs the same scan loop so that they
// differ only in the priority queue (see docs/design/baselines.md):
//   - every arc of a scanned vertex costs one addition (checked for integers) and one
//     comparison, including arcs back to already-settled vertices and self-loops;
//   - a vertex's parent / parent_edge are updated on every strict improvement only.

#include <stdexcept>
#include <string>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::baseline::detail {

/// Validates the source and returns a result with dist = infinity, no parents, dist[s] = 0.
template <class W>
SsspResult<W> init_result(const CsrGraph<W>& g, Vertex source) {
  const std::size_t n = g.num_vertices();
  if (source >= n) {
    throw std::out_of_range("source " + std::to_string(source) +
                            " is not a vertex (n = " + std::to_string(n) + ")");
  }
  SsspResult<W> r;
  r.dist.assign(n, WeightTraits<W>::infinity());
  r.parent.assign(n, kNoVertex);
  r.parent_edge.assign(n, kNoEdge);
  r.dist[source] = WeightTraits<W>::zero();
  return r;
}

/// Dijkstra over an indexed decrease-key heap (d-ary, Fibonacci, radix). The heap is sized to n
/// once; each vertex is pushed at most once and popped at most once.
template <class W, class Heap>
SsspResult<W> dijkstra_indexed(const CsrGraph<W>& g, Vertex source) {
  SsspResult<W> r = init_result(g, source);
  [[maybe_unused]] OpCounters& ops = r.diag.ops;
  Heap heap(g.num_vertices());
  heap.push(source, r.dist[source]);
  SSSP_COUNT(++ops.queue_pushes);
  while (!heap.empty()) {
    const auto [du, u] = heap.pop();
    SSSP_COUNT(++ops.queue_pops);
    SSSP_COUNT(++ops.vertex_scans);
    for (const auto& a : g.arcs(u)) {
      SSSP_COUNT(++ops.edge_scans);
      const W nd = WeightTraits<W>::checked_add(du, a.weight);
      SSSP_COUNT(++ops.additions);
      SSSP_COUNT(++ops.comparisons);
      if (nd < r.dist[a.to]) {
        r.dist[a.to] = nd;
        r.parent[a.to] = u;
        r.parent_edge[a.to] = a.edge;
        SSSP_COUNT(++ops.relaxations);
        if (heap.contains(a.to)) {
          heap.decrease_key(a.to, nd);
          SSSP_COUNT(++ops.decrease_keys);
        } else {
          heap.push(a.to, nd);
          SSSP_COUNT(++ops.queue_pushes);
        }
      }
    }
  }
  SSSP_COUNT(ops.comparisons += heap.comparisons());
  return r;
}

}  // namespace sssp::baseline::detail
