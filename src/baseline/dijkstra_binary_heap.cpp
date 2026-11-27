// Dijkstra with a lazy-deletion binary heap: no decrease-key; an improved vertex is pushed
// again and outdated entries are skipped when popped. O((n + m) log m) time, O(m) heap space.

#include <cstdint>

#include "baseline/dijkstra_common.hpp"
#include "baseline/heaps/binary_heap.hpp"
#include "sssp/baseline.hpp"

namespace sssp::baseline {

template <class W>
SsspResult<W> dijkstra_binary_heap(const CsrGraph<W>& g, Vertex source) {
  SsspResult<W> r = detail::init_result(g, source);
  [[maybe_unused]] OpCounters& ops = r.diag.ops;
  // At most one entry per successful relaxation plus the source: reserve once, never grow.
  heaps::BinaryHeap<W, Vertex> heap(g.num_arcs() + 1);
  heap.push(r.dist[source], source);
  SSSP_COUNT(++ops.queue_pushes);
  while (!heap.empty()) {
    const auto [du, u] = heap.pop();
    SSSP_COUNT(++ops.queue_pops);
    SSSP_COUNT(++ops.comparisons);
    if (r.dist[u] < du) continue;  // stale entry: u was improved after this push
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
        heap.push(nd, a.to);
        SSSP_COUNT(++ops.queue_pushes);
      }
    }
  }
  SSSP_COUNT(ops.comparisons += heap.comparisons());
  return r;
}

template SsspResult<std::uint64_t> dijkstra_binary_heap<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, Vertex);
template SsspResult<double> dijkstra_binary_heap<double>(const CsrGraph<double>&, Vertex);

}  // namespace sssp::baseline
