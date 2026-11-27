// Round-based Bellman–Ford: the independent correctness oracle.
//
// Deliberately shares no code with the Dijkstra family (no priority queue, no shared relax
// helper). Each round scans every vertex with a finite distance and relaxes all of its arcs;
// the algorithm stops after the first round without an improvement. With non-negative weights
// shortest paths are simple, so at most n-1 improving rounds occur.
//
// Floating point: relaxation uses strict `<` and fl(a + w) >= a for w >= 0, so the parent graph
// is always a tree (no cycle can be closed by a strict decrease), and the result is the exact
// minimum over paths of the left-to-right rounded path sums (see docs/design/numerics.md).
//
// Integer overflow: a candidate sum that would reach infinity() is skipped. If, after
// convergence, a vertex adjacent to a finite one is still at infinity, its true distance is not
// representable and std::overflow_error is thrown.

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "sssp/baseline.hpp"

namespace sssp::baseline {
namespace {

/// Returns false if a + b is not representable as a finite distance.
template <class W>
bool add_finite(W a, W b, W& out) noexcept {
  if constexpr (std::is_integral_v<W>) {
    if (a >= WeightTraits<W>::infinity() - b) {
      return false;
    }
    out = static_cast<W>(a + b);
    return true;
  } else {
    out = a + b;
    return out != WeightTraits<W>::infinity();
  }
}

}  // namespace

template <class W>
SsspResult<W> bellman_ford(const CsrGraph<W>& g, Vertex source) {
  const std::size_t n = g.num_vertices();
  if (source >= n) {
    throw std::out_of_range("bellman_ford: source " + std::to_string(source) +
                            " is not a vertex of a graph with " + std::to_string(n) + " vertices");
  }
  const W inf = WeightTraits<W>::infinity();
  SsspResult<W> r;
  r.dist.assign(n, inf);
  r.parent.assign(n, kNoVertex);
  r.parent_edge.assign(n, kNoEdge);
  r.dist[source] = WeightTraits<W>::zero();
  auto& ops = r.diag.ops;

  std::uint64_t rounds = 0;
  bool changed = true;
  while (changed) {
    if (rounds > n) {
      throw std::logic_error("bellman_ford: no convergence after n rounds");
    }
    changed = false;
    ++rounds;
    for (Vertex u = 0; u < n; ++u) {
      const W du = r.dist[u];
      if (du == inf) {
        continue;
      }
      SSSP_COUNT(++ops.vertex_scans);
      for (const auto& a : g.arcs(u)) {
        SSSP_COUNT(++ops.edge_scans);
        W cand{};
        SSSP_COUNT(++ops.additions);
        if (!add_finite(du, a.weight, cand)) {
          continue;
        }
        SSSP_COUNT(++ops.comparisons);
        if (cand < r.dist[a.to]) {
          r.dist[a.to] = cand;
          r.parent[a.to] = u;
          r.parent_edge[a.to] = a.edge;
          changed = true;
          SSSP_COUNT(++ops.relaxations);
        }
      }
    }
  }
  SSSP_COUNT(r.diag.counters["rounds"] = rounds);

  for (Vertex u = 0; u < n; ++u) {
    if (r.dist[u] == inf) {
      continue;
    }
    for (const auto& a : g.arcs(u)) {
      if (r.dist[a.to] == inf) {
        throw std::overflow_error("bellman_ford: shortest-path distance of vertex " +
                                  std::to_string(a.to) + " exceeds the representable range");
      }
    }
  }
  return r;
}

template SsspResult<std::uint64_t> bellman_ford(const CsrGraph<std::uint64_t>&, Vertex);
template SsspResult<double> bellman_ford(const CsrGraph<double>&, Vertex);

}  // namespace sssp::baseline
