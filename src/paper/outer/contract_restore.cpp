#include "paper/outer/contract_restore.hpp"

#include <stdexcept>
#include <string>
#include <utility>

#include "paper/model/length.hpp"

namespace sssp::paper::outer {

namespace {

using model::Len;
using model::LenOps;

/// w == 0 as one counted comparison of lengths (Sec 2 model). An integer weight equal to the
/// infinity sentinel is a valid input weight; it is positive, hence never zero.
template <class W>
bool is_zero_weight(const LenOps<W>& ops, W w) {
  return ops.equal(Len<W>::from_weight(w), Len<W>::zero());
}

}  // namespace

template <class W>
ZeroContraction<W> contract_zero_edges(const CsrGraph<W>& g, Vertex source, Diagnostics* diag) {
  const std::size_t n = g.num_vertices();
  if (source >= n) throw std::out_of_range("contract_zero_edges: source is not a vertex");
  OpCounters* const sink = diag != nullptr ? &diag->ops : nullptr;
  const LenOps<W> ops(sink);

  ZeroContraction<W> c;
  c.vertex_map.assign(n, kNoVertex);
  c.zero_parent.assign(n, kNoVertex);
  c.zero_parent_edge.assign(n, kNoEdge);

  // Alg 8 l.2, one traversal of the source component (Theorem 1 proof, step 1: O(m + n)).
  // Roots are discovered over positive arcs; each unassigned root opens a new class, which is
  // filled by a BFS over zero arcs (its zero-weight spanning tree). Classes are therefore whole
  // zero-edge components, and the source class is opened first with s as representative.
  std::vector<Vertex> roots{source};
  std::vector<char> queued(n, 0);
  queued[source] = 1;
  std::vector<Vertex> members;
  for (std::size_t ri = 0; ri < roots.size(); ++ri) {
    const Vertex r = roots[ri];
    if (c.vertex_map[r] != kNoVertex) continue;  // absorbed by an earlier class via zero edges
    const auto id = static_cast<Vertex>(c.representative.size());
    c.representative.push_back(r);
    c.vertex_map[r] = id;
    members.assign(1, r);
    for (std::size_t qi = 0; qi < members.size(); ++qi) {
      const Vertex u = members[qi];
      SSSP_COUNT(if (sink != nullptr) { ++sink->vertex_scans; });
      for (const auto& a : g.arcs(u)) {
        SSSP_COUNT(if (sink != nullptr) { ++sink->edge_scans; });
        const Vertex v = a.to;
        if (is_zero_weight(ops, a.weight)) {
          if (c.vertex_map[v] == kNoVertex) {
            c.vertex_map[v] = id;
            c.zero_parent[v] = u;
            c.zero_parent_edge[v] = a.edge;
            members.push_back(v);
          }
        } else if (c.vertex_map[v] == kNoVertex && queued[v] == 0) {
          queued[v] = 1;
          roots.push_back(v);
        }
      }
    }
  }
  c.source = c.vertex_map[source];

  // Alg 8 l.2: map every edge of the source component to its classes, discard the loops (zero
  // edges, positive edges inside a class, original self-loops) and keep original edge ids.
  // Parallel edges stay; the lightest-per-pair filter is Alg 8 l.5, not l.2.
  std::vector<Edge<W>> edges;
  std::uint64_t dropped_loops = 0;
  std::uint64_t unreached_edges = 0;
  std::uint64_t zero_edges = 0;
  for (std::size_t i = 0; i < g.num_edges(); ++i) {
    const auto id = static_cast<EdgeId>(i);
    const auto& e = g.edge(id);
    const Vertex cu = c.vertex_map[e.u];
    if (cu == kNoVertex) {  // both endpoints lie outside the source component
      ++unreached_edges;
      continue;
    }
    const Vertex cv = c.vertex_map[e.v];
    if (cu == cv) {
      ++dropped_loops;
      if (diag != nullptr && is_zero_weight(ops, e.weight)) ++zero_edges;
      continue;
    }
    edges.push_back(Edge<W>{cu, cv, e.weight});
    c.original_edge.push_back(id);
  }
  const std::size_t big_n = c.representative.size();
  c.graph = CsrGraph<W>::from_edges(big_n, std::move(edges));

  if (diag != nullptr) {
    std::uint64_t reached = 0;
    for (const Vertex x : c.vertex_map) reached += x != kNoVertex ? 1 : 0;
    const std::uint64_t merged = reached - big_n;
    const std::uint64_t unreached = n - reached;
    auto& k = diag->counters;
    k["contract_input_vertices"] += n;
    k["contract_input_edges"] += g.num_edges();
    k["contract_retained_vertices"] += big_n;
    k["contract_retained_edges"] += c.graph.num_edges();
    k["contract_merged_vertices"] += merged;
    k["contract_unreached_vertices"] += unreached;
    k["contract_removed_vertices"] += merged + unreached;
    k["contract_dropped_loops"] += dropped_loops;
    k["contract_unreached_edges"] += unreached_edges;
    k["contract_removed_edges"] += dropped_loops + unreached_edges;
    k["contract_zero_edges"] += zero_edges;
  }
  return c;
}

template <class W>
SsspResult<W> restore_predecessor_tree(const CsrGraph<W>& g, Vertex source,
                                       const ZeroContraction<W>& c,
                                       const std::vector<W>& contracted_dist) {
  const std::size_t n = g.num_vertices();
  const std::size_t big_n = c.graph.num_vertices();
  if (source >= n || c.vertex_map.size() != n || c.source >= big_n ||
      c.vertex_map[source] != c.source || c.representative.size() != big_n ||
      c.representative[c.source] != source) {
    throw std::invalid_argument("restore_predecessor_tree: contraction does not match (g, s)");
  }
  if (contracted_dist.size() != big_n) {
    throw std::invalid_argument("restore_predecessor_tree: " + std::to_string(big_n) +
                                " contracted distances expected, got " +
                                std::to_string(contracted_dist.size()));
  }

  SsspResult<W> res;
  const LenOps<W> ops(&res.diag.ops);

  // Lemma 8.5: zero edges identify vertices at equal distance, so every member of a class gets
  // its class label; outside the source component the label is infinity (Alg 8 l.6).
  std::vector<Len<W>> label(n, Len<W>::infinity());
  for (std::size_t v = 0; v < n; ++v) {
    const Vertex x = c.vertex_map[v];
    if (x != kNoVertex) label[v] = Len<W>::from_weight(contracted_dist[x]);
  }
  if (!ops.equal(label[source], Len<W>::zero())) {
    throw std::invalid_argument("restore_predecessor_tree: source label is not 0");
  }

  // Lemma 8.5: traverse the finite tight original arcs u -> v, d(u) + w(uv) = d(v) < inf, from s
  // and keep the first-discovery arcs (a tree even with zero cycles). out[u] is the sum of
  // original weights down the tree; it equals d(u) for every discovered u by induction from the
  // root, so the tightness test and the output sum are the same addition.
  std::vector<Len<W>> out(n, Len<W>::infinity());
  res.dist.assign(n, WeightTraits<W>::infinity());
  res.parent.assign(n, kNoVertex);
  res.parent_edge.assign(n, kNoEdge);
  std::vector<char> found(n, 0);
  std::vector<Vertex> queue{source};
  found[source] = 1;
  out[source] = Len<W>::zero();
  for (std::size_t qi = 0; qi < queue.size(); ++qi) {
    const Vertex u = queue[qi];
    SSSP_COUNT(++res.diag.ops.vertex_scans);
    for (const auto& a : g.arcs(u)) {
      SSSP_COUNT(++res.diag.ops.edge_scans);
      const Vertex v = a.to;
      if (found[v] != 0 || label[v].is_infinite()) continue;
      const Len<W> w = Len<W>::from_weight(a.weight);
      if (w.is_infinite()) continue;  // integer weight at the sentinel: no finite label uses it
      Len<W> sum;
      try {
        sum = ops.add(out[u], w);
      } catch (const std::overflow_error&) {
        continue;  // the sum exceeds every finite length, so u -> v is not tight
      }
      if (!ops.equal(sum, label[v])) continue;
      found[v] = 1;
      out[v] = sum;
      res.parent[v] = u;
      res.parent_edge[v] = a.edge;
      queue.push_back(v);
    }
  }

  // Every vertex with a finite exact label has a tight shortest path from s (Lemma 8.5), so an
  // undiscovered finite label means the input distances were not exact.
  for (std::size_t v = 0; v < n; ++v) {
    if (found[v] != 0) {
      res.dist[v] = out[v].to_weight();
    } else if (label[v].is_finite()) {
      throw std::invalid_argument("restore_predecessor_tree: label of vertex " + std::to_string(v) +
                                  " is not attained by a tight arc (distances not exact)");
    }
  }
  res.diag.counters["restore_tree_vertices"] += queue.size();
  res.diag.counters["restore_unreachable_vertices"] += n - queue.size();
  return res;
}

template struct ZeroContraction<std::uint64_t>;
template struct ZeroContraction<double>;
template ZeroContraction<std::uint64_t> contract_zero_edges(const CsrGraph<std::uint64_t>&, Vertex,
                                                            Diagnostics*);
template ZeroContraction<double> contract_zero_edges(const CsrGraph<double>&, Vertex, Diagnostics*);
template SsspResult<std::uint64_t> restore_predecessor_tree(const CsrGraph<std::uint64_t>&, Vertex,
                                                            const ZeroContraction<std::uint64_t>&,
                                                            const std::vector<std::uint64_t>&);
template SsspResult<double> restore_predecessor_tree(const CsrGraph<double>&, Vertex,
                                                     const ZeroContraction<double>&,
                                                     const std::vector<double>&);

}  // namespace sssp::paper::outer
