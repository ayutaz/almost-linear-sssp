// Sec 3.2 minimum spanning trees: GrowthPass (Algorithm 1, Lemma 3.4), the three-pass MST
// (Lemma 3.5), and a Kruskal reference. The lightest edge per vertex pair is in edge_pairs.cpp.
#include "paper/primitives/mst.hpp"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "paper/primitives/fib_heap.hpp"

namespace sssp::paper {

namespace {

constexpr Vertex kNoGrowth = kNoVertex;

struct QuotientArc {
  Vertex to;
  std::uint32_t edge;  ///< index into QuotientGraph::edges
};

/// Adjacency of a quotient graph (each non-loop edge appears at both endpoints). O(k' + l').
template <class W>
struct QuotientAdjacency {
  std::vector<std::size_t> offsets;
  std::vector<QuotientArc> arcs;

  explicit QuotientAdjacency(const QuotientGraph<W>& h) : offsets(h.num_vertices + 1, 0) {
    for (const auto& e : h.edges) {
      if (e.a == e.b) continue;  // a loop lies inside its quotient vertex
      ++offsets[e.a + 1];
      ++offsets[e.b + 1];
    }
    for (std::size_t v = 0; v < h.num_vertices; ++v) offsets[v + 1] += offsets[v];
    arcs.resize(offsets[h.num_vertices]);
    std::vector<std::size_t> cursor(offsets.begin(), offsets.end() - 1);
    for (std::size_t i = 0; i < h.edges.size(); ++i) {
      const auto& e = h.edges[i];
      if (e.a == e.b) continue;
      const auto idx = static_cast<std::uint32_t>(i);
      arcs[cursor[e.a]++] = QuotientArc{e.b, idx};
      arcs[cursor[e.b]++] = QuotientArc{e.a, idx};
    }
  }
};

template <class W>
void validate_quotient(const QuotientGraph<W>& h) {
  if (h.num_vertices >= static_cast<std::size_t>(kNoVertex)) {
    throw std::length_error("quotient graph: too many vertices");
  }
  if (h.edges.size() >= static_cast<std::size_t>(kNoEdge)) {
    throw std::length_error("quotient graph: too many edges");
  }
  for (const auto& e : h.edges) {
    if (e.a >= h.num_vertices || e.b >= h.num_vertices) {
      throw std::out_of_range("quotient graph: edge endpoint out of range");
    }
  }
}

}  // namespace

template <class W>
QuotientGraph<W> singleton_quotient(std::size_t n, std::span<const Edge<W>> edges) {
  validate_edges<W>(n, edges);
  QuotientGraph<W> q;
  q.num_vertices = n;
  q.edges.reserve(edges.size());
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u == e.v) continue;
    q.edges.push_back(QuotientEdge<W>{e.u, e.v, e.weight, static_cast<EdgeId>(i)});
  }
  return q;
}

template <class W>
GrowthPassResult<W> growth_pass(const QuotientGraph<W>& h, std::size_t limit) {
  if (limit == 0) throw std::invalid_argument("growth_pass: the limit B must be >= 1");
  validate_quotient(h);
  const std::size_t k = h.num_vertices;
  const QuotientAdjacency<W> adj(h);

  GrowthPassResult<W> out;
  GrowthPassStats& st = out.stats;
  st.limit = limit;
  st.vertices_in = k;
  st.edges_in = adj.arcs.size() / 2;

  // growth_of[v]: the growth that scanned v (Alg 1 l.3 "unvisited" = kNoGrowth).
  std::vector<Vertex> growth_of(k, kNoGrowth);
  // offer[x]: quotient edge named by x's heap entry (its lightest cut edge, Sec 3.2).
  std::vector<std::uint32_t> offer(k, 0);
  std::vector<std::uint32_t> selected_q;  // F as quotient-edge indices
  FibHeap<MstEdgeKey<W>> heap(k);         // one node pool per pass; discard() per growth
  [[maybe_unused]] std::uint64_t boundary_comparisons = 0;

  for (Vertex start = 0; start < k; ++start) {  // Alg 1 l.4
    if (growth_of[start] != kNoGrowth) continue;
    const auto gid = static_cast<Vertex>(st.growths++);
    heap.discard();  // Alg 1 l.5: fresh heap and membership identifiers
    Vertex u = start;
    while (true) {                                                         // Alg 1 l.6
      growth_of[u] = gid;                                                  // Alg 1 l.7
      for (std::size_t i = adj.offsets[u]; i < adj.offsets[u + 1]; ++i) {  // Alg 1 l.7-8
        const auto& arc = adj.arcs[i];
        const Vertex x = arc.to;
        if (growth_of[x] == gid) continue;  // inside this growth
        const auto& e = h.edges[arc.edge];
        const MstEdgeKey<W> key{e.weight, e.id};
        if (!heap.contains(x)) {
          heap.insert(x, key);
          offer[x] = arc.edge;
          ++st.heap_inserts;
        } else {
          SSSP_COUNT(++boundary_comparisons);
          if (key < heap.key(x)) {  // improve x's boundary edge
            heap.decrease(x, key);
            offer[x] = arc.edge;
            ++st.heap_decreases;
          }
        }
      }
      if (heap.empty()) {  // Alg 1 l.9-10
        ++st.stopped_empty;
        break;
      }
      if (heap.size() > limit) {  // Alg 1 l.9-10; size = distinct outside vertices
        ++st.stopped_large;
        break;
      }
      st.max_extract_heap_size = std::max(st.max_extract_heap_size, heap.size());
      const auto top = heap.extract_min();  // Alg 1 l.11
      ++st.extractions;
      const Vertex v = top.id;
      selected_q.push_back(offer[v]);
      out.selected.push_back(top.key.id);
      if (growth_of[v] != kNoGrowth) {  // Alg 1 l.12-13: v belongs to an earlier growth
        ++st.stopped_attached;
        break;
      }
      u = v;  // Alg 1 l.14
    }
  }
  st.comparisons = heap.comparisons() + boundary_comparisons;

  // Alg 1 l.15 / proof of Lemma 3.4: graph search in F, then one scan of all edge endpoints.
  std::vector<std::size_t> f_off(k + 1, 0);
  for (const auto qe : selected_q) {
    ++f_off[h.edges[qe].a + 1];
    ++f_off[h.edges[qe].b + 1];
  }
  for (std::size_t v = 0; v < k; ++v) f_off[v + 1] += f_off[v];
  std::vector<Vertex> f_adj(f_off[k]);
  {
    std::vector<std::size_t> cursor(f_off.begin(), f_off.end() - 1);
    for (const auto qe : selected_q) {
      const auto& e = h.edges[qe];
      f_adj[cursor[e.a]++] = e.b;
      f_adj[cursor[e.b]++] = e.a;
    }
  }
  out.component_of.assign(k, kNoVertex);
  std::vector<Vertex> stack;
  Vertex next_comp = 0;
  for (Vertex s = 0; s < k; ++s) {
    if (out.component_of[s] != kNoVertex) continue;
    out.component_of[s] = next_comp;
    stack.push_back(s);
    while (!stack.empty()) {
      const Vertex x = stack.back();
      stack.pop_back();
      for (std::size_t i = f_off[x]; i < f_off[x + 1]; ++i) {
        const Vertex y = f_adj[i];
        if (out.component_of[y] != kNoVertex) continue;
        out.component_of[y] = next_comp;
        stack.push_back(y);
      }
    }
    ++next_comp;
  }
  out.quotient.num_vertices = next_comp;
  for (const auto& e : h.edges) {
    const Vertex a = out.component_of[e.a];
    const Vertex b = out.component_of[e.b];
    if (a == b) continue;  // loop-free quotient; parallel edges are kept
    out.quotient.edges.push_back(QuotientEdge<W>{a, b, e.weight, e.id});
  }
  st.vertices_out = out.quotient.num_vertices;
  st.edges_out = out.quotient.edges.size();
  return out;
}

GrowthLimits three_pass_limits(std::size_t k) noexcept {
  static_assert(sizeof(std::size_t) == sizeof(std::uint64_t), "64-bit size_t expected");
  // L_k = ceil(log2(k + 2)): the least c with 2^c >= k + 2, i.e. 2^c >= k and 2^c - k >= 2,
  // which avoids forming k + 2 (Sec 2: integer word operations, no floating-point log).
  std::size_t lk = 0;
  while (lk < 64) {
    const std::uint64_t p = std::uint64_t{1} << lk;
    if (p >= k && p - k >= 2) break;
    ++lk;
  }
  std::size_t log_lk = 0;  // ceil(log2 L_k)
  while ((std::size_t{1} << log_lk) < lk) ++log_lk;
  return GrowthLimits{std::max<std::size_t>(2, log_lk), lk};
}

namespace {

void record_pass(Diagnostics& d, std::size_t index, const GrowthPassStats& st) {
  const std::string p = "mst.pass" + std::to_string(index) + ".";
  if (st.limit != kUnboundedGrowth) d.counters[p + "limit"] += st.limit;
  d.counters[p + "vertices_in"] += st.vertices_in;
  d.counters[p + "edges_in"] += st.edges_in;
  d.counters[p + "vertices_out"] += st.vertices_out;
  d.counters[p + "edges_out"] += st.edges_out;
  d.counters[p + "growths"] += st.growths;
  d.counters[p + "extractions"] += st.extractions;
  d.ops.vertex_scans += st.vertices_in;  // each quotient vertex is scanned once (Lemma 3.4)
  d.ops.edge_scans += 2 * st.edges_in;
  d.ops.comparisons += st.comparisons;
  d.ops.queue_pushes += st.heap_inserts;
  d.ops.decrease_keys += st.heap_decreases;
  d.ops.queue_pops += st.extractions;
}

/// Lemma 3.5: passes with B1, B2 (chosen from the original k) and +infinity, retaining all
/// selected indices across contractions.
template <class W>
SpanningForestResult three_passes(std::size_t n, std::span<const Edge<W>> edges, Diagnostics* diag,
                                  bool require_connected) {
  if (require_connected && n == 0) {
    throw std::invalid_argument("mst_three_pass: the graph must have k >= 1 vertices");
  }
  QuotientGraph<W> q = singleton_quotient<W>(n, edges);
  const GrowthLimits lim = three_pass_limits(n);
  const std::size_t limits[3] = {lim.b1, lim.b2, kUnboundedGrowth};
  SpanningForestResult out;
  std::vector<bool> chosen(edges.size(), false);
  for (const std::size_t limit : limits) {
    // "A singleton requires no pass": stop once no quotient edge is left (for a connected
    // input this happens exactly when the quotient is a single vertex).
    if (q.edges.empty()) break;
    GrowthPassResult<W> r = growth_pass(q, limit);
    for (const EdgeId id : r.selected) chosen[id] = true;
    out.passes.push_back(r.stats);
    q = std::move(r.quotient);
  }
  if (require_connected && q.num_vertices != 1) {
    throw std::invalid_argument("mst_three_pass: the graph is disconnected");
  }
  // Increasing indices by one scan of the edge list (no comparison sort).
  for (std::size_t i = 0; i < chosen.size(); ++i) {
    if (chosen[i]) out.edges.push_back(static_cast<EdgeId>(i));
  }
  if (diag != nullptr) {
    diag->counters["mst.passes"] += out.passes.size();
    for (std::size_t i = 0; i < out.passes.size(); ++i) record_pass(*diag, i + 1, out.passes[i]);
  }
  return out;
}

}  // namespace

template <class W>
SpanningForestResult mst_three_pass(std::size_t n, std::span<const Edge<W>> edges,
                                    Diagnostics* diag) {
  return three_passes<W>(n, edges, diag, true);
}

template <class W>
SpanningForestResult mst_three_pass(const CsrGraph<W>& g, Diagnostics* diag) {
  return three_passes<W>(g.num_vertices(), g.edges(), diag, true);
}

template <class W>
SpanningForestResult spanning_forest_three_pass(std::size_t n, std::span<const Edge<W>> edges,
                                                Diagnostics* diag) {
  return three_passes<W>(n, edges, diag, false);
}

template <class W>
SpanningForestResult spanning_forest_three_pass(const CsrGraph<W>& g, Diagnostics* diag) {
  return three_passes<W>(g.num_vertices(), g.edges(), diag, false);
}

template <class W>
std::vector<EdgeId> mst_kruskal_ref(std::size_t n, std::span<const Edge<W>> edges) {
  validate_edges<W>(n, edges);
  std::vector<EdgeId> order(edges.size());
  std::iota(order.begin(), order.end(), EdgeId{0});
  std::sort(order.begin(), order.end(), [&](EdgeId a, EdgeId b) {
    return MstEdgeKey<W>{edges[a].weight, a} < MstEdgeKey<W>{edges[b].weight, b};
  });
  std::vector<Vertex> parent(n);
  std::vector<std::uint32_t> size(n, 1);
  std::iota(parent.begin(), parent.end(), Vertex{0});
  auto find = [&](Vertex x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];  // path halving
      x = parent[x];
    }
    return x;
  };
  std::vector<EdgeId> forest;
  for (const EdgeId id : order) {
    Vertex a = find(edges[id].u);
    Vertex b = find(edges[id].v);
    if (a == b) continue;
    if (size[a] < size[b]) std::swap(a, b);
    parent[b] = a;
    size[a] += size[b];
    forest.push_back(id);
  }
  std::sort(forest.begin(), forest.end());
  return forest;
}

// Explicit instantiations for the two supported weight types.
template QuotientGraph<std::uint64_t> singleton_quotient<std::uint64_t>(
    std::size_t, std::span<const Edge<std::uint64_t>>);
template QuotientGraph<double> singleton_quotient<double>(std::size_t,
                                                          std::span<const Edge<double>>);
template GrowthPassResult<std::uint64_t> growth_pass<std::uint64_t>(
    const QuotientGraph<std::uint64_t>&, std::size_t);
template GrowthPassResult<double> growth_pass<double>(const QuotientGraph<double>&, std::size_t);
template std::vector<EdgeId> mst_kruskal_ref<std::uint64_t>(std::size_t,
                                                            std::span<const Edge<std::uint64_t>>);
template std::vector<EdgeId> mst_kruskal_ref<double>(std::size_t, std::span<const Edge<double>>);
template SpanningForestResult mst_three_pass<std::uint64_t>(std::size_t,
                                                            std::span<const Edge<std::uint64_t>>,
                                                            Diagnostics*);
template SpanningForestResult mst_three_pass<double>(std::size_t, std::span<const Edge<double>>,
                                                     Diagnostics*);
template SpanningForestResult mst_three_pass<std::uint64_t>(const CsrGraph<std::uint64_t>&,
                                                            Diagnostics*);
template SpanningForestResult mst_three_pass<double>(const CsrGraph<double>&, Diagnostics*);
template SpanningForestResult spanning_forest_three_pass<std::uint64_t>(
    std::size_t, std::span<const Edge<std::uint64_t>>, Diagnostics*);
template SpanningForestResult spanning_forest_three_pass<double>(std::size_t,
                                                                 std::span<const Edge<double>>,
                                                                 Diagnostics*);
template SpanningForestResult spanning_forest_three_pass<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, Diagnostics*);
template SpanningForestResult spanning_forest_three_pass<double>(const CsrGraph<double>&,
                                                                 Diagnostics*);

}  // namespace sssp::paper
