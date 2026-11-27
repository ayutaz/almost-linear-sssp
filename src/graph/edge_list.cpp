// Edge-list utilities; contract in include/sssp/edge_list.hpp.
#include "sssp/edge_list.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace sssp {

namespace {

constexpr std::uint32_t kNoLabel = std::numeric_limits<std::uint32_t>::max();

/// Union-find with union by size and path halving over the vertices of one edge list.
class VertexSets {
 public:
  explicit VertexSets(std::size_t n) : parent_(n), size_(n, 1) {
    for (std::size_t v = 0; v < n; ++v) parent_[v] = static_cast<Vertex>(v);
  }

  Vertex find(Vertex x) {
    while (parent_[x] != x) {
      parent_[x] = parent_[parent_[x]];
      x = parent_[x];
    }
    return x;
  }

  void unite(Vertex a, Vertex b) {
    a = find(a);
    b = find(b);
    if (a == b) return;
    if (size_[a] < size_[b]) std::swap(a, b);
    parent_[b] = a;
    size_[a] += size_[b];
  }

 private:
  std::vector<Vertex> parent_;
  std::vector<std::size_t> size_;
};

/// Stable counting sort of the edge indices `order` by key(edge), keys in [0, n).
template <class Key>
std::vector<EdgeId> counting_sort(const std::vector<EdgeId>& order, std::size_t n, Key key) {
  std::vector<std::size_t> start(n + 1, 0);
  for (const EdgeId i : order) ++start[key(i) + 1];
  for (std::size_t k = 0; k < n; ++k) start[k + 1] += start[k];
  std::vector<EdgeId> out(order.size());
  for (const EdgeId i : order) out[start[key(i)]++] = i;
  return out;
}

/// Appends every edge whose endpoints both have a new id (new_id[x] != kNoVertex), renamed, in
/// input order, to `edges`, and its input index to `original_edge`.
template <class W>
void map_edges(const EdgeList<W>& g, std::span<const Vertex> new_id, std::vector<Edge<W>>& edges,
               std::vector<EdgeId>& original_edge) {
  for (std::size_t i = 0; i < g.edges.size(); ++i) {
    const Edge<W>& e = g.edges[i];
    if (new_id[e.u] == kNoVertex || new_id[e.v] == kNoVertex) continue;
    edges.push_back(Edge<W>{new_id[e.u], new_id[e.v], e.weight});
    original_edge.push_back(static_cast<EdgeId>(i));
  }
}

/// The subgraph induced by the vertices with in(v), renumbered in increasing order.
template <class W, class In>
InducedSubgraph<W> induced(const EdgeList<W>& g, In in) {
  InducedSubgraph<W> s;
  std::vector<Vertex> new_id(g.num_vertices, kNoVertex);
  for (std::size_t v = 0; v < g.num_vertices; ++v) {
    if (!in(static_cast<Vertex>(v))) continue;
    new_id[v] = static_cast<Vertex>(s.original_vertex.size());
    s.original_vertex.push_back(static_cast<Vertex>(v));
  }
  s.graph.num_vertices = s.original_vertex.size();
  map_edges(g, std::span<const Vertex>(new_id), s.graph.edges, s.original_edge);
  return s;
}

}  // namespace

template <class W>
FilteredEdgeList<W> keep_lightest_per_pair(const EdgeList<W>& g, LoopPolicy loops) {
  validate_edges<W>(g.num_vertices, g.edges);
  const std::size_t n = g.num_vertices;
  const std::size_t m = g.edges.size();
  auto lo = [&](EdgeId i) -> std::size_t {
    const Edge<W>& e = g.edges[i];
    return e.u < e.v ? e.u : e.v;
  };
  auto hi = [&](EdgeId i) -> std::size_t {
    const Edge<W>& e = g.edges[i];
    return e.u < e.v ? e.v : e.u;
  };

  std::vector<EdgeId> order;
  order.reserve(m);
  for (std::size_t i = 0; i < m; ++i) {
    const Edge<W>& e = g.edges[i];
    if (loops == LoopPolicy::kDrop && e.u == e.v) continue;
    order.push_back(static_cast<EdgeId>(i));
  }
  // Two stable passes (by hi, then by lo) group equal unordered pairs; within a group the
  // edges stay in increasing input index.
  order = counting_sort(order, n, hi);
  order = counting_sort(order, n, lo);

  std::vector<bool> keep(m, false);
  for (std::size_t a = 0; a < order.size();) {
    EdgeId best = order[a];
    std::size_t b = a + 1;
    for (; b < order.size() && lo(order[b]) == lo(best) && hi(order[b]) == hi(best); ++b) {
      // Strictly lighter only: an equal weight later in the group has a larger index.
      if (g.edges[order[b]].weight < g.edges[best].weight) best = order[b];
    }
    keep[best] = true;
    a = b;
  }

  FilteredEdgeList<W> r;
  r.graph.num_vertices = n;
  for (std::size_t i = 0; i < m; ++i) {
    if (!keep[i]) continue;
    r.graph.edges.push_back(g.edges[i]);
    r.original_edge.push_back(static_cast<EdgeId>(i));
  }
  return r;
}

template <class W>
FilteredEdgeList<W> relabel_vertices(const EdgeList<W>& g, std::span<const Vertex> new_id,
                                     std::size_t new_num_vertices) {
  validate_edges<W>(g.num_vertices, g.edges);
  if (new_id.size() != g.num_vertices) {
    throw std::invalid_argument("relabel_vertices: the map has " + std::to_string(new_id.size()) +
                                " entries for " + std::to_string(g.num_vertices) + " vertices");
  }
  if (new_num_vertices >= static_cast<std::size_t>(kNoVertex)) {
    throw std::length_error("relabel_vertices: too many vertices");
  }
  std::vector<bool> used(new_num_vertices, false);
  for (std::size_t v = 0; v < new_id.size(); ++v) {
    const Vertex x = new_id[v];
    if (x == kNoVertex) continue;
    if (x >= new_num_vertices) {
      throw std::out_of_range("relabel_vertices: vertex " + std::to_string(v) + " maps to " +
                              std::to_string(x) + " >= " + std::to_string(new_num_vertices));
    }
    if (used[x]) {
      throw std::invalid_argument("relabel_vertices: two vertices map to " + std::to_string(x));
    }
    used[x] = true;
  }

  FilteredEdgeList<W> r;
  r.graph.num_vertices = new_num_vertices;
  map_edges(g, new_id, r.graph.edges, r.original_edge);
  return r;
}

template <class W>
ComponentLabels connected_components(const EdgeList<W>& g) {
  validate_edges<W>(g.num_vertices, g.edges);
  const std::size_t n = g.num_vertices;
  VertexSets sets(n);
  for (const Edge<W>& e : g.edges) sets.unite(e.u, e.v);

  ComponentLabels c;
  c.label.assign(n, kNoLabel);
  std::vector<std::uint32_t> label_of_root(n, kNoLabel);
  // Scanning vertices in increasing order numbers components by their smallest vertex.
  for (std::size_t v = 0; v < n; ++v) {
    const Vertex root = sets.find(static_cast<Vertex>(v));
    if (label_of_root[root] == kNoLabel) {
      label_of_root[root] = static_cast<std::uint32_t>(c.count++);
      c.size.push_back(0);
    }
    c.label[v] = label_of_root[root];
    ++c.size[c.label[v]];
  }
  return c;
}

template <class W>
InducedSubgraph<W> extract_component(const EdgeList<W>& g, const ComponentLabels& components,
                                     std::uint32_t label) {
  validate_edges<W>(g.num_vertices, g.edges);
  if (components.label.size() != g.num_vertices) {
    throw std::invalid_argument("extract_component: labels for " +
                                std::to_string(components.label.size()) + " vertices, graph has " +
                                std::to_string(g.num_vertices));
  }
  if (label >= components.count) {
    throw std::out_of_range("extract_component: no component " + std::to_string(label));
  }
  return induced(g, [&](Vertex v) { return components.label[v] == label; });
}

template <class W>
InducedSubgraph<W> component_of(const EdgeList<W>& g, Vertex v) {
  if (v >= g.num_vertices) {
    throw std::out_of_range("component_of: vertex " + std::to_string(v) + " is not in the graph");
  }
  const ComponentLabels c = connected_components(g);
  return extract_component(g, c, c.label[v]);
}

template <class W>
InducedSubgraph<W> largest_component(const EdgeList<W>& g) {
  const ComponentLabels c = connected_components(g);
  if (c.count == 0) return InducedSubgraph<W>{};
  std::uint32_t best = 0;
  for (std::uint32_t l = 1; l < c.count; ++l) {
    if (c.size[l] > c.size[best]) best = l;  // strict: ties keep the smaller label
  }
  return extract_component(g, c, best);
}

#define SSSP_EDGE_LIST_INSTANTIATE(W)                                                           \
  template FilteredEdgeList<W> keep_lightest_per_pair<W>(const EdgeList<W>&, LoopPolicy);       \
  template FilteredEdgeList<W> relabel_vertices<W>(const EdgeList<W>&, std::span<const Vertex>, \
                                                   std::size_t);                                \
  template ComponentLabels connected_components<W>(const EdgeList<W>&);                         \
  template InducedSubgraph<W> extract_component<W>(const EdgeList<W>&, const ComponentLabels&,  \
                                                   std::uint32_t);                              \
  template InducedSubgraph<W> component_of<W>(const EdgeList<W>&, Vertex);                      \
  template InducedSubgraph<W> largest_component<W>(const EdgeList<W>&);

SSSP_EDGE_LIST_INSTANTIATE(std::uint64_t)
SSSP_EDGE_LIST_INSTANTIATE(double)

#undef SSSP_EDGE_LIST_INSTANTIATE

}  // namespace sssp
