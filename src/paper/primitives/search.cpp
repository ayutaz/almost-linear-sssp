// Sec 3.1, Lemma 3.1: searches and topological sorting in O(k + l).
#include "paper/primitives/search.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace sssp::paper {

namespace {

// Uniform adjacency access so that one search routine serves both graph types.
struct DigraphView {
  const Digraph& g;
  [[nodiscard]] std::size_t size() const noexcept { return g.num_vertices(); }
  [[nodiscard]] std::size_t degree(Vertex v) const noexcept { return g.heads(v).size(); }
  [[nodiscard]] std::pair<Vertex, EdgeId> arc(Vertex v, std::size_t i) const noexcept {
    return {g.heads(v)[i], g.arc_ids(v)[i]};
  }
};

template <class W>
struct CsrView {
  const CsrGraph<W>& g;
  [[nodiscard]] std::size_t size() const noexcept { return g.num_vertices(); }
  [[nodiscard]] std::size_t degree(Vertex v) const noexcept { return g.degree(v); }
  [[nodiscard]] std::pair<Vertex, EdgeId> arc(Vertex v, std::size_t i) const noexcept {
    const auto& a = g.arcs(v)[i];
    return {a.to, a.edge};
  }
};

void check_root(std::size_t k, Vertex root) {
  if (root >= k) {
    throw std::out_of_range("search root " + std::to_string(root) + " out of range");
  }
}

Traversal empty_traversal(std::size_t k) {
  Traversal t;
  t.parent.assign(k, kNoVertex);
  t.parent_edge.assign(k, kNoEdge);
  return t;
}

// Sec 3.1: mark at discovery; each vertex's adjacency list is scanned once, resumed from the
// saved position when the search returns to it (explicit stack instead of recursion).
template <class View>
Traversal dfs_impl(const View& g, Vertex root, OpCounters* ops) {
  const std::size_t k = g.size();
  check_root(k, root);
  Traversal t = empty_traversal(k);
  std::vector<char> discovered(k, 0);
  struct Frame {
    Vertex v;
    std::size_t next;
  };
  std::vector<Frame> stack;
  discovered[root] = 1;
  t.order.push_back(root);
  stack.push_back({root, 0});
  if (ops) SSSP_COUNT(++ops->vertex_scans);
  while (!stack.empty()) {
    const Vertex v = stack.back().v;
    const std::size_t i = stack.back().next;
    if (i == g.degree(v)) {
      stack.pop_back();
      continue;
    }
    ++stack.back().next;
    if (ops) SSSP_COUNT(++ops->edge_scans);
    const auto [x, id] = g.arc(v, i);
    if (discovered[x]) continue;
    discovered[x] = 1;
    t.parent[x] = v;
    t.parent_edge[x] = id;
    t.order.push_back(x);
    stack.push_back({x, 0});
    if (ops) SSSP_COUNT(++ops->vertex_scans);
  }
  return t;
}

// BFS: the discovery order itself serves as the FIFO queue.
template <class View>
Traversal bfs_impl(const View& g, Vertex root, OpCounters* ops) {
  const std::size_t k = g.size();
  check_root(k, root);
  Traversal t = empty_traversal(k);
  std::vector<char> discovered(k, 0);
  discovered[root] = 1;
  t.order.push_back(root);
  for (std::size_t head = 0; head < t.order.size(); ++head) {
    const Vertex v = t.order[head];
    if (ops) SSSP_COUNT(++ops->vertex_scans);
    for (std::size_t i = 0; i < g.degree(v); ++i) {
      if (ops) SSSP_COUNT(++ops->edge_scans);
      const auto [x, id] = g.arc(v, i);
      if (discovered[x]) continue;
      discovered[x] = 1;
      t.parent[x] = v;
      t.parent_edge[x] = id;
      t.order.push_back(x);
    }
  }
  return t;
}

}  // namespace

Digraph Digraph::from_arcs(std::size_t k, std::span<const DiArc> arcs) {
  if (k >= static_cast<std::size_t>(kNoVertex)) throw std::length_error("too many vertices");
  if (arcs.size() >= static_cast<std::size_t>(kNoEdge)) throw std::length_error("too many arcs");
  Digraph g;
  g.offsets_.assign(k + 1, 0);
  g.indegree_.assign(k, 0);
  for (std::size_t i = 0; i < arcs.size(); ++i) {
    if (arcs[i].tail >= k || arcs[i].head >= k) {
      throw std::out_of_range("arc " + std::to_string(i) + " has an endpoint out of range");
    }
    ++g.offsets_[arcs[i].tail + 1];
    ++g.indegree_[arcs[i].head];
  }
  // Lemma 3.1 counting sort on the tail: prefix sums give start positions, records are written
  // in input order (stable).
  for (std::size_t v = 0; v < k; ++v) g.offsets_[v + 1] += g.offsets_[v];
  g.heads_.resize(arcs.size());
  g.ids_.resize(arcs.size());
  std::vector<std::size_t> cursor(g.offsets_.begin(), g.offsets_.end() - 1);
  for (std::size_t i = 0; i < arcs.size(); ++i) {
    const std::size_t pos = cursor[arcs[i].tail]++;
    g.heads_[pos] = arcs[i].head;
    g.ids_[pos] = static_cast<EdgeId>(i);
  }
  return g;
}

Traversal dfs_order(const Digraph& g, Vertex root, OpCounters* ops) {
  return dfs_impl(DigraphView{g}, root, ops);
}

Traversal bfs_order(const Digraph& g, Vertex root, OpCounters* ops) {
  return bfs_impl(DigraphView{g}, root, ops);
}

template <class W>
Traversal dfs_order(const CsrGraph<W>& g, Vertex root, OpCounters* ops) {
  return dfs_impl(CsrView<W>{g}, root, ops);
}

template <class W>
Traversal bfs_order(const CsrGraph<W>& g, Vertex root, OpCounters* ops) {
  return bfs_impl(CsrView<W>{g}, root, ops);
}

template <class W>
Components connected_components(const CsrGraph<W>& g, OpCounters* ops) {
  const std::size_t k = g.num_vertices();
  Components c;
  c.label.assign(k, kNoVertex);
  std::vector<Vertex> queue;
  queue.reserve(k);
  for (Vertex s = 0; s < k; ++s) {
    if (c.label[s] != kNoVertex) continue;
    const auto id = static_cast<Vertex>(c.count++);
    queue.clear();
    queue.push_back(s);
    c.label[s] = id;
    for (std::size_t head = 0; head < queue.size(); ++head) {
      const Vertex v = queue[head];
      if (ops) SSSP_COUNT(++ops->vertex_scans);
      for (const auto& a : g.arcs(v)) {
        if (ops) SSSP_COUNT(++ops->edge_scans);
        if (c.label[a.to] != kNoVertex) continue;
        c.label[a.to] = id;
        queue.push_back(a.to);
      }
    }
  }
  return c;
}

std::vector<Vertex> topological_sort(const Digraph& g, OpCounters* ops) {
  const std::size_t k = g.num_vertices();
  std::vector<std::size_t> indeg(k);
  std::vector<Vertex> order;
  order.reserve(k);
  for (Vertex v = 0; v < k; ++v) {
    indeg[v] = g.indegree(v);
    if (indeg[v] == 0) order.push_back(v);
  }
  // Lemma 3.1 proof: each removed vertex is queued once and each adjacency entry examined once.
  for (std::size_t head = 0; head < order.size(); ++head) {
    const Vertex v = order[head];
    if (ops) SSSP_COUNT(++ops->vertex_scans);
    for (const Vertex x : g.heads(v)) {
      if (ops) SSSP_COUNT(++ops->edge_scans);
      if (--indeg[x] == 0) order.push_back(x);
    }
  }
  // A nonempty finite acyclic digraph always has a vertex of indegree zero (Lemma 3.1 proof), so
  // a stall before all k vertices are removed certifies a cycle.
  if (order.size() != k) throw std::invalid_argument("topological_sort: the digraph has a cycle");
  return order;
}

template Traversal dfs_order<std::uint64_t>(const CsrGraph<std::uint64_t>&, Vertex, OpCounters*);
template Traversal dfs_order<double>(const CsrGraph<double>&, Vertex, OpCounters*);
template Traversal bfs_order<std::uint64_t>(const CsrGraph<std::uint64_t>&, Vertex, OpCounters*);
template Traversal bfs_order<double>(const CsrGraph<double>&, Vertex, OpCounters*);
template Components connected_components<std::uint64_t>(const CsrGraph<std::uint64_t>&,
                                                        OpCounters*);
template Components connected_components<double>(const CsrGraph<double>&, OpCounters*);

}  // namespace sssp::paper
