#pragma once
// Alternative graph representations and conversions to/from the shared CsrGraph.
//
// Every SSSP implementation runs on CsrGraph so that benchmarks compare identical memory
// layouts. These formats exist for input/output, generators, and for experiments that need a
// mutable or pointer-based layout. Conversions are free functions with no hidden caching, so a
// benchmark can time them separately from the algorithm (see sssp::timed below).
//
// Ordering contract: all conversions preserve the edge order of the edge list (EdgeId i in one
// representation is EdgeId i in every other), and AdjacencyList arcs appear in exactly the same
// per-vertex order as CsrGraph::arcs(v) (for each edge in id order: the u-side arc, then the
// v-side arc; a self-loop {u,u} is one arc).

#include <chrono>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "sssp/graph.hpp"

namespace sssp {

/// Plain undirected edge list: the interchange format used by IO and generators.
template <class W>
struct EdgeList {
  std::size_t num_vertices = 0;
  std::vector<Edge<W>> edges;
};

/// Vector-of-vectors adjacency list with the same Arc layout as CsrGraph.
template <class W>
class AdjacencyList {
 public:
  using weight_type = W;

  AdjacencyList() = default;

  /// Builds the adjacency list from an edge list (validated like CsrGraph::from_edges).
  static AdjacencyList from_edges(std::size_t n, std::vector<Edge<W>> edges) {
    validate_edges<W>(n, edges);
    AdjacencyList a;
    a.edges_ = std::move(edges);
    a.adj_.assign(n, {});
    std::vector<std::size_t> deg(n, 0);
    for (const auto& e : a.edges_) {
      ++deg[e.u];
      if (e.u != e.v) ++deg[e.v];
    }
    for (std::size_t v = 0; v < n; ++v) a.adj_[v].reserve(deg[v]);
    for (std::size_t i = 0; i < a.edges_.size(); ++i) {
      const auto& e = a.edges_[i];
      const auto id = static_cast<EdgeId>(i);
      a.adj_[e.u].push_back(Arc<W>{e.v, e.weight, id});
      if (e.u != e.v) a.adj_[e.v].push_back(Arc<W>{e.u, e.weight, id});
      a.num_arcs_ += (e.u != e.v) ? 2 : 1;
    }
    return a;
  }

  [[nodiscard]] std::size_t num_vertices() const noexcept { return adj_.size(); }
  [[nodiscard]] std::size_t num_edges() const noexcept { return edges_.size(); }
  [[nodiscard]] std::size_t num_arcs() const noexcept { return num_arcs_; }
  [[nodiscard]] std::size_t degree(Vertex v) const noexcept { return adj_[v].size(); }
  [[nodiscard]] std::span<const Arc<W>> arcs(Vertex v) const noexcept { return adj_[v]; }
  [[nodiscard]] const Edge<W>& edge(EdgeId id) const noexcept { return edges_[id]; }
  [[nodiscard]] std::span<const Edge<W>> edges() const noexcept { return edges_; }

 private:
  std::vector<std::vector<Arc<W>>> adj_;
  std::vector<Edge<W>> edges_;
  std::size_t num_arcs_ = 0;
};

// ------------------------------------------------------------------ conversions

template <class W>
CsrGraph<W> to_csr(const EdgeList<W>& el) {
  return CsrGraph<W>::from_edges(el.num_vertices, el.edges);
}

/// Move overload: avoids copying the edge vector.
template <class W>
CsrGraph<W> to_csr(EdgeList<W>&& el) {
  return CsrGraph<W>::from_edges(el.num_vertices, std::move(el.edges));
}

template <class W>
CsrGraph<W> to_csr(const AdjacencyList<W>& a) {
  return CsrGraph<W>::from_edges(a.num_vertices(),
                                 std::vector<Edge<W>>(a.edges().begin(), a.edges().end()));
}

template <class W>
EdgeList<W> to_edge_list(const CsrGraph<W>& g) {
  return EdgeList<W>{g.num_vertices(), std::vector<Edge<W>>(g.edges().begin(), g.edges().end())};
}

template <class W>
EdgeList<W> to_edge_list(const AdjacencyList<W>& a) {
  return EdgeList<W>{a.num_vertices(), std::vector<Edge<W>>(a.edges().begin(), a.edges().end())};
}

template <class W>
AdjacencyList<W> to_adjacency_list(const EdgeList<W>& el) {
  return AdjacencyList<W>::from_edges(el.num_vertices, el.edges);
}

template <class W>
AdjacencyList<W> to_adjacency_list(const CsrGraph<W>& g) {
  return AdjacencyList<W>::from_edges(g.num_vertices(),
                                      std::vector<Edge<W>>(g.edges().begin(), g.edges().end()));
}

// ------------------------------------------------------------------ timing helper

/// Value plus the wall-clock time (steady_clock) spent producing it.
template <class T>
struct Timed {
  T value;
  double seconds = 0.0;
};

/// Runs f() and measures it, so conversion cost can be reported separately from algorithm cost:
///   auto [g, t] = sssp::timed([&] { return sssp::to_csr(edge_list); });
template <class F>
auto timed(F&& f) -> Timed<decltype(std::forward<F>(f)())> {
  const auto start = std::chrono::steady_clock::now();
  auto value = std::forward<F>(f)();
  const auto stop = std::chrono::steady_clock::now();
  return {std::move(value), std::chrono::duration<double>(stop - start).count()};
}

}  // namespace sssp
