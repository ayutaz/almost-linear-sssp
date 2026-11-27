#pragma once
// Test helpers for the Sec 3.2 MST primitives: seeded random multigraphs and naive checks that
// do not use the code under test.

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

#include "sssp/graph.hpp"

namespace mst_test {

/// Weight r for integer types; r * 0.1 for double (not dyadic: only comparisons are made).
template <class W>
W make_weight(std::uint64_t r) {
  if constexpr (std::is_floating_point_v<W>) {
    return static_cast<W>(r) * W{0.1};
  } else {
    return static_cast<W>(r);
  }
}

/// Random multigraph on n vertices: if `connected`, a random spanning tree first, then
/// `extra` random edges (parallel edges and loops allowed). Weights in [0, max_weight];
/// a small max_weight produces many ties. The edge list is shuffled so ids are not ordered.
template <class W>
std::vector<sssp::Edge<W>> random_edges(std::size_t n, std::size_t extra, std::uint64_t seed,
                                        std::uint64_t max_weight, bool connected) {
  std::mt19937_64 rng(seed);
  std::vector<sssp::Edge<W>> edges;
  auto weight = [&] { return make_weight<W>(rng() % (max_weight + 1)); };
  if (connected) {
    for (std::size_t v = 1; v < n; ++v) {
      const auto u = static_cast<sssp::Vertex>(rng() % v);
      edges.push_back({u, static_cast<sssp::Vertex>(v), weight()});
    }
  }
  if (n > 0) {
    for (std::size_t i = 0; i < extra; ++i) {
      const auto u = static_cast<sssp::Vertex>(rng() % n);
      const auto v = static_cast<sssp::Vertex>(rng() % n);
      edges.push_back({u, v, weight()});
    }
  }
  std::shuffle(edges.begin(), edges.end(), rng);
  return edges;
}

/// Naive union-find with full relabelling (O(n) per union), independent of the library.
class NaiveLabels {
 public:
  explicit NaiveLabels(std::size_t n) : label_(n) { std::iota(label_.begin(), label_.end(), 0); }
  [[nodiscard]] std::size_t label(std::size_t v) const { return label_[v]; }
  bool unite(std::size_t a, std::size_t b) {
    const std::size_t la = label_[a];
    const std::size_t lb = label_[b];
    if (la == lb) return false;
    for (auto& l : label_) {
      if (l == lb) l = la;
    }
    return true;
  }

 private:
  std::vector<std::size_t> label_;
};

template <class W>
std::size_t count_components(std::size_t n, const std::vector<sssp::Edge<W>>& edges) {
  NaiveLabels d(n);
  std::size_t c = n;
  for (const auto& e : edges) {
    if (d.unite(e.u, e.v)) --c;
  }
  return c;
}

/// Matrix of optional weights: `has[i][j]` false means "no walk".
template <class W>
struct PairTable {
  std::vector<std::vector<W>> value;
  std::vector<std::vector<bool>> has;
  explicit PairTable(std::size_t n)
      : value(n, std::vector<W>(n, W{0})), has(n, std::vector<bool>(n, false)) {}
};

/// Naive minimax walk values for all pairs (Floyd–Warshall over (min, max)); O(n^3).
/// value[x][y] = least possible maximum edge length of an x–y walk. Uses comparisons only.
template <class W>
PairTable<W> minimax_all_pairs(std::size_t n, const std::vector<sssp::Edge<W>>& edges) {
  PairTable<W> t(n);
  for (const auto& e : edges) {
    if (e.u == e.v) continue;
    for (const auto& [a, b] : {std::pair{e.u, e.v}, std::pair{e.v, e.u}}) {
      if (!t.has[a][b] || e.weight < t.value[a][b]) {
        t.value[a][b] = e.weight;
        t.has[a][b] = true;
      }
    }
  }
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      if (!t.has[i][k]) continue;
      for (std::size_t j = 0; j < n; ++j) {
        if (!t.has[k][j] || i == j) continue;
        const W via = std::max(t.value[i][k], t.value[k][j]);
        if (!t.has[i][j] || via < t.value[i][j]) {
          t.value[i][j] = via;
          t.has[i][j] = true;
        }
      }
    }
  }
  return t;
}

/// b_T(x, y): heaviest edge on the x–y path of the forest `tree` (edge ids into `edges`),
/// by a search from every vertex; O(n^2).
template <class W>
PairTable<W> tree_bottlenecks(std::size_t n, const std::vector<sssp::Edge<W>>& edges,
                              const std::vector<sssp::EdgeId>& tree) {
  std::vector<std::vector<std::pair<std::size_t, W>>> adj(n);
  for (const auto id : tree) {
    adj[edges[id].u].emplace_back(edges[id].v, edges[id].weight);
    adj[edges[id].v].emplace_back(edges[id].u, edges[id].weight);
  }
  PairTable<W> t(n);
  for (std::size_t s = 0; s < n; ++s) {
    std::vector<std::size_t> stack{s};
    std::vector<bool> seen(n, false);
    seen[s] = true;
    while (!stack.empty()) {
      const std::size_t x = stack.back();
      stack.pop_back();
      for (const auto& [y, w] : adj[x]) {
        if (seen[y]) continue;
        seen[y] = true;
        t.value[s][y] = (x == s) ? w : std::max(t.value[s][x], w);
        t.has[s][y] = true;
        stack.push_back(y);
      }
    }
  }
  return t;
}

}  // namespace mst_test
