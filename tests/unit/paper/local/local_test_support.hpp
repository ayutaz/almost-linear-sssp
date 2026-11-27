#pragma once
// Test support for the Sec 5 local constructions: seeded connected graph families and
// brute-force oracles written independently of src/paper (tree paths by parent walking,
// bottlenecks and sums along them, an O(n^2) Dijkstra on undirected graphs).
//
// Weights are integers (std::uint64_t) or quarter-integers (double) below 2^42, so every sum
// formed by the tests below 2^53 is exact in double as well.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "sssp/graph.hpp"
#include "support/prng.hpp"

namespace sssp::test::local {

enum class Family : std::uint8_t {
  kPath,
  kSpider,      ///< a centre with several legs
  kRandomTree,  ///< random recursive tree, no extra edges
  kSparse,      ///< random tree plus about n extra edges (parallel edges and loops allowed)
  kGrid,
  kComplete,
};

enum class Weights : std::uint8_t {
  kTies,         ///< uniform in [1, 3]: many equal weights
  kUniform,      ///< uniform in [1, 1000]
  kExponential,  ///< 2^U[0, 40]: radius ratios far above K^2, near and far edges
};

inline const char* name(Family f) {
  switch (f) {
    case Family::kPath:
      return "path";
    case Family::kSpider:
      return "spider";
    case Family::kRandomTree:
      return "random-tree";
    case Family::kSparse:
      return "sparse";
    case Family::kGrid:
      return "grid";
    case Family::kComplete:
      return "complete";
  }
  return "?";
}

inline const char* name(Weights w) {
  switch (w) {
    case Weights::kTies:
      return "ties";
    case Weights::kUniform:
      return "uniform";
    case Weights::kExponential:
      return "exponential";
  }
  return "?";
}

/// Integer weight drawn for distribution `w`.
inline std::uint64_t draw_weight(Weights w, SplitMix64& rng) {
  switch (w) {
    case Weights::kTies:
      return rng.between(1, 3);
    case Weights::kUniform:
      return rng.between(1, 1000);
    case Weights::kExponential:
      return std::uint64_t{1} << rng.between(0, 40);
  }
  return 1;
}

template <class W>
W to_weight(std::uint64_t x, SplitMix64& rng) {
  if constexpr (std::is_floating_point_v<W>) {
    // Quarter-integers: exact in every sum below 2^53.
    return static_cast<W>(x) + static_cast<W>(rng.between(0, 3)) * W{0.25};
  } else {
    return static_cast<W>(x);
  }
}

/// A connected graph with n >= 1 vertices of family f and weight distribution wd.
template <class W>
std::vector<Edge<W>> make_graph(Family f, Weights wd, std::size_t n, std::uint64_t seed) {
  SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  const auto add = [&](std::size_t u, std::size_t v) {
    e.push_back(
        {static_cast<Vertex>(u), static_cast<Vertex>(v), to_weight<W>(draw_weight(wd, rng), rng)});
  };
  switch (f) {
    case Family::kPath:
      for (std::size_t v = 1; v < n; ++v) add(v - 1, v);
      break;
    case Family::kSpider: {
      const std::size_t legs = 2 + rng.below(4);
      for (std::size_t v = 1; v < n; ++v) add(v <= legs ? 0 : v - legs, v);
      break;
    }
    case Family::kRandomTree:
      for (std::size_t v = 1; v < n; ++v) add(rng.below(v), v);
      break;
    case Family::kSparse:
      for (std::size_t v = 1; v < n; ++v) add(rng.below(v), v);
      for (std::size_t i = 0; i < n; ++i) add(rng.below(n), rng.below(n));
      break;
    case Family::kGrid: {
      std::size_t side = 1;
      while ((side + 1) * (side + 1) <= n) ++side;
      for (std::size_t v = 0; v < n; ++v) {
        if ((v + 1) % side != 0 && v + 1 < n) add(v, v + 1);
        if (v + side < n) add(v, v + side);
      }
      break;
    }
    case Family::kComplete:
      for (std::size_t u = 0; u < n; ++u) {
        for (std::size_t v = u + 1; v < n; ++v) add(u, v);
      }
      break;
  }
  return e;
}

/// A rooted spanning tree given by its edge ids; parent arrays from a BFS at `root`.
template <class W>
struct RootedTree {
  std::vector<Vertex> parent;
  std::vector<EdgeId> parent_edge;
  std::vector<std::uint32_t> depth;
};

template <class W>
RootedTree<W> root_tree(std::size_t n, const std::vector<Edge<W>>& edges,
                        const std::vector<EdgeId>& tree, Vertex root = 0) {
  std::vector<std::vector<std::pair<Vertex, EdgeId>>> adj(n);
  for (const EdgeId id : tree) {
    adj[edges[id].u].emplace_back(edges[id].v, id);
    adj[edges[id].v].emplace_back(edges[id].u, id);
  }
  RootedTree<W> t{std::vector<Vertex>(n, kNoVertex), std::vector<EdgeId>(n, kNoEdge),
                  std::vector<std::uint32_t>(n, 0)};
  std::vector<bool> seen(n, false);
  std::vector<Vertex> queue{root};
  seen[root] = true;
  for (std::size_t i = 0; i < queue.size(); ++i) {
    const Vertex x = queue[i];
    for (const auto& [y, id] : adj[x]) {
      if (seen[y]) continue;
      seen[y] = true;
      t.parent[y] = x;
      t.parent_edge[y] = id;
      t.depth[y] = t.depth[x] + 1;
      queue.push_back(y);
    }
  }
  return t;
}

/// Edge ids on the tree path between u and v (in no particular order).
template <class W>
std::vector<EdgeId> tree_path(const RootedTree<W>& t, Vertex u, Vertex v) {
  std::vector<EdgeId> out;
  while (u != v) {
    if (t.depth[u] >= t.depth[v]) {
      out.push_back(t.parent_edge[u]);
      u = t.parent[u];
    } else {
      out.push_back(t.parent_edge[v]);
      v = t.parent[v];
    }
  }
  return out;
}

template <class W>
W path_max(const std::vector<Edge<W>>& edges, const std::vector<EdgeId>& path) {
  W m{0};
  for (const EdgeId id : path) m = std::max(m, edges[id].weight);
  return m;
}

template <class W>
W path_sum(const std::vector<Edge<W>>& edges, const std::vector<EdgeId>& path) {
  W s{0};
  for (const EdgeId id : path) s += edges[id].weight;
  return s;
}

/// Total weight of a set of edges.
template <class W>
W total_weight(const std::vector<Edge<W>>& edges, const std::vector<EdgeId>& ids) {
  return path_sum(edges, ids);
}

/// O(n^2 + m) Dijkstra on the undirected multigraph `edges` from several weighted sources.
template <class W>
std::vector<W> dijkstra_naive(std::size_t n, const std::vector<Edge<W>>& edges,
                              const std::vector<std::pair<Vertex, W>>& sources) {
  constexpr W kInf = WeightTraits<W>::infinity();
  std::vector<std::vector<std::pair<Vertex, W>>> adj(n);
  for (const auto& e : edges) {
    adj[e.u].emplace_back(e.v, e.weight);
    adj[e.v].emplace_back(e.u, e.weight);
  }
  std::vector<W> d(n, kInf);
  for (const auto& [s, c] : sources) d[s] = std::min(d[s], c);
  std::vector<bool> done(n, false);
  for (;;) {
    std::size_t best = n;
    for (std::size_t v = 0; v < n; ++v) {
      if (!done[v] && d[v] != kInf && (best == n || d[v] < d[best])) best = v;
    }
    if (best == n) break;
    done[best] = true;
    for (const auto& [y, w] : adj[best]) d[y] = std::min(d[y], d[best] + w);
  }
  return d;
}

}  // namespace sssp::test::local
