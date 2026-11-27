#pragma once
// Test-only generators of rooted forests given as parent arrays (kNoVertex marks a root) and
// naive reference queries that walk parent pointers. Shared by the LCA / level-ancestor tests.

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <optional>
#include <random>
#include <vector>

#include "sssp/graph.hpp"

namespace forest_gen {

using sssp::kNoVertex;
using sssp::Vertex;
using Parents = std::vector<Vertex>;

/// Applies a vertex relabelling: new id of old vertex v is perm[v].
inline Parents relabel(const Parents& parent, const std::vector<Vertex>& perm) {
  Parents out(parent.size(), kNoVertex);
  for (std::size_t v = 0; v < parent.size(); ++v) {
    out[perm[v]] = parent[v] == kNoVertex ? kNoVertex : perm[parent[v]];
  }
  return out;
}

inline Parents shuffle_labels(const Parents& parent, std::mt19937_64& rng) {
  std::vector<Vertex> perm(parent.size());
  std::iota(perm.begin(), perm.end(), Vertex{0});
  std::shuffle(perm.begin(), perm.end(), rng);
  return relabel(parent, perm);
}

/// Path 0 - 1 - ... - (n-1) rooted at 0.
inline Parents path(std::size_t n) {
  Parents p(n, kNoVertex);
  for (std::size_t v = 1; v < n; ++v) p[v] = static_cast<Vertex>(v - 1);
  return p;
}

/// Star with centre 0.
inline Parents star(std::size_t n) {
  Parents p(n, kNoVertex);
  for (std::size_t v = 1; v < n; ++v) p[v] = 0;
  return p;
}

/// Caterpillar: a spine path of `spine` vertices, remaining vertices are leaves hung on the
/// spine round-robin.
inline Parents caterpillar(std::size_t n, std::size_t spine) {
  spine = std::max<std::size_t>(1, std::min(spine, n));
  Parents p(n, kNoVertex);
  for (std::size_t v = 1; v < spine; ++v) p[v] = static_cast<Vertex>(v - 1);
  for (std::size_t v = spine; v < n; ++v) p[v] = static_cast<Vertex>((v - spine) % spine);
  return p;
}

/// Complete binary tree in heap order.
inline Parents binary(std::size_t n) {
  Parents p(n, kNoVertex);
  for (std::size_t v = 1; v < n; ++v) p[v] = static_cast<Vertex>((v - 1) / 2);
  return p;
}

/// Random recursive tree: vertex v attaches to a uniform earlier vertex.
inline Parents random_tree(std::size_t n, std::mt19937_64& rng) {
  Parents p(n, kNoVertex);
  for (std::size_t v = 1; v < n; ++v) {
    p[v] = static_cast<Vertex>(std::uniform_int_distribution<std::size_t>(0, v - 1)(rng));
  }
  return p;
}

/// Random deep tree: vertex v attaches to one of the last `window` vertices (long paths).
inline Parents random_deep_tree(std::size_t n, std::size_t window, std::mt19937_64& rng) {
  Parents p(n, kNoVertex);
  for (std::size_t v = 1; v < n; ++v) {
    const std::size_t lo = v > window ? v - window : 0;
    p[v] = static_cast<Vertex>(std::uniform_int_distribution<std::size_t>(lo, v - 1)(rng));
  }
  return p;
}

/// Random forest with `trees` roots (each vertex v >= trees attaches to an earlier vertex).
inline Parents random_forest(std::size_t n, std::size_t trees, std::mt19937_64& rng) {
  trees = std::max<std::size_t>(1, std::min(trees, n));
  Parents p(n, kNoVertex);
  for (std::size_t v = trees; v < n; ++v) {
    p[v] = static_cast<Vertex>(std::uniform_int_distribution<std::size_t>(0, v - 1)(rng));
  }
  return p;
}

/// Naive depth by walking to the root.
inline std::uint32_t naive_depth(const Parents& parent, Vertex v) {
  std::uint32_t d = 0;
  while (parent[v] != kNoVertex) {
    v = parent[v];
    ++d;
  }
  return d;
}

inline Vertex naive_root(const Parents& parent, Vertex v) {
  while (parent[v] != kNoVertex) v = parent[v];
  return v;
}

/// Naive LCA by walking parent pointers; nullopt if u and v lie in different trees.
inline std::optional<Vertex> naive_lca(const Parents& parent, Vertex u, Vertex v) {
  auto du = naive_depth(parent, u);
  auto dv = naive_depth(parent, v);
  while (du > dv) {
    u = parent[u];
    --du;
  }
  while (dv > du) {
    v = parent[v];
    --dv;
  }
  while (u != v) {
    if (parent[u] == kNoVertex) return std::nullopt;
    u = parent[u];
    v = parent[v];
  }
  return u;
}

/// Naive level ancestor: the ancestor of v at depth j (0 <= j <= depth(v)).
inline Vertex naive_level_ancestor(const Parents& parent, Vertex v, std::uint32_t j) {
  auto d = naive_depth(parent, v);
  while (d > j) {
    v = parent[v];
    --d;
  }
  return v;
}

}  // namespace forest_gen
