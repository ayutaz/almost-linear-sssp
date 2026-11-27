#pragma once
// Minimal, self-contained graph builders for tests.
//
// These deliberately do not depend on the library's graph generators so that a bug in a
// generator cannot hide a bug in the algorithm it is used to test.

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "sssp/graph.hpp"
#include "support/prng.hpp"

namespace sssp::test {

template <class W>
struct GraphCase {
  std::string name;
  std::size_t n = 0;
  std::vector<Edge<W>> edges;
  Vertex source = 0;
  std::uint64_t seed = 0;  ///< PRNG seed for random cases (0 = deterministic construction)

  [[nodiscard]] CsrGraph<W> graph() const { return CsrGraph<W>::from_edges(n, edges); }
};

template <class W>
constexpr W weight_of(std::uint64_t x) noexcept {
  return static_cast<W>(x);
}

// ------------------------------------------------------------------ deterministic shapes

/// 0 - 1 - ... - (n-1), edge i has weight weight_of(i).
template <class W, class WeightFn>
std::vector<Edge<W>> path_edges(std::size_t n, WeightFn weight) {
  std::vector<Edge<W>> e;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    e.push_back({static_cast<Vertex>(i), static_cast<Vertex>(i + 1), weight(i)});
  }
  return e;
}

template <class W>
std::vector<Edge<W>> uniform_path_edges(std::size_t n, W w) {
  return path_edges<W>(n, [w](std::size_t) { return w; });
}

/// Centre 0, leaves 1..n-1, leaf i has weight weight_of(i).
template <class W>
std::vector<Edge<W>> star_edges(std::size_t n) {
  std::vector<Edge<W>> e;
  for (std::size_t i = 1; i < n; ++i) e.push_back({0, static_cast<Vertex>(i), weight_of<W>(i)});
  return e;
}

template <class W>
std::vector<Edge<W>> cycle_edges(std::size_t n, W w) {
  auto e = uniform_path_edges<W>(n, w);
  if (n >= 3) e.push_back({static_cast<Vertex>(n - 1), 0, w});
  return e;
}

/// Complete graph K_n with weights uniform in [lo, hi] drawn from `seed`.
template <class W>
std::vector<Edge<W>> complete_edges(std::size_t n, std::uint64_t lo, std::uint64_t hi,
                                    std::uint64_t seed) {
  SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  for (std::size_t u = 0; u < n; ++u) {
    for (std::size_t v = u + 1; v < n; ++v) {
      e.push_back(
          {static_cast<Vertex>(u), static_cast<Vertex>(v), weight_of<W>(rng.between(lo, hi))});
    }
  }
  return e;
}

/// Random recursive tree: vertex i > 0 attaches to a uniform vertex in [0, i).
template <class W>
std::vector<Edge<W>> random_tree_edges(std::size_t n, std::uint64_t lo, std::uint64_t hi,
                                       std::uint64_t seed) {
  SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  for (std::size_t i = 1; i < n; ++i) {
    const auto p = static_cast<Vertex>(rng.below(i));
    e.push_back({p, static_cast<Vertex>(i), weight_of<W>(rng.between(lo, hi))});
  }
  return e;
}

/// rows x cols grid with all weights equal to w.
template <class W>
std::vector<Edge<W>> grid_edges(std::size_t rows, std::size_t cols, W w) {
  std::vector<Edge<W>> e;
  const auto id = [cols](std::size_t r, std::size_t c) {
    return static_cast<Vertex>(r * cols + c);
  };
  for (std::size_t r = 0; r < rows; ++r) {
    for (std::size_t c = 0; c < cols; ++c) {
      if (c + 1 < cols) e.push_back({id(r, c), id(r, c + 1), w});
      if (r + 1 < rows) e.push_back({id(r, c), id(r + 1, c), w});
    }
  }
  return e;
}

/// `layers` layers of `width` vertices, complete bipartite between consecutive layers, plus a
/// source vertex 0 joined to the whole first layer. All weights w, so the number of distinct
/// shortest paths to a vertex in layer i is width^i.
template <class W>
std::vector<Edge<W>> layered_edges(std::size_t layers, std::size_t width, W w) {
  std::vector<Edge<W>> e;
  const auto id = [width](std::size_t l, std::size_t i) {
    return static_cast<Vertex>(1 + l * width + i);
  };
  for (std::size_t i = 0; i < width; ++i) e.push_back({0, id(0, i), w});
  for (std::size_t l = 0; l + 1 < layers; ++l) {
    for (std::size_t i = 0; i < width; ++i) {
      for (std::size_t j = 0; j < width; ++j) e.push_back({id(l, i), id(l + 1, j), w});
    }
  }
  return e;
}

/// m edges with uniform endpoints (self-loops and multi-edges possible) and weights in [lo, hi].
template <class W>
std::vector<Edge<W>> random_edges(std::size_t n, std::size_t m, std::uint64_t lo, std::uint64_t hi,
                                  std::uint64_t seed) {
  SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  e.reserve(m);
  for (std::size_t i = 0; i < m; ++i) {
    const auto u = static_cast<Vertex>(rng.below(n));
    const auto v = static_cast<Vertex>(rng.below(n));
    e.push_back({u, v, weight_of<W>(rng.between(lo, hi))});
  }
  return e;
}

// ------------------------------------------------------------------ required catalogue

/// The mandatory graph families every SSSP implementation must be tested on.
/// Weights are chosen so that every shortest-path distance fits in uint64_t.
template <class W>
std::vector<GraphCase<W>> required_graphs() {
  std::vector<GraphCase<W>> cs;
  const auto add = [&cs](std::string name, std::size_t n, std::vector<Edge<W>> e, Vertex s = 0,
                         std::uint64_t seed = 0) {
    cs.push_back(GraphCase<W>{std::move(name), n, std::move(e), s, seed});
  };
  const W one = weight_of<W>(1);

  add("single_vertex", 1, {});
  add("no_edges", 5, {}, 2);
  add("two_vertices", 2, {{0, 1, weight_of<W>(3)}});
  add("path", 6, path_edges<W>(6, [](std::size_t i) { return weight_of<W>(i + 1); }));
  add("path_middle_source", 6, path_edges<W>(6, [](std::size_t i) { return weight_of<W>(i + 1); }),
      3);
  add("tree", 15, random_tree_edges<W>(15, 1, 9, 11), 0, 11);
  add("star", 8, star_edges<W>(8));
  add("star_leaf_source", 8, star_edges<W>(8), 5);
  add("cycle", 7, cycle_edges<W>(7, one));
  add("complete", 6, complete_edges<W>(6, 1, 20, 12), 0, 12);
  add("disconnected", 6, {{0, 1, weight_of<W>(2)}, {1, 2, weight_of<W>(2)}, {3, 4, one}});
  add("multi_edges", 3,
      {{0, 1, weight_of<W>(5)},
       {0, 1, weight_of<W>(2)},
       {1, 0, weight_of<W>(7)},
       {1, 2, one},
       {2, 1, one},
       {2, 2, weight_of<W>(0)},
       {0, 0, weight_of<W>(4)}});
  add("zero_weights", 5,
      {{0, 1, weight_of<W>(0)},
       {1, 2, weight_of<W>(0)},
       {2, 0, weight_of<W>(0)},
       {2, 3, weight_of<W>(0)},
       {3, 4, one},
       {0, 4, one}});
  add("uniform_weights", 16, grid_edges<W>(4, 4, weight_of<W>(3)));
  {
    SplitMix64 rng(13);
    auto e = random_edges<W>(12, 30, 1, 1, 13);
    for (auto& edge : e) {
      static constexpr std::uint64_t kMagnitudes[] = {1, 1000, 1000000};
      edge.weight = weight_of<W>(kMagnitudes[rng.below(3)]);
    }
    add("mixed_magnitudes", 12, std::move(e), 0, 13);
  }
  {
    // 2^60 + 1 is not representable in double: exposes absorption under rounding.
    const W big = weight_of<W>(std::uint64_t{1} << 60);
    add("extreme_ratio", 5, {{0, 1, big}, {1, 2, one}, {2, 3, one}, {0, 3, big}, {3, 4, one}});
  }
  {
    // 0-1-...-(n-1) with weight 1, plus chords i -- i+2 of weight 3 (never shorter).
    constexpr std::size_t kN = 300;
    auto e = uniform_path_edges<W>(kN, one);
    for (std::size_t i = 0; i + 2 < kN; ++i) {
      e.push_back({static_cast<Vertex>(i), static_cast<Vertex>(i + 2), weight_of<W>(3)});
    }
    add("long_shortest_path", kN, std::move(e));
  }
  add("many_equal_paths", 1 + 6 * 3, layered_edges<W>(6, 3, one));
  add("many_equal_paths_zero", 1 + 4 * 3, layered_edges<W>(4, 3, weight_of<W>(0)));
  for (std::uint64_t seed = 1; seed <= 5; ++seed) {
    add("random_small_" + std::to_string(seed), 10, random_edges<W>(10, 25, 0, 10, seed),
        static_cast<Vertex>(seed % 10), seed);
  }
  return cs;
}

}  // namespace sssp::test
