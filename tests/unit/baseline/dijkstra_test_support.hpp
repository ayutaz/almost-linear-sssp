#pragma once
// Shared checks for every Dijkstra baseline. Expected values of the hand-made graphs were
// computed on paper; random graphs are compared with a naive Floyd–Warshall written here (test
// only, independent of the code under test). Shortest-path trees are verified structurally,
// because several different trees can be correct.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::test {

template <class W>
using Matrix = std::vector<std::vector<W>>;

/// All-pairs distances by Floyd–Warshall, O(n^3). Self-loops never shorten a path.
template <class W>
Matrix<W> floyd_warshall(const CsrGraph<W>& g) {
  const auto inf = WeightTraits<W>::infinity();
  const std::size_t n = g.num_vertices();
  Matrix<W> d(n, std::vector<W>(n, inf));
  for (std::size_t i = 0; i < n; ++i) d[i][i] = W{0};
  for (const auto& e : g.edges()) {
    if (e.weight < d[e.u][e.v]) d[e.u][e.v] = d[e.v][e.u] = e.weight;
  }
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      if (d[i][k] == inf) continue;
      for (std::size_t j = 0; j < n; ++j) {
        if (d[k][j] == inf) continue;
        const W through = d[i][k] + d[k][j];
        if (through < d[i][j]) d[i][j] = through;
      }
    }
  }
  return d;
}

/// Structural check of a result: sizes, source, unreachable vertices, and for every reachable
/// v != s that parent_edge exists, joins parent[v] and v, and dist[v] == dist[parent] + w
/// exactly; following parents from any vertex reaches the source (no zero-weight cycles).
template <class W>
void check_shortest_path_tree(const CsrGraph<W>& g, Vertex s, const SsspResult<W>& r) {
  const auto inf = WeightTraits<W>::infinity();
  const std::size_t n = g.num_vertices();
  REQUIRE(r.dist.size() == n);
  REQUIRE(r.parent.size() == n);
  REQUIRE(r.parent_edge.size() == n);
  REQUIRE(r.dist[s] == W{0});
  REQUIRE(r.parent[s] == kNoVertex);
  REQUIRE(r.parent_edge[s] == kNoEdge);
  for (Vertex v = 0; v < n; ++v) {
    if (v == s) continue;
    CAPTURE(s, v);
    if (r.dist[v] == inf) {
      REQUIRE(r.parent[v] == kNoVertex);
      REQUIRE(r.parent_edge[v] == kNoEdge);
      continue;
    }
    const Vertex p = r.parent[v];
    const EdgeId id = r.parent_edge[v];
    REQUIRE(p < n);
    REQUIRE(id < g.num_edges());
    const auto& e = g.edge(id);
    REQUIRE(((e.u == p && e.v == v) || (e.u == v && e.v == p)));
    REQUIRE(r.dist[p] != inf);
    REQUIRE(r.dist[v] == r.dist[p] + e.weight);
    Vertex x = v;
    std::size_t steps = 0;
    while (x != s && steps <= n) {
      x = r.parent[x];
      ++steps;
    }
    REQUIRE(x == s);
  }
}

template <class W>
void check_against_floyd_warshall(const CsrGraph<W>& g, const Matrix<W>& fw, Vertex s,
                                  const SsspResult<W>& r, bool exact) {
  check_shortest_path_tree(g, s, r);
  for (Vertex v = 0; v < g.num_vertices(); ++v) {
    CAPTURE(s, v);
    if (exact || fw[s][v] == WeightTraits<W>::infinity()) {
      REQUIRE(r.dist[v] == fw[s][v]);
    } else {
      REQUIRE_THAT(static_cast<double>(r.dist[v]),
                   Catch::Matchers::WithinRel(static_cast<double>(fw[s][v]), 1e-12) ||
                       Catch::Matchers::WithinAbs(static_cast<double>(fw[s][v]), 1e-12));
    }
  }
}

// ------------------------------------------------------------------ hand-made graphs

/// Classic 6-vertex example. Paper calculation from source 0:
///   dist = {0, 7, 9, 20, 20, 11}; 7 successful relaxations
///   (0->1, 0->2, 0->5, 1->3 @22, 2->3 @20, 2->5 @11, 5->4 @20).
template <class W>
CsrGraph<W> classic_graph() {
  return CsrGraph<W>::from_edges(6, {{0, 1, W{7}},
                                     {0, 2, W{9}},
                                     {0, 5, W{14}},
                                     {1, 2, W{10}},
                                     {1, 3, W{15}},
                                     {2, 3, W{11}},
                                     {2, 5, W{2}},
                                     {3, 4, W{6}},
                                     {4, 5, W{9}}});
}

template <class W, class Algo>
void expect_classic_graph(Algo algo) {
  const auto g = classic_graph<W>();
  const auto r = algo(g, Vertex{0});
  const std::vector<W> expected{W{0}, W{7}, W{9}, W{20}, W{20}, W{11}};
  REQUIRE(r.dist == expected);
  check_shortest_path_tree(g, 0, r);
  // Every shortest path here is unique, so the tree is too.
  REQUIRE(r.parent == std::vector<Vertex>{kNoVertex, 0, 0, 2, 5, 2});
  REQUIRE(r.parent_edge == std::vector<EdgeId>{kNoEdge, 0, 1, 5, 8, 6});
}

/// Zero weights, a parallel edge, two self-loops (one zero-weight), a second component {4,5}
/// and an isolated vertex 6. Paper calculation:
///   from 0: {0, 0, 0, 4, inf, inf, inf}, tree 0 -e0- 1 -e1- 2 -e5- 3 (unique)
///   from 4: {inf, inf, inf, inf, 0, 1, inf}
///   from 6: only dist[6] = 0
template <class W>
CsrGraph<W> special_graph() {
  return CsrGraph<W>::from_edges(7, {{0, 1, W{0}},
                                     {1, 2, W{0}},
                                     {0, 2, W{5}},
                                     {0, 1, W{3}},
                                     {2, 2, W{4}},
                                     {2, 3, W{4}},
                                     {3, 3, W{0}},
                                     {1, 3, W{9}},
                                     {4, 5, W{1}}});
}

template <class W, class Algo>
void expect_special_graph(Algo algo) {
  const auto inf = WeightTraits<W>::infinity();
  const auto g = special_graph<W>();

  const auto r0 = algo(g, Vertex{0});
  REQUIRE(r0.dist == std::vector<W>{W{0}, W{0}, W{0}, W{4}, inf, inf, inf});
  check_shortest_path_tree(g, 0, r0);
  REQUIRE(r0.parent == std::vector<Vertex>{kNoVertex, 0, 1, 2, kNoVertex, kNoVertex, kNoVertex});
  REQUIRE(r0.parent_edge == std::vector<EdgeId>{kNoEdge, 0, 1, 5, kNoEdge, kNoEdge, kNoEdge});

  const auto r4 = algo(g, Vertex{4});
  REQUIRE(r4.dist == std::vector<W>{inf, inf, inf, inf, W{0}, W{1}, inf});
  check_shortest_path_tree(g, 4, r4);
  REQUIRE(r4.parent[5] == 4);
  REQUIRE(r4.parent_edge[5] == 8);

  const auto r6 = algo(g, Vertex{6});
  REQUIRE(r6.dist == std::vector<W>{inf, inf, inf, inf, inf, inf, W{0}});
  check_shortest_path_tree(g, 6, r6);
}

/// A zero-weight cycle 0-1-2-0 hanging off the source must not create a parent cycle.
template <class W, class Algo>
void expect_zero_weight_cycle(Algo algo) {
  const auto g = CsrGraph<W>::from_edges(
      4, {{0, 1, W{0}}, {1, 2, W{0}}, {2, 0, W{0}}, {2, 3, W{0}}, {3, 1, W{0}}});
  for (Vertex s = 0; s < 4; ++s) {
    const auto r = algo(g, s);
    REQUIRE(r.dist == std::vector<W>(4, W{0}));
    check_shortest_path_tree(g, s, r);
  }
}

template <class W, class Algo>
void expect_single_vertex(Algo algo) {
  const auto lone = CsrGraph<W>::from_edges(1, {});
  const auto r = algo(lone, Vertex{0});
  REQUIRE(r.dist == std::vector<W>{W{0}});
  REQUIRE(r.parent == std::vector<Vertex>{kNoVertex});
  REQUIRE(r.parent_edge == std::vector<EdgeId>{kNoEdge});

  const auto looped = CsrGraph<W>::from_edges(1, {{0, 0, W{3}}});
  const auto r2 = algo(looped, Vertex{0});
  REQUIRE(r2.dist == std::vector<W>{W{0}});
  REQUIRE(r2.parent == std::vector<Vertex>{kNoVertex});
}

/// Dyadic weights are exact in binary floating point: {0, 0.5, 0.75, 0.875}.
template <class Algo>
void expect_dyadic_doubles(Algo algo) {
  const auto g = CsrGraph<double>::from_edges(
      4, {{0, 1, 0.5}, {1, 2, 0.25}, {0, 2, 1.0}, {2, 3, 0.125}, {3, 0, 2.0}});
  const auto r = algo(g, Vertex{0});
  REQUIRE(r.dist == std::vector<double>{0.0, 0.5, 0.75, 0.875});
  REQUIRE(r.parent == std::vector<Vertex>{kNoVertex, 0, 1, 2});
  check_shortest_path_tree(g, 0, r);
}

template <class W, class Algo>
void expect_invalid_source_rejected(Algo algo) {
  const auto g = CsrGraph<W>::from_edges(3, {{0, 1, W{1}}});
  REQUIRE_THROWS_AS(algo(g, Vertex{3}), std::out_of_range);
  REQUIRE_THROWS_AS(algo(g, kNoVertex), std::out_of_range);
  const auto empty = CsrGraph<W>::from_edges(0, {});
  REQUIRE_THROWS_AS(algo(empty, Vertex{0}), std::out_of_range);
}

template <class Algo>
void expect_integer_overflow_detected(Algo algo) {
  using W = std::uint64_t;
  constexpr W kMax = std::numeric_limits<W>::max();
  // dist[1] = 2^64 - 2 is representable; dist[1] + 5 is not.
  const auto g = CsrGraph<W>::from_edges(3, {{0, 1, kMax - 1}, {1, 2, 5}});
  REQUIRE_THROWS_AS(algo(g, Vertex{0}), std::overflow_error);
  // 2^62 + 2^62 = 2^63 is fine.
  const auto ok = CsrGraph<W>::from_edges(3, {{0, 1, W{1} << 62}, {1, 2, W{1} << 62}});
  REQUIRE(algo(ok, Vertex{0}).dist == std::vector<W>{0, W{1} << 62, W{1} << 63});
}

#if SSSP_DIAGNOSTICS
/// Counters on the classic graph from source 0 (paper calculation, see classic_graph()):
/// 6 vertices scanned, 18 arcs inspected (one addition and one comparison each), 7 relaxations.
/// lazy heap: 1 + 7 pushes, all 8 popped (2 stale); indexed heap: 6 pushes, 2 decrease-keys
/// (3: 22 -> 20 and 5: 14 -> 11), 6 pops.
template <class W, class Algo>
void expect_classic_counters(Algo algo, bool decrease_key_heap) {
  const auto r = algo(classic_graph<W>(), Vertex{0});
  const auto& ops = r.diag.ops;
  REQUIRE(ops.vertex_scans == 6);
  REQUIRE(ops.edge_scans == 18);
  REQUIRE(ops.additions == 18);
  REQUIRE(ops.relaxations == 7);
  if (decrease_key_heap) {
    REQUIRE(ops.queue_pushes == 6);
    REQUIRE(ops.decrease_keys == 2);
    REQUIRE(ops.queue_pops == 6);
    REQUIRE(ops.comparisons >= 18);
  } else {
    REQUIRE(ops.queue_pushes == 8);
    REQUIRE(ops.decrease_keys == 0);
    REQUIRE(ops.queue_pops == 8);
    REQUIRE(ops.comparisons >= 18 + 8);  // + one staleness check per pop
  }
}
#endif

// ------------------------------------------------------------------ random graphs

inline constexpr std::uint64_t kGraphSeeds[] = {1, 2, 3, 42, 2026};

/// Random multigraph with self-loops; integer weights in [0, max_w] with ~20% zeros.
template <class W, class WeightGen>
CsrGraph<W> random_graph(std::mt19937_64& rng, std::size_t n, std::size_t m, WeightGen gen) {
  std::uniform_int_distribution<Vertex> pick(0, static_cast<Vertex>(n - 1));
  std::bernoulli_distribution zero(0.2);
  std::vector<Edge<W>> edges;
  edges.reserve(m);
  for (std::size_t i = 0; i < m; ++i) {
    const Vertex u = pick(rng);
    const Vertex v = pick(rng);
    edges.push_back({u, v, zero(rng) ? W{0} : gen(rng)});
  }
  return CsrGraph<W>::from_edges(n, std::move(edges));
}

/// Runs algo from every source of many small random graphs and from source 0 of one larger
/// graph per seed, comparing with Floyd–Warshall. Integer and dyadic weights compare exactly;
/// general doubles with a 1e-12 relative tolerance (different summation order).
template <class W, class Algo, class WeightGen>
void expect_matches_floyd_warshall(Algo algo, WeightGen gen, bool exact) {
  for (const auto seed : kGraphSeeds) {
    std::mt19937_64 rng(seed);
    for (int trial = 0; trial < 25; ++trial) {
      const std::size_t n = std::uniform_int_distribution<std::size_t>(1, 30)(rng);
      const std::size_t m = std::uniform_int_distribution<std::size_t>(0, 3 * n)(rng);
      const auto g = random_graph<W>(rng, n, m, gen);
      const auto fw = floyd_warshall(g);
      for (Vertex s = 0; s < n; ++s) {
        CAPTURE(seed, trial, n, m);
        check_against_floyd_warshall(g, fw, s, algo(g, s), exact);
      }
    }
    const auto g = random_graph<W>(rng, 200, 800, gen);
    const auto fw = floyd_warshall(g);
    CAPTURE(seed);
    check_against_floyd_warshall(g, fw, 0, algo(g, Vertex{0}), exact);
  }
}

inline auto int_weights(std::uint64_t max_w) {
  return [max_w](std::mt19937_64& r) {
    return std::uniform_int_distribution<std::uint64_t>(0, max_w)(r);
  };
}
/// Multiples of 1/16 in [0, 64]: sums of a few hundred stay exact in a double.
inline auto dyadic_weights() {
  return [](std::mt19937_64& r) {
    return static_cast<double>(std::uniform_int_distribution<int>(0, 1024)(r)) / 16.0;
  };
}
inline auto real_weights() {
  return [](std::mt19937_64& r) { return std::uniform_real_distribution<double>(0.0, 10.0)(r); };
}

}  // namespace sssp::test
