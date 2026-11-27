// Sec 3.1 / Lemma 3.2: Dijkstra's search with a radius cutoff (stop when the smallest finite
// value exceeds R, so vertices at distance exactly R are settled), an indexed binary heap, the
// scan-based O(k^2 + l) version, and the search-identifier arena reused across searches.
// Reference: Floyd-Warshall written in this file.
#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

#include "paper/primitives/dijkstra.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::WeightTraits;

namespace {

// Lemma 3.2: the radius search settles exactly the REACHABLE vertices at distance <= R.
template <class W>
bool within(W d, W R) {
  return d != WeightTraits<W>::infinity() && d <= R;
}

template <class W>
std::vector<std::vector<W>> floyd_warshall(const CsrGraph<W>& g) {
  const std::size_t k = g.num_vertices();
  const W inf = WeightTraits<W>::infinity();
  std::vector<std::vector<W>> d(k, std::vector<W>(k, inf));
  for (std::size_t v = 0; v < k; ++v) d[v][v] = W{0};
  for (const auto& e : g.edges()) {
    d[e.u][e.v] = std::min(d[e.u][e.v], e.weight);
    d[e.v][e.u] = std::min(d[e.v][e.u], e.weight);
  }
  for (std::size_t m = 0; m < k; ++m)
    for (std::size_t u = 0; u < k; ++u) {
      if (d[u][m] == inf) continue;
      for (std::size_t v = 0; v < k; ++v) {
        if (d[m][v] == inf) continue;
        d[u][v] = std::min(d[u][v], static_cast<W>(d[u][m] + d[m][v]));
      }
    }
  return d;
}

// Multigraph with zero weights, parallel edges, self-loops and (often) several components.
template <class W>
CsrGraph<W> random_graph(std::mt19937_64& rng, std::size_t k, std::size_t l, std::uint64_t wmax) {
  std::uniform_int_distribution<Vertex> pick(0, static_cast<Vertex>(k - 1));
  std::uniform_int_distribution<std::uint64_t> wt(0, wmax);
  std::vector<Edge<W>> e;
  for (std::size_t i = 0; i < l; ++i) e.push_back({pick(rng), pick(rng), static_cast<W>(wt(rng))});
  return CsrGraph<W>::from_edges(k, std::move(e));
}

// Shortest-path-tree verification against exact reference distances.
template <class W>
void check_tree(const CsrGraph<W>& g, Vertex s, const sssp::SsspResult<W>& r,
                const std::vector<W>& expect) {
  const W inf = WeightTraits<W>::infinity();
  REQUIRE(r.dist == expect);
  REQUIRE(r.parent[s] == kNoVertex);
  for (Vertex v = 0; v < g.num_vertices(); ++v) {
    if (v == s) continue;
    if (r.dist[v] == inf) {
      REQUIRE(r.parent[v] == kNoVertex);
      REQUIRE(r.parent_edge[v] == kNoEdge);
      continue;
    }
    const Vertex p = r.parent[v];
    REQUIRE(p != kNoVertex);
    const auto& e = g.edge(r.parent_edge[v]);
    REQUIRE(((e.u == p && e.v == v) || (e.v == p && e.u == v)));
    REQUIRE(r.dist[p] + e.weight == r.dist[v]);  // tight tree arc
  }
}

template <class W>
void check_against_fw(std::uint64_t seed_base) {
  for (std::uint64_t seed = 0; seed < 40; ++seed) {
    std::mt19937_64 rng(seed_base + seed);
    const std::size_t k = 1 + seed % 24;
    const std::size_t l = (seed * 3) % 70;
    const auto g = random_graph<W>(rng, k, l, seed % 4 == 0 ? 1 : 20);
    const auto fw = floyd_warshall(g);
    for (Vertex s = 0; s < k; ++s) {
      check_tree(g, s, sssp::paper::dijkstra(g, s), fw[s]);
      check_tree(g, s, sssp::paper::dijkstra_scan(g, s), fw[s]);
    }
  }
}

}  // namespace

TEST_CASE("indexed binary heap matches a multiset model", "[paper][dijkstra][heap]") {
  for (std::uint64_t seed = 0; seed < 20; ++seed) {
    std::mt19937_64 rng(seed);
    const std::size_t n = 1 + seed * 7;
    sssp::paper::IndexedBinaryHeap<std::uint64_t> heap;
    heap.reset(n);
    std::map<Vertex, std::uint64_t> model;
    for (int step = 0; step < 2000; ++step) {
      const auto v = static_cast<Vertex>(rng() % n);
      const auto op = rng() % 3;
      if (op == 0 && !model.contains(v)) {
        const std::uint64_t key = rng() % 50;
        heap.push(v, key);
        model[v] = key;
      } else if (op == 1 && model.contains(v) && model[v] > 0) {
        const std::uint64_t key = rng() % model[v];
        heap.decrease(v, key);
        model[v] = key;
      } else if (op == 2 && !model.empty()) {
        const auto [key, u] = heap.pop_min();
        std::uint64_t best = std::numeric_limits<std::uint64_t>::max();
        for (const auto& [x, kx] : model) best = std::min(best, kx);
        REQUIRE(key == best);
        REQUIRE(model.at(u) == key);
        model.erase(u);
      }
      REQUIRE(heap.size() == model.size());
      for (Vertex x = 0; x < n; ++x) REQUIRE(heap.contains(x) == model.contains(x));
    }
  }
}

TEST_CASE("binary-heap and scan Dijkstra equal Floyd-Warshall (uint64)", "[paper][dijkstra]") {
  check_against_fw<std::uint64_t>(0);
}

TEST_CASE("binary-heap and scan Dijkstra equal Floyd-Warshall (double)", "[paper][dijkstra]") {
  // Integer-valued doubles: every sum is exact, so equality with the reference is exact.
  check_against_fw<double>(500);
}

TEST_CASE("Dijkstra boundary inputs", "[paper][dijkstra]") {
  SECTION("empty graph has no valid source") {
    const auto g = CsrGraph<std::uint64_t>::from_edges(0, {});
    REQUIRE_THROWS_AS(sssp::paper::dijkstra(g, 0), std::out_of_range);
    REQUIRE_THROWS_AS(sssp::paper::dijkstra_scan(g, 0), std::out_of_range);
  }
  SECTION("single vertex with a self-loop") {
    const auto g = CsrGraph<std::uint64_t>::from_edges(1, {{0, 0, 5}});
    REQUIRE(sssp::paper::dijkstra(g, 0).dist == std::vector<std::uint64_t>{0});
    REQUIRE(sssp::paper::dijkstra_scan(g, 0).dist == std::vector<std::uint64_t>{0});
  }
  SECTION("zero weights and many equal weights") {
    // A 6x6 grid of unit edges plus a zero-weight diagonal 0-7-14-21.
    std::vector<Edge<std::uint64_t>> e;
    for (Vertex r = 0; r < 6; ++r)
      for (Vertex c = 0; c < 6; ++c) {
        if (c + 1 < 6) e.push_back({r * 6 + c, r * 6 + c + 1, 1});
        if (r + 1 < 6) e.push_back({r * 6 + c, (r + 1) * 6 + c, 1});
      }
    e.push_back({0, 7, 0});
    e.push_back({7, 14, 0});
    e.push_back({14, 21, 0});
    const auto g = CsrGraph<std::uint64_t>::from_edges(36, e);
    const auto fw = floyd_warshall(g);
    check_tree(g, 0, sssp::paper::dijkstra(g, 0), fw[0]);
    check_tree(g, 0, sssp::paper::dijkstra_scan(g, 0), fw[0]);
    REQUIRE(sssp::paper::dijkstra(g, 0).dist[21] == 0);
  }
  SECTION("unreachable vertices stay at infinity") {
    const auto g = CsrGraph<double>::from_edges(4, {{0, 1, 2.5}, {2, 3, 1.0}});
    const auto r = sssp::paper::dijkstra(g, 0);
    REQUIRE(r.dist[1] == 2.5);
    REQUIRE(std::isinf(r.dist[2]));
    REQUIRE(r.parent[3] == kNoVertex);
  }
  SECTION("a NaN radius is rejected") {
    const auto g = CsrGraph<double>::from_edges(2, {{0, 1, 1.0}});
    REQUIRE_THROWS_AS(sssp::paper::dijkstra(g, 0, std::nan("")), std::invalid_argument);
    REQUIRE_THROWS_AS(sssp::paper::dijkstra_scan(g, 0, std::nan("")), std::invalid_argument);
  }
}

TEST_CASE("radius cutoff settles exactly the vertices at distance at most R",
          "[paper][dijkstra][radius]") {
  const auto inf = WeightTraits<std::uint64_t>::infinity();
  for (std::uint64_t seed = 0; seed < 40; ++seed) {
    std::mt19937_64 rng(900 + seed);
    const std::size_t k = 2 + seed % 20;
    const auto g = random_graph<std::uint64_t>(rng, k, 2 * k, seed % 3 == 0 ? 2 : 9);
    const auto fw = floyd_warshall(g);
    for (Vertex s = 0; s < k; ++s) {
      // Every exact distance (boundary R = d), one above and one below it, R = 0 and R = inf.
      std::set<std::uint64_t> radii{0, inf};
      for (const auto d : fw[s]) {
        if (d == inf) continue;
        radii.insert(d);
        radii.insert(d + 1);
        if (d > 0) radii.insert(d - 1);
      }
      for (const auto R : radii) {
        sssp::paper::SearchArena<std::uint64_t> arena(k);
        const std::size_t settled = sssp::paper::dijkstra_search(g, s, R, arena);
        sssp::paper::SearchArena<std::uint64_t> arena2(k);
        const std::size_t settled2 = sssp::paper::dijkstra_scan_search(g, s, R, arena2);
        std::size_t expect = 0;
        for (Vertex v = 0; v < k; ++v) {
          const bool in_ball = within(fw[s][v], R);
          expect += in_ball ? 1 : 0;
          REQUIRE(arena.settled(v) == in_ball);
          REQUIRE(arena2.settled(v) == in_ball);
          if (in_ball) {
            REQUIRE(arena.dist(v) == fw[s][v]);
            REQUIRE(arena2.dist(v) == fw[s][v]);
          }
        }
        REQUIRE(settled == expect);
        REQUIRE(settled2 == expect);
        REQUIRE(arena.settled_order().size() == expect);
        // Settled in nondecreasing distance order.
        for (std::size_t i = 1; i < arena.settled_order().size(); ++i) {
          REQUIRE(arena.dist(arena.settled_order()[i - 1]) <= arena.dist(arena.settled_order()[i]));
        }
        // The convenience form reports unsettled vertices as unreachable.
        const auto res = sssp::paper::dijkstra(g, s, R);
        for (Vertex v = 0; v < k; ++v) {
          REQUIRE(res.dist[v] == (within(fw[s][v], R) ? fw[s][v] : inf));
          if (fw[s][v] > R) REQUIRE(res.parent[v] == kNoVertex);
        }
      }
    }
  }
}

TEST_CASE("radius boundary: distance exactly R is settled, the search stops above R",
          "[paper][dijkstra][radius]") {
  // Path 0 -1- 1 -2- 2 -0- 3 -1- 4: distances 0, 1, 3, 3, 4.
  const auto g =
      CsrGraph<double>::from_edges(5, {{0, 1, 1.0}, {1, 2, 2.0}, {2, 3, 0.0}, {3, 4, 1.0}});
  sssp::paper::SearchArena<double> arena(5);
  REQUIRE(sssp::paper::dijkstra_search(g, 0, 3.0, arena) == 4);
  REQUIRE(arena.settled(2));
  REQUIRE(arena.settled(3));
  REQUIRE_FALSE(arena.settled(4));
  REQUIRE(arena.touched(4));  // offered (tentative 4 > R) but never settled
  REQUIRE(arena.dist(4) == 4.0);
  REQUIRE(sssp::paper::dijkstra_search(g, 0, 0.0, arena) == 1);
  REQUIRE(sssp::paper::dijkstra_search(g, 0, 2.999, arena) == 2);
  // R = 0 also settles zero-distance vertices.
  REQUIRE(sssp::paper::dijkstra_search(g, 2, 0.0, arena) == 2);
  REQUIRE(arena.settled(3));
  // Near the integer limit the cutoff stops before an overflowing offer is formed.
  const std::uint64_t half = std::uint64_t{1} << 63;
  const auto big = CsrGraph<std::uint64_t>::from_edges(3, {{0, 1, half}, {1, 2, half}});
  REQUIRE(sssp::paper::dijkstra(big, 0, half - 1).dist[1] ==
          WeightTraits<std::uint64_t>::infinity());
  REQUIRE_THROWS_AS(sssp::paper::dijkstra(big, 0), std::overflow_error);
}

TEST_CASE("a reused arena initializes only the entries a search touches",
          "[paper][dijkstra][arena]") {
  std::mt19937_64 rng(77);
  const auto g1 = random_graph<std::uint64_t>(rng, 30, 45, 9);
  const auto g2 = random_graph<std::uint64_t>(rng, 60, 50, 9);  // larger: the arena grows
  const auto fw1 = floyd_warshall(g1);
  const auto fw2 = floyd_warshall(g2);
  const auto inf = WeightTraits<std::uint64_t>::infinity();
  sssp::paper::SearchArena<std::uint64_t> arena;
  for (int round = 0; round < 6; ++round) {
    for (Vertex s = 0; s < 60; ++s) {
      const bool first = (s + round) % 2 == 0 && s < 30;
      const auto& g = first ? g1 : g2;
      const auto& fw = first ? fw1 : fw2;
      const std::uint64_t R = (s % 3 == 0) ? inf : 8;
      const std::size_t n = sssp::paper::dijkstra_search(g, s, R, arena);
      std::size_t reached = 0;
      for (Vertex v = 0; v < g.num_vertices(); ++v) {
        const bool in_ball = within(fw[s][v], R);
        reached += in_ball ? 1 : 0;
        REQUIRE(arena.settled(v) == in_ball);
        if (in_ball) REQUIRE(arena.dist(v) == fw[s][v]);
        // Nothing leaks from an earlier search: untouched entries read as infinity.
        if (!arena.touched(v)) {
          REQUIRE(arena.dist(v) == inf);
          REQUIRE(arena.parent(v) == kNoVertex);
        }
      }
      REQUIRE(n == reached);
      // The touched list holds exactly the touched entries, each once.
      std::set<Vertex> touched(arena.touched_vertices().begin(), arena.touched_vertices().end());
      REQUIRE(touched.size() == arena.touched_vertices().size());
      for (Vertex v = 0; v < g.num_vertices(); ++v)
        REQUIRE(touched.contains(v) == arena.touched(v));
    }
  }
  // A search from an isolated vertex touches one entry.
  const auto iso = CsrGraph<std::uint64_t>::from_edges(60, {{0, 1, 1}});
  REQUIRE(sssp::paper::dijkstra_search(iso, 59, inf, arena) == 1);
  REQUIRE(arena.touched_vertices().size() == 1);
}

TEST_CASE("search identifiers survive wrap-around", "[paper][dijkstra][arena]") {
  // An 8-bit identifier wraps after 255 searches; the arena then clears its stamps once.
  std::mt19937_64 rng(4242);
  const auto g = random_graph<std::uint64_t>(rng, 25, 40, 9);
  const auto fw = floyd_warshall(g);
  sssp::paper::SearchArena<std::uint64_t, std::uint8_t> arena(25);
  for (int i = 0; i < 800; ++i) {
    const auto s = static_cast<Vertex>((i * 7) % 25);
    const std::uint64_t R = (i % 2 == 0) ? WeightTraits<std::uint64_t>::infinity() : 5;
    if (i % 3 == 0) {
      (void)sssp::paper::dijkstra_scan_search(g, s, R, arena);
    } else {
      (void)sssp::paper::dijkstra_search(g, s, R, arena);
    }
    for (Vertex v = 0; v < 25; ++v) {
      REQUIRE(arena.settled(v) == within(fw[s][v], R));
      if (within(fw[s][v], R)) REQUIRE(arena.dist(v) == fw[s][v]);
    }
  }
  REQUIRE(arena.full_resets() >= 3);
}

#if SSSP_DIAGNOSTICS  // counters are compiled out when SSSP_DIAGNOSTICS=0
TEST_CASE("operation counts respect Lemma 3.2", "[paper][dijkstra][lemma3.2]") {
  std::mt19937_64 rng(31337);
  // Connected: a spanning path plus random edges.
  const std::size_t k = 400;
  std::vector<Edge<std::uint64_t>> e;
  for (Vertex v = 1; v < k; ++v) e.push_back({v - 1, v, rng() % 50});
  for (int i = 0; i < 1600; ++i) {
    e.push_back({static_cast<Vertex>(rng() % k), static_cast<Vertex>(rng() % k), rng() % 50});
  }
  const auto g = CsrGraph<std::uint64_t>::from_edges(k, std::move(e));
  const std::uint64_t l = g.num_arcs();
  SECTION("scan version: O(k^2 + l)") {
    const auto r = sssp::paper::dijkstra_scan(g, 0);
    const auto& ops = r.diag.ops;
    REQUIRE(ops.vertex_scans == k);
    REQUIRE(ops.edge_scans == l);
    REQUIRE(ops.additions == l);
    // At most k - 1 comparisons per minimum scan, one radius test per extraction, one offer
    // comparison per arc.
    REQUIRE(ops.comparisons <= k * k + l);
  }
  SECTION("binary heap: O(1 + log(k + 1)) per heap operation") {
    const auto r = sssp::paper::dijkstra(g, 0);
    const auto& ops = r.diag.ops;
    REQUIRE(ops.vertex_scans == k);
    REQUIRE(ops.edge_scans == l);
    REQUIRE(ops.queue_pops == k);
    REQUIRE(ops.queue_pushes == k);
    const std::uint64_t heap_ops = ops.queue_pushes + ops.decrease_keys + ops.queue_pops;
    const std::uint64_t height = static_cast<std::uint64_t>(std::bit_width(k));
    // One root-leaf path per operation: <= 2 comparisons per level; plus one offer comparison
    // per arc and one radius test per extraction.
    REQUIRE(ops.comparisons <= heap_ops * 2 * (height + 1) + l + ops.queue_pops);
  }
}
#endif
