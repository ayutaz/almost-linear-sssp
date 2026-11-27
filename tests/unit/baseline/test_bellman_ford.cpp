// Bellman–Ford is the independent oracle; its tests use hand-computed distances and a
// Floyd–Warshall implementation that lives in the test tree, never another library algorithm.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "support/floyd_warshall.hpp"
#include "support/graphs.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::baseline::bellman_ford;

using U64 = std::uint64_t;

TEMPLATE_TEST_CASE("Bellman-Ford rejects an invalid source", "[baseline][bellman_ford]", U64,
                   double) {
  const auto g = CsrGraph<TestType>::from_edges(3, {{0, 1, TestType{1}}});
  REQUIRE_THROWS_AS(bellman_ford(g, 3), std::out_of_range);
  REQUIRE_THROWS_AS(bellman_ford(CsrGraph<TestType>::from_edges(0, {}), 0), std::out_of_range);
}

TEMPLATE_TEST_CASE("Bellman-Ford on a single vertex", "[baseline][bellman_ford]", U64, double) {
  const auto r = bellman_ford(CsrGraph<TestType>::from_edges(1, {}), 0);
  REQUIRE(r.dist == std::vector<TestType>{TestType{0}});
  REQUIRE(r.parent == std::vector<sssp::Vertex>{kNoVertex});
  REQUIRE(r.parent_edge == std::vector<sssp::EdgeId>{kNoEdge});
}

TEMPLATE_TEST_CASE("Bellman-Ford computes distances and a tree on a small graph",
                   "[baseline][bellman_ford]", U64, double) {
  using W = TestType;
  //    0 --4-- 1 --1-- 2
  //     \             /
  //      ----- 7 -----        3 isolated
  const std::vector<Edge<W>> e{{0, 1, W{4}}, {1, 2, W{1}}, {0, 2, W{7}}};
  const auto r = bellman_ford(CsrGraph<W>::from_edges(4, e), 0);
  const W inf = sssp::WeightTraits<W>::infinity();
  REQUIRE(r.dist == std::vector<W>{W{0}, W{4}, W{5}, inf});
  REQUIRE(r.parent == std::vector<sssp::Vertex>{kNoVertex, 0, 1, kNoVertex});
  REQUIRE(r.parent_edge == std::vector<sssp::EdgeId>{kNoEdge, 0, 1, kNoEdge});
}

TEMPLATE_TEST_CASE("Bellman-Ford uses the lighter of parallel edges and ignores self-loops",
                   "[baseline][bellman_ford]", U64, double) {
  using W = TestType;
  const std::vector<Edge<W>> e{{0, 1, W{5}}, {1, 1, W{0}}, {1, 0, W{2}}, {0, 1, W{9}}};
  const auto r = bellman_ford(CsrGraph<W>::from_edges(2, e), 0);
  REQUIRE(r.dist[1] == W{2});
  REQUIRE(r.parent[1] == 0);
  REQUIRE(r.parent_edge[1] == 2);
}

TEMPLATE_TEST_CASE("Bellman-Ford handles zero-weight cycles without a parent cycle",
                   "[baseline][bellman_ford]", U64, double) {
  using W = TestType;
  const std::vector<Edge<W>> e{{0, 1, W{0}}, {1, 2, W{0}}, {2, 0, W{0}}, {2, 3, W{1}}};
  const auto r = bellman_ford(CsrGraph<W>::from_edges(4, e), 1);
  REQUIRE(r.dist == std::vector<W>{W{0}, W{0}, W{0}, W{1}});
  // Every parent chain must reach the source within n steps.
  for (sssp::Vertex v = 0; v < 4; ++v) {
    sssp::Vertex x = v;
    for (int steps = 0; steps < 4 && x != 1; ++steps) x = r.parent[x];
    REQUIRE(x == 1);
  }
}

TEMPLATE_TEST_CASE("Bellman-Ford matches Floyd-Warshall on the required graph catalogue",
                   "[baseline][bellman_ford]", U64) {
  for (const auto& c : sssp::test::required_graphs<TestType>()) {
    CAPTURE(c.name);
    const auto fw = sssp::test::floyd_warshall<TestType>(c.n, c.edges);
    const auto r = bellman_ford(c.graph(), c.source);
    REQUIRE(r.dist == fw[c.source]);
  }
}

TEST_CASE("Bellman-Ford with double weights matches Floyd-Warshall where sums are exact",
          "[baseline][bellman_ford]") {
  for (const auto& c : sssp::test::required_graphs<double>()) {
    if (c.name == "extreme_ratio") continue;  // rounding-dependent, covered below
    CAPTURE(c.name);
    const auto fw = sssp::test::floyd_warshall<double>(c.n, c.edges);
    REQUIRE(bellman_ford(c.graph(), c.source).dist == fw[c.source]);
  }
}

TEST_CASE("Bellman-Ford with double weights is exact in the IEEE-754 path algebra",
          "[baseline][bellman_ford][numerics]") {
  // 0 --2^60-- 1 --1-- 2: in double, 2^60 + 1 rounds to 2^60 (absorption).
  const double big = 0x1.0p60;
  const auto g = CsrGraph<double>::from_edges(3, {{0, 1, big}, {1, 2, 1.0}});
  const auto r = bellman_ford(g, 0);
  REQUIRE(r.dist[2] == big);  // fl(2^60 + 1) == 2^60, not the real-number 2^60 + 1
  REQUIRE(r.parent[2] == 1);
}

TEST_CASE("Bellman-Ford detects integer distance overflow", "[baseline][bellman_ford]") {
  const U64 max = std::numeric_limits<U64>::max();
  SECTION("a true shortest distance exceeds the representable range") {
    const auto g = CsrGraph<U64>::from_edges(3, {{0, 1, max - 5}, {1, 2, 10}});
    REQUIRE_THROWS_AS(bellman_ford(g, 0), std::overflow_error);
  }
  SECTION("a distance of exactly infinity() - 1 is representable") {
    const auto g = CsrGraph<U64>::from_edges(3, {{0, 1, max - 3}, {1, 2, 2}});
    REQUIRE(bellman_ford(g, 0).dist[2] == max - 1);
  }
  SECTION("an overflowing sum on a non-shortest path is not an error") {
    // dist[1] + w(1,2) = 1 + (max - 1) overflows, but 2 is reached directly with weight 1.
    const auto g = CsrGraph<U64>::from_edges(3, {{0, 1, 1}, {0, 2, 1}, {1, 2, max - 1}});
    const auto r = bellman_ford(g, 0);
    REQUIRE(r.dist == std::vector<U64>{0, 1, 1});
  }
}

TEMPLATE_TEST_CASE("Bellman-Ford fills comparison-addition counters", "[baseline][bellman_ford]",
                   U64, double) {
  using W = TestType;
  const auto r = bellman_ford(CsrGraph<W>::from_edges(3, {{0, 1, W{1}}, {1, 2, W{1}}}), 0);
#if SSSP_DIAGNOSTICS
  REQUIRE(r.diag.ops.relaxations == 2);
  REQUIRE(r.diag.ops.additions > 0);
  REQUIRE(r.diag.ops.comparisons >= r.diag.ops.additions);
  REQUIRE(r.diag.ops.edge_scans >= r.diag.ops.additions);
  REQUIRE(r.diag.counters.at("rounds") >= 2);
#else
  REQUIRE(r.diag.ops.relaxations == 0);
#endif
}
