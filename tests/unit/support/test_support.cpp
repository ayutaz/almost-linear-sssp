// Self-tests of the test-support helpers (tests/support/). A wrong helper would make every
// test that uses it meaningless, so the helpers are checked against hand-derived facts.

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "support/enumerate.hpp"
#include "support/floyd_warshall.hpp"
#include "support/graphs.hpp"
#include "support/prng.hpp"

using sssp::test::SplitMix64;

TEST_CASE("SplitMix64 matches the reference sequence", "[support][prng]") {
  // Reference values of splitmix64.c (Vigna) for seed 0.
  SplitMix64 rng(0);
  REQUIRE(rng.next() == 0xE220A8397B1DCDAFULL);
  REQUIRE(rng.next() == 0x6E789E6AA1B965F4ULL);
  REQUIRE(rng.next() == 0x06C45D188009454FULL);
}

TEST_CASE("SplitMix64 bounded draws stay in range", "[support][prng]") {
  SplitMix64 rng(42);
  for (int i = 0; i < 10000; ++i) {
    REQUIRE(rng.below(7) < 7);
    const auto x = rng.between(3, 5);
    REQUIRE((x >= 3 && x <= 5));
    const double u = rng.unit();
    REQUIRE((u >= 0.0 && u < 1.0));
  }
}

TEST_CASE("small-graph enumeration is complete and duplicate-free", "[support][enumerate]") {
  REQUIRE(sssp::test::count_small_graphs(0, 3) == 1);
  REQUIRE(sssp::test::count_small_graphs(1, 3) == 1);
  REQUIRE(sssp::test::count_small_graphs(3, 2) == 27);
  REQUIRE(sssp::test::count_small_graphs(4, 3) == 4096);

  const std::vector<std::uint64_t> w{0, 1, 2};
  std::set<std::vector<std::tuple<sssp::Vertex, sssp::Vertex, std::uint64_t>>> seen;
  sssp::test::for_each_small_graph<std::uint64_t>(4, w, [&](const auto& edges) {
    std::vector<std::tuple<sssp::Vertex, sssp::Vertex, std::uint64_t>> key;
    for (const auto& e : edges) {
      REQUIRE(e.u < e.v);  // simple graph: no self-loops, each pair once
      key.emplace_back(e.u, e.v, e.weight);
    }
    seen.insert(key);
  });
  REQUIRE(seen.size() == 4096);
  REQUIRE(seen.count({}) == 1);  // the empty graph is included
}

TEST_CASE("required catalogue contains every mandatory family", "[support][graphs]") {
  const auto cs = sssp::test::required_graphs<std::uint64_t>();
  std::set<std::string> names;
  for (const auto& c : cs) {
    names.insert(c.name);
    REQUIRE(c.source < c.n);
    REQUIRE_NOTHROW(c.graph());
  }
  for (const char* required :
       {"single_vertex", "no_edges", "two_vertices", "path", "tree", "star", "cycle", "complete",
        "disconnected", "multi_edges", "zero_weights", "uniform_weights", "mixed_magnitudes",
        "extreme_ratio", "long_shortest_path", "many_equal_paths", "random_small_1"}) {
    CAPTURE(required);
    REQUIRE(names.count(required) == 1);
  }
}

TEST_CASE("catalogue shapes have their defining properties", "[support][graphs]") {
  const auto cs = sssp::test::required_graphs<std::uint64_t>();
  const auto find = [&](const std::string& name) {
    for (const auto& c : cs) {
      if (c.name == name) return c;
    }
    FAIL("missing " << name);
    return cs.front();
  };
  SECTION("tree has n-1 edges and is connected") {
    const auto c = find("tree");
    REQUIRE(c.edges.size() == c.n - 1);
    const auto fw = sssp::test::floyd_warshall<std::uint64_t>(c.n, c.edges);
    for (const auto d : fw[0]) REQUIRE(d != sssp::WeightTraits<std::uint64_t>::infinity());
  }
  SECTION("long shortest path ends n-1 hops away at distance n-1") {
    const auto c = find("long_shortest_path");
    const auto fw = sssp::test::floyd_warshall<std::uint64_t>(c.n, c.edges);
    REQUIRE(fw[0][c.n - 1] == c.n - 1);
  }
  SECTION("complete graph has n(n-1)/2 edges") {
    const auto c = find("complete");
    REQUIRE(c.edges.size() == c.n * (c.n - 1) / 2);
  }
  SECTION("disconnected graph has unreachable vertices") {
    const auto c = find("disconnected");
    const auto fw = sssp::test::floyd_warshall<std::uint64_t>(c.n, c.edges);
    REQUIRE(fw[0][3] == sssp::WeightTraits<std::uint64_t>::infinity());
    REQUIRE(fw[0][5] == sssp::WeightTraits<std::uint64_t>::infinity());
  }
  SECTION("random cases are reproducible from their seed") {
    REQUIRE(sssp::test::random_edges<std::uint64_t>(10, 25, 0, 10, 3).size() == 25);
    const auto a = sssp::test::random_edges<std::uint64_t>(10, 25, 0, 10, 3);
    const auto b = sssp::test::random_edges<std::uint64_t>(10, 25, 0, 10, 3);
    for (std::size_t i = 0; i < a.size(); ++i) {
      REQUIRE(std::tie(a[i].u, a[i].v, a[i].weight) == std::tie(b[i].u, b[i].v, b[i].weight));
    }
  }
}

TEST_CASE("Floyd-Warshall on a hand-computed graph", "[support][floyd_warshall]") {
  // 0 -4- 1 -1- 2, 0 -7- 2, 3 isolated; parallel edge 0 -3- 1 is lighter.
  const std::vector<sssp::Edge<std::uint64_t>> e{{0, 1, 4}, {1, 2, 1}, {0, 2, 7}, {1, 0, 3}};
  const auto d = sssp::test::floyd_warshall<std::uint64_t>(4, e);
  const std::uint64_t inf = sssp::WeightTraits<std::uint64_t>::infinity();
  REQUIRE(d[0] == std::vector<std::uint64_t>{0, 3, 4, inf});
  REQUIRE(d[2] == std::vector<std::uint64_t>{4, 1, 0, inf});
  REQUIRE(d[3] == std::vector<std::uint64_t>{inf, inf, inf, 0});
}
