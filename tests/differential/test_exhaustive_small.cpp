// Exhaustive differential test: Bellman–Ford vs an in-test Floyd–Warshall on EVERY weighted
// simple graph with few vertices, from every source, and verify_sssp on every result.

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <string>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/verify.hpp"
#include "support/enumerate.hpp"
#include "support/floyd_warshall.hpp"
#include "support/graphs.hpp"

namespace {

template <class W>
std::string describe(std::size_t n, const std::vector<sssp::Edge<W>>& edges, sssp::Vertex s) {
  std::string out = "n=" + std::to_string(n) + " source=" + std::to_string(s) + " edges:";
  for (const auto& e : edges) {
    out += " {" + std::to_string(e.u) + "," + std::to_string(e.v) + "," + std::to_string(e.weight) +
           "}";
  }
  return out;
}

/// Returns the number of (graph, source) pairs checked; stops at the first counterexample.
template <class W>
std::uint64_t check_all(std::size_t n, const std::vector<W>& weights) {
  std::uint64_t graphs = 0;
  std::uint64_t checked = 0;
  bool failed = false;
  sssp::test::for_each_small_graph<W>(n, weights, [&](const std::vector<sssp::Edge<W>>& edges) {
    ++graphs;
    if (failed) return;
    const auto g = sssp::CsrGraph<W>::from_edges(n, edges);
    const auto fw = sssp::test::floyd_warshall<W>(n, edges);
    for (sssp::Vertex s = 0; s < n; ++s) {
      const auto r = sssp::baseline::bellman_ford(g, s);
      const auto rep = sssp::verify_sssp(g, s, r, &fw[s]);
      ++checked;
      if (!rep.ok || r.dist != fw[s]) {
        failed = true;
        FAIL_CHECK("counterexample: " << describe(n, edges, s) << "\n" << rep.summary());
        return;
      }
    }
  });
  REQUIRE_FALSE(failed);
  REQUIRE(graphs == sssp::test::count_small_graphs(n, weights.size()));
  return checked;
}

}  // namespace

TEST_CASE("exhaustive: all simple graphs with n <= 4 and weights {0,1,2}",
          "[differential][exhaustive]") {
  const std::vector<std::uint64_t> w{0, 1, 2};
  std::uint64_t total = 0;
  for (std::size_t n = 1; n <= 4; ++n) total += check_all<std::uint64_t>(n, w);
  // sum over n of n * 4^(n(n-1)/2): 1*1 + 2*4 + 3*64 + 4*4096
  REQUIRE(total == 1 + 8 + 192 + 16384);
}

TEST_CASE("exhaustive: all simple graphs with n = 5 and weights {0,1}",
          "[differential][exhaustive]") {
  REQUIRE(check_all<std::uint64_t>(5, {0, 1}) == 5 * 59049);
}

TEST_CASE("exhaustive: double weights {0, 0.5, 1} on n <= 4 (dyadic sums are exact)",
          "[differential][exhaustive][numerics]") {
  const std::vector<double> w{0.0, 0.5, 1.0};
  for (std::size_t n = 1; n <= 4; ++n) check_all<double>(n, w);
}

TEST_CASE("random multigraphs: Bellman-Ford vs Floyd-Warshall", "[differential][random]") {
  for (std::uint64_t seed = 1; seed <= 300; ++seed) {
    sssp::test::SplitMix64 rng(seed);
    const std::size_t n = 1 + rng.below(25);
    const std::size_t m = rng.below(3 * n + 1);
    const std::uint64_t max_w = rng.below(2) == 0 ? 3 : 1000;
    const auto edges = sssp::test::random_edges<std::uint64_t>(n, m, 0, max_w, seed);
    const auto g = sssp::CsrGraph<std::uint64_t>::from_edges(n, edges);
    const auto fw = sssp::test::floyd_warshall<std::uint64_t>(n, edges);
    const auto s = static_cast<sssp::Vertex>(rng.below(n));
    CAPTURE(seed, n, m, s);
    const auto r = sssp::baseline::bellman_ford(g, s);
    const auto rep = sssp::verify_sssp(g, s, r, &fw[s]);
    INFO(rep.summary());
    REQUIRE(rep.ok);
  }
}
