// Regression 20261127-regular-switch-dead-end (record:
// cases/20261127-regular-switch-dead-end.json).
//
// include/sssp/generators.hpp promises that MultiEdgePolicy::kSwitchToSimple returns a simple,
// exactly d-regular graph whenever n*d is even and d < n. For (n, d, seed) = (3, 2, 20) the
// configuration model pairs every stub with a stub of the same vertex: three self-loops
// {0,0}, {1,1}, {2,2}. Every proposed swap of two self-loops {a,a}, {b,b} yields {a,b}, {a,b},
// which the switching loop rejects (k1 == k2), so it never made progress and threw
// "random_regular: switching did not converge". About 7% of seeds failed for n = 3, d = 2, and
// similar rates for d = n - 1 (n <= 7); found by tests/differential at SSSP_DIFF_SEEDS=300.
//
// FIXED (P2-DIFF-REGULAR-SWITCH): when the swaps run out of their proposal budget, the
// generator now starts over from a fresh configuration drawn from the same random stream
// (src/graph/generators.cpp: random_regular), so the result is still determined by the
// arguments. Arguments that succeeded before the fix still give the same graph.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "sssp/generators.hpp"

namespace {

/// Independent check of the documented postcondition: simple and every degree exactly d.
bool is_simple_d_regular(const sssp::EdgeList<std::uint64_t>& g, std::size_t n, std::size_t d) {
  if (g.num_vertices != n || g.edges.size() * 2 != n * d) return false;
  std::vector<std::size_t> deg(n, 0);
  std::set<std::pair<sssp::Vertex, sssp::Vertex>> seen;
  for (const auto& e : g.edges) {
    if (e.u == e.v || e.u >= n || e.v >= n) return false;
    if (!seen.insert({std::min(e.u, e.v), std::max(e.u, e.v)}).second) return false;
    ++deg[e.u];
    ++deg[e.v];
  }
  for (const auto x : deg) {
    if (x != d) return false;
  }
  return true;
}

sssp::EdgeList<std::uint64_t> regular_switch(std::size_t n, std::size_t d, std::uint64_t seed) {
  sssp::gen::GenOptions o;
  o.seed = seed;
  return sssp::gen::random_regular<std::uint64_t>(n, d, sssp::gen::MultiEdgePolicy::kSwitchToSimple,
                                                  o);
}

}  // namespace

TEST_CASE("20261127-regular-switch-dead-end: switching from three self-loops yields the triangle",
          "[regression][generators]") {
  // The only simple 2-regular graph on 3 vertices is the triangle.
  const auto g = regular_switch(3, 2, 20);
  CHECK(is_simple_d_regular(g, 3, 2));
}

TEST_CASE("20261127-regular-switch-dead-end: every small valid (n, d) succeeds for seeds 1..50",
          "[regression][generators]") {
  for (std::size_t n = 2; n <= 8; ++n) {
    for (std::size_t d = 1; d < n; ++d) {
      if ((n * d) % 2 != 0) continue;
      for (std::uint64_t seed = 1; seed <= 50; ++seed) {
        INFO("n=" << n << " d=" << d << " seed=" << seed);
        sssp::EdgeList<std::uint64_t> g;
        REQUIRE_NOTHROW(g = regular_switch(n, d, seed));
        REQUIRE(is_simple_d_regular(g, n, d));
      }
    }
  }
}
