// Property tests for the Sec 3.2 MST (Lemmas 3.4 and 3.5) on seeded random graph families.
#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include "paper/primitives/mst.hpp"
#include "sssp/graph.hpp"
#include "unit/paper/primitives/mst_test_support.hpp"

using sssp::Edge;
using sssp::Vertex;
using sssp::paper::kUnboundedGrowth;
using sssp::paper::mst_kruskal_ref;
using sssp::paper::mst_three_pass;
using sssp::paper::spanning_forest_three_pass;
using sssp::paper::three_pass_limits;

namespace {

enum class Family { kRandom, kPath, kStar, kGrid, kComplete, kDisconnected };

template <class W>
std::vector<Edge<W>> make_family(Family f, std::size_t n, std::uint64_t seed,
                                 std::uint64_t max_weight) {
  std::mt19937_64 rng(seed);
  auto wt = [&] { return mst_test::make_weight<W>(rng() % (max_weight + 1)); };
  std::vector<Edge<W>> e;
  switch (f) {
    case Family::kRandom:
      return mst_test::random_edges<W>(n, 3 * n, seed, max_weight, true);
    case Family::kDisconnected:
      return mst_test::random_edges<W>(n, n / 2, seed, max_weight, false);
    case Family::kPath:
      for (Vertex v = 1; v < n; ++v) e.push_back({v - 1, v, wt()});
      break;
    case Family::kStar:
      for (Vertex v = 1; v < n; ++v) e.push_back({0, v, wt()});
      for (Vertex v = 2; v < n; v += 3) e.push_back({v - 1, v, wt()});
      break;
    case Family::kGrid: {
      std::size_t side = 1;
      while ((side + 1) * (side + 1) <= n) ++side;
      for (std::size_t i = 0; i < n; ++i) {
        if ((i + 1) % side != 0 && i + 1 < n) {
          e.push_back({static_cast<Vertex>(i), static_cast<Vertex>(i + 1), wt()});
        }
        if (i + side < n) {
          e.push_back({static_cast<Vertex>(i), static_cast<Vertex>(i + side), wt()});
        }
      }
      break;
    }
    case Family::kComplete:
      for (Vertex u = 0; u < n; ++u) {
        for (Vertex v = u + 1; v < n; ++v) e.push_back({u, v, wt()});
      }
      break;
  }
  std::shuffle(e.begin(), e.end(), rng);
  return e;
}

}  // namespace

TEMPLATE_TEST_CASE("three-pass spanning forest equals Kruskal on seeded graph families",
                   "[paper][mst][property]", std::uint64_t, double) {
  using T = TestType;
  const std::vector<Family> families{Family::kRandom, Family::kPath,     Family::kStar,
                                     Family::kGrid,   Family::kComplete, Family::kDisconnected};
  std::size_t checked_passes = 0;
  for (const Family f : families) {
    for (std::uint64_t seed = 0; seed < 25; ++seed) {
      const std::size_t n = f == Family::kComplete ? 1 + seed % 40 : 1 + (seed * 977 + 13) % 2500;
      const std::uint64_t max_weight = seed % 3 == 0 ? 1 : (seed % 3 == 1 ? 20 : 1'000'000);
      const auto edges = make_family<T>(f, n, seed, max_weight);
      INFO("family " << static_cast<int>(f) << " seed " << seed << " n " << n);
      const std::size_t comps = mst_test::count_components(n, edges);
      const auto ref = mst_kruskal_ref<T>(n, edges);
      const auto r = spanning_forest_three_pass<T>(n, edges);
      REQUIRE(r.edges == ref);
      REQUIRE(r.edges.size() == n - comps);
      if (comps == 1) {
        REQUIRE(mst_three_pass<T>(n, edges).edges == ref);
      } else {
        REQUIRE_THROWS_AS(mst_three_pass<T>(n, edges), std::invalid_argument);
      }
      // Lemma 3.4 per pass (with one empty-heap growth per component) and Lemma 3.5's chain.
      const auto lim = three_pass_limits(n);
      const std::vector<std::size_t> limits{lim.b1, lim.b2, kUnboundedGrowth};
      std::size_t vertices = n;
      std::size_t extractions = 0;
      for (std::size_t i = 0; i < r.passes.size(); ++i) {
        const auto& p = r.passes[i];
        REQUIRE(p.limit == limits[i]);
        REQUIRE(p.vertices_in == vertices);
        REQUIRE(p.max_extract_heap_size <= std::min(p.vertices_in, p.limit));
        if (p.limit != kUnboundedGrowth) {
          REQUIRE(p.vertices_out * p.limit <= comps * p.limit + 2 * p.edges_in);
        } else {
          REQUIRE(p.vertices_out == comps);
        }
        vertices = p.vertices_out;
        extractions += p.extractions;
        ++checked_passes;
      }
      REQUIRE(vertices == comps);
      REQUIRE(extractions == n - comps);
    }
  }
  REQUIRE(checked_passes > 0);
}

TEMPLATE_TEST_CASE("three-pass MST on a large sparse graph runs all passes within bounds",
                   "[paper][mst][property]", std::uint64_t, double) {
  using T = TestType;
  for (std::uint64_t seed : {11ULL, 12ULL, 13ULL}) {
    const std::size_t n = 20000;
    const auto edges = mst_test::random_edges<T>(n, 2 * n, seed, 1u << 20, true);
    const auto r = mst_three_pass<T>(n, edges);
    REQUIRE(r.edges == mst_kruskal_ref<T>(n, edges));
    REQUIRE(r.passes.size() == 3);  // B1 = 4 and B2 = 15 stop growths on this input
    const auto lim = three_pass_limits(n);
    REQUIRE(lim.b1 == 4);
    REQUIRE(lim.b2 == 15);
    // Lemma 3.5 chain: k, then <= 1 + 2l/B1, then <= 1 + 2l/B2 (l = non-loop edges).
    const std::size_t l = r.passes[0].edges_in;
    REQUIRE(r.passes[0].vertices_out * lim.b1 <= lim.b1 + 2 * l);
    REQUIRE(r.passes[1].vertices_out * lim.b2 <= lim.b2 + 2 * l);
    REQUIRE(r.passes[2].vertices_out == 1);
  }
}

TEMPLATE_TEST_CASE("bottleneck property holds on seeded random graphs", "[paper][mst][property]",
                   std::uint64_t, double) {
  using T = TestType;
  for (std::uint64_t seed = 500; seed < 540; ++seed) {
    const std::size_t n = 2 + seed % 45;
    const auto edges =
        mst_test::random_edges<T>(n, (seed % 4) * n, seed, seed % 2 == 0 ? 3 : 10'000, true);
    const auto tree = mst_three_pass<T>(n, edges).edges;
    const auto minimax = mst_test::minimax_all_pairs<T>(n, edges);
    const auto bt = mst_test::tree_bottlenecks<T>(n, edges, tree);
    for (std::size_t x = 0; x < n; ++x) {
      for (std::size_t y = 0; y < n; ++y) {
        if (x == y) continue;
        REQUIRE(bt.has[x][y]);
        REQUIRE(bt.value[x][y] == minimax.value[x][y]);
      }
    }
  }
}
