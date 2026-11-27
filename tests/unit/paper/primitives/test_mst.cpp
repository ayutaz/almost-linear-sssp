// Tests for the three-pass MST of Sec 3.2 (Lemma 3.5) and lightest_edge_per_pair.
#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "paper/primitives/mst.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "unit/paper/primitives/mst_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::EdgeId;
using sssp::paper::kUnboundedGrowth;
using sssp::paper::lightest_edge_per_pair;
using sssp::paper::mst_kruskal_ref;
using sssp::paper::mst_three_pass;
using sssp::paper::spanning_forest_three_pass;
using sssp::paper::three_pass_limits;

namespace {

template <class W>
W w(std::uint64_t r) {
  return mst_test::make_weight<W>(r);
}

}  // namespace

TEST_CASE("three_pass_limits follow Lemma 3.5: B1 = max(2, ceil(log2 L_k)), B2 = L_k",
          "[paper][mst]") {
  // L_k = ceil(log2(k + 2)), evaluated by hand.
  struct Row {
    std::size_t k, b1, b2;
  };
  const std::vector<Row> rows{
      {1, 2, 2},       // L = 2 (log2 3), ceil(log2 2) = 1 -> 2
      {2, 2, 2},       // L = 2 (log2 4)
      {3, 2, 3},       // L = 3, ceil(log2 3) = 2
      {14, 2, 4},      // L = 4 (log2 16)
      {15, 3, 5},      // L = 5 (log2 17), ceil(log2 5) = 3
      {254, 3, 8},     // L = 8 (log2 256)
      {255, 4, 9},     // L = 9, ceil(log2 9) = 4
      {65534, 4, 16},  // L = 16
      {65535, 5, 17},  // L = 17
      {std::size_t{1} << 32, 6, 33},
      {std::numeric_limits<std::size_t>::max(), 6, 64},  // k + 2 > 2^64 - 1: L = 64, no wrap
  };
  for (const auto& r : rows) {
    const auto lim = three_pass_limits(r.k);
    INFO("k = " << r.k);
    REQUIRE(lim.b1 == r.b1);
    REQUIRE(lim.b2 == r.b2);
  }
}

TEMPLATE_TEST_CASE("mst_three_pass equals the Kruskal MST on a hand example", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  // 0-1 (4), 1-2 (1), 2-3 (3), 3-0 (2), 0-2 (5), 1-3 (3), 4-0 (6), loop 2-2 (0), parallel 1-2 (1)
  const std::vector<Edge<T>> edges{{0, 1, w<T>(4)}, {1, 2, w<T>(1)}, {2, 3, w<T>(3)},
                                   {3, 0, w<T>(2)}, {0, 2, w<T>(5)}, {1, 3, w<T>(3)},
                                   {4, 0, w<T>(6)}, {2, 2, w<T>(0)}, {2, 1, w<T>(1)}};
  // Order: 7(loop) 1 8 3 2 5 0 4 6 -> picks 1, 3, 2, 6 (8 parallel to 1, 5 closes 1-2-3).
  const std::vector<EdgeId> expected{1, 2, 3, 6};
  const auto r = mst_three_pass<T>(5, edges);
  REQUIRE(r.edges == expected);
  REQUIRE(mst_kruskal_ref<T>(5, edges) == expected);
}

TEMPLATE_TEST_CASE("mst_three_pass breaks equal lengths by input index", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  // Square 0-1-2-3-0 plus diagonal 0-2, all of length 7: the MST is the three lowest indices
  // that do not close a cycle.
  const std::vector<Edge<T>> edges{
      {0, 1, w<T>(7)}, {1, 2, w<T>(7)}, {2, 3, w<T>(7)}, {3, 0, w<T>(7)}, {0, 2, w<T>(7)}};
  REQUIRE(mst_three_pass<T>(4, edges).edges == std::vector<EdgeId>{0, 1, 2});
  // Reversing the indices changes the tree accordingly.
  const std::vector<Edge<T>> rev{
      {0, 2, w<T>(7)}, {3, 0, w<T>(7)}, {2, 3, w<T>(7)}, {1, 2, w<T>(7)}, {0, 1, w<T>(7)}};
  REQUIRE(mst_three_pass<T>(4, rev).edges == std::vector<EdgeId>{0, 1, 3});
  // Zero lengths are ordinary lengths.
  const std::vector<Edge<T>> zeros{{0, 1, w<T>(0)}, {1, 2, w<T>(0)}, {0, 2, w<T>(0)}};
  REQUIRE(mst_three_pass<T>(3, zeros).edges == std::vector<EdgeId>{0, 1});
}

TEMPLATE_TEST_CASE("mst_three_pass handles a single vertex, loops and multi-edges", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  const auto single = mst_three_pass<T>(1, std::vector<Edge<T>>{});
  REQUIRE(single.edges.empty());
  REQUIRE(single.passes.empty());  // "A singleton requires no pass"
  const auto loops = mst_three_pass<T>(1, std::vector<Edge<T>>{{0, 0, w<T>(1)}, {0, 0, w<T>(0)}});
  REQUIRE(loops.edges.empty());
  const std::vector<Edge<T>> multi{
      {0, 1, w<T>(9)}, {1, 0, w<T>(2)}, {0, 1, w<T>(2)}, {1, 1, w<T>(0)}, {0, 1, w<T>(3)}};
  REQUIRE(mst_three_pass<T>(2, multi).edges == std::vector<EdgeId>{1});
}

TEMPLATE_TEST_CASE("mst_three_pass requires a connected input; the forest variant does not",
                   "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  // Components {0,1,2}, {3}, {4,5}.
  const std::vector<Edge<T>> edges{
      {0, 1, w<T>(3)}, {1, 2, w<T>(1)}, {2, 0, w<T>(2)}, {4, 5, w<T>(8)}, {3, 3, w<T>(1)}};
  REQUIRE_THROWS_AS(mst_three_pass<T>(6, edges), std::invalid_argument);
  REQUIRE_THROWS_AS(mst_three_pass<T>(2, std::vector<Edge<T>>{}), std::invalid_argument);
  REQUIRE_THROWS_AS(mst_three_pass<T>(0, std::vector<Edge<T>>{}), std::invalid_argument);
  const auto forest = spanning_forest_three_pass<T>(6, edges);
  REQUIRE(forest.edges == std::vector<EdgeId>{1, 2, 3});
  REQUIRE(forest.edges == mst_kruskal_ref<T>(6, edges));
  REQUIRE(spanning_forest_three_pass<T>(0, std::vector<Edge<T>>{}).edges.empty());
  REQUIRE(spanning_forest_three_pass<T>(3, std::vector<Edge<T>>{}).edges.empty());
}

TEMPLATE_TEST_CASE("mst_three_pass validates edges", "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  REQUIRE_THROWS_AS(mst_three_pass<T>(2, std::vector<Edge<T>>{{0, 2, w<T>(1)}}), std::out_of_range);
  if constexpr (std::is_floating_point_v<T>) {
    REQUIRE_THROWS_AS(mst_three_pass<T>(2, std::vector<Edge<T>>{{0, 1, T{-1}}}),
                      std::invalid_argument);
  }
}

TEMPLATE_TEST_CASE("mst_three_pass equals Kruskal and stays within the per-pass bounds",
                   "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  bool saw_three_passes = false;
  for (std::uint64_t seed = 1; seed <= 60; ++seed) {
    const std::size_t n = 1 + (seed * 37) % 300;
    const std::size_t extra = (seed % 5) * n;
    const auto edges = mst_test::random_edges<T>(n, extra, seed, seed % 2 == 0 ? 3 : 1000, true);
    const auto r = mst_three_pass<T>(n, edges);
    INFO("seed " << seed << " n " << n);
    REQUIRE(r.edges == mst_kruskal_ref<T>(n, edges));
    REQUIRE(r.edges.size() == n - 1);
    const auto lim = three_pass_limits(n);
    const std::vector<std::size_t> limits{lim.b1, lim.b2, kUnboundedGrowth};
    REQUIRE(r.passes.size() <= 3);
    std::size_t vertices = n;
    for (std::size_t i = 0; i < r.passes.size(); ++i) {
      const auto& p = r.passes[i];
      REQUIRE(p.limit == limits[i]);
      REQUIRE(p.vertices_in == vertices);
      if (p.limit != kUnboundedGrowth) {
        // Lemma 3.4 for a connected quotient: at most 1 + 2 l'/B vertices.
        REQUIRE(p.vertices_out * p.limit <= p.limit + 2 * p.edges_in);
      }
      vertices = p.vertices_out;
    }
    REQUIRE(vertices == 1);
    if (r.passes.size() == 3) saw_three_passes = true;
  }
  REQUIRE(saw_three_passes);  // the inputs exercise every pass
}

TEMPLATE_TEST_CASE("MST tree paths realise the minimax walk value (Lemma 3.5 bottleneck)",
                   "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  for (std::uint64_t seed = 100; seed < 130; ++seed) {
    const std::size_t n = 2 + seed % 25;
    const auto edges = mst_test::random_edges<T>(n, 2 * n, seed, seed % 3 == 0 ? 2 : 40, true);
    const auto tree = mst_three_pass<T>(n, edges).edges;
    const auto minimax = mst_test::minimax_all_pairs<T>(n, edges);
    const auto bt = mst_test::tree_bottlenecks<T>(n, edges, tree);
    for (std::size_t x = 0; x < n; ++x) {
      for (std::size_t y = 0; y < n; ++y) {
        if (x == y) continue;
        REQUIRE(bt.has[x][y]);
        REQUIRE(minimax.has[x][y]);
        REQUIRE(bt.value[x][y] == minimax.value[x][y]);
      }
    }
    for (const auto& e : edges) {
      if (e.u == e.v) continue;
      REQUIRE_FALSE(e.weight < bt.value[e.u][e.v]);  // b_T(x, y) <= w(xy)
    }
  }
}

TEMPLATE_TEST_CASE("mst_three_pass records per-pass quotient sizes in Diagnostics", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  const std::size_t n = 400;
  const auto edges = mst_test::random_edges<T>(n, 3 * n, 7, 1000, true);
  Diagnostics diag;
  const auto r = mst_three_pass<T>(n, edges, &diag);
  REQUIRE(diag.counters.at("mst.passes") == r.passes.size());
  std::uint64_t pops = 0;
  for (std::size_t i = 0; i < r.passes.size(); ++i) {
    const std::string p = "mst.pass" + std::to_string(i + 1) + ".";
    REQUIRE(diag.counters.at(p + "vertices_in") == r.passes[i].vertices_in);
    REQUIRE(diag.counters.at(p + "edges_in") == r.passes[i].edges_in);
    REQUIRE(diag.counters.at(p + "vertices_out") == r.passes[i].vertices_out);
    REQUIRE(diag.counters.at(p + "growths") == r.passes[i].growths);
    pops += r.passes[i].extractions;
  }
  REQUIRE(pops == n - 1);  // one extraction per selected edge
  REQUIRE(diag.ops.queue_pops == pops);
#if SSSP_DIAGNOSTICS
  REQUIRE(diag.ops.comparisons > 0);
#else
  REQUIRE(diag.ops.comparisons == 0);  // comparison counting is compiled out in timed builds
#endif
  REQUIRE(diag.fallbacks.empty());  // no direct-search substitute anywhere
}

TEMPLATE_TEST_CASE("mst_three_pass accepts a CsrGraph", "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  auto edges = mst_test::random_edges<T>(50, 100, 3, 5, true);
  const auto g = CsrGraph<T>::from_edges(50, edges);
  REQUIRE(mst_three_pass(g).edges == mst_kruskal_ref<T>(50, edges));
  REQUIRE(spanning_forest_three_pass(g).edges == mst_kruskal_ref<T>(50, edges));
}

TEMPLATE_TEST_CASE("lightest_edge_per_pair keeps the (length, index)-least edge per pair",
                   "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  const std::vector<Edge<T>> edges{{0, 1, w<T>(5)}, {1, 0, w<T>(3)}, {0, 1, w<T>(3)},
                                   {2, 2, w<T>(1)}, {1, 2, w<T>(4)}, {2, 1, w<T>(4)},
                                   {0, 2, w<T>(9)}};
  REQUIRE(lightest_edge_per_pair<T>(3, edges) == std::vector<EdgeId>{1, 4, 6});
  REQUIRE(lightest_edge_per_pair<T>(0, std::vector<Edge<T>>{}).empty());
  REQUIRE(lightest_edge_per_pair<T>(1, std::vector<Edge<T>>{{0, 0, w<T>(1)}}).empty());
  REQUIRE_THROWS_AS(lightest_edge_per_pair<T>(1, std::vector<Edge<T>>{{0, 1, w<T>(1)}}),
                    std::out_of_range);
  // On random multigraphs: one survivor per unordered pair, it is the least, and the MST is
  // unchanged by dropping the others (Alg 4 l.4 relies on this).
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    const std::size_t n = 2 + seed % 12;
    const auto g = mst_test::random_edges<T>(n, 4 * n, seed, 4, true);
    const auto kept = lightest_edge_per_pair<T>(n, g);
    std::vector<std::vector<int>> count(n, std::vector<int>(n, 0));
    std::vector<Edge<T>> sub;
    for (const auto id : kept) {
      const auto& e = g[id];
      REQUIRE(e.u != e.v);
      ++count[std::min(e.u, e.v)][std::max(e.u, e.v)];
      for (std::size_t j = 0; j < g.size(); ++j) {
        const auto& f = g[j];
        const bool same = (f.u == e.u && f.v == e.v) || (f.u == e.v && f.v == e.u);
        if (same) {
          REQUIRE_FALSE((sssp::paper::MstEdgeKey<T>{f.weight, static_cast<EdgeId>(j)} <
                         sssp::paper::MstEdgeKey<T>{e.weight, id}));
        }
      }
    }
    for (const auto& e : g) {
      if (e.u != e.v) REQUIRE(count[std::min(e.u, e.v)][std::max(e.u, e.v)] == 1);
    }
    const auto full = mst_three_pass<T>(n, g).edges;
    for (const auto id : full) {
      REQUIRE(std::find(kept.begin(), kept.end(), id) != kept.end());
    }
  }
}
