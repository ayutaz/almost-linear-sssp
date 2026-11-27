// Outer base case of Algorithm 8 (arXiv 2611.34911v1, Sec 8.3, Alg 8 l.4-5; Sec 2 txt l.286-288;
// Lemma 3.2): retain a lightest edge per vertex pair, then compute exact source distances by
// Dijkstra's search. `outer_base_case` uses the scan version of Lemma 3.2 (paper-faithful),
// `outer_base_case_d1` the binary-heap version of the same lemma (deviation D-1, ADR 0008).
// Both record a FallbackRecord with the ADR 0008 reason and the size of the graph they solved.
//
// Expected values are computed by hand or by independent oracles (Bellman–Ford, baseline
// binary-heap Dijkstra, verify_sssp), never copied from an implementation run.
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/outer/base_case.hpp"
#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::SsspResult;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::outer::kOuterBaseCaseD1Reason;
using sssp::paper::outer::kOuterBaseCaseScanReason;
using sssp::paper::outer::outer_base_case;
using sssp::paper::outer::outer_base_case_d1;

namespace {

enum class Variant { kScan, kD1 };

template <class W>
SsspResult<W> run(Variant v, const CsrGraph<W>& g, Vertex s) {
  return v == Variant::kScan ? outer_base_case(g, s) : outer_base_case_d1(g, s);
}

std::string expected_reason(Variant v) {
  return v == Variant::kScan ? "outer: N <= L_N^64 (paper base case, scan Dijkstra)"
                             : "outer: N <= L_N^64 (paper base case, binary-heap Dijkstra, "
                               "deviation D-1)";
}

template <class W>
CsrGraph<W> make(std::size_t n, std::vector<Edge<W>> e) {
  return CsrGraph<W>::from_edges(n, std::move(e));
}

template <class W>
void require_single_record(const SsspResult<W>& r, Variant v, std::uint64_t n, std::uint64_t m) {
  REQUIRE(r.diag.fallbacks.size() == 1);
  CHECK(r.diag.fallbacks[0].reason == expected_reason(v));
  CHECK(r.diag.fallbacks[0].vertices == n);
  CHECK(r.diag.fallbacks[0].edges == m);
  CHECK(r.diag.counters.at("outer_base_cases") == 1);
}

constexpr Variant kVariants[] = {Variant::kScan, Variant::kD1};

}  // namespace

TEST_CASE("outer base case reasons are exactly those of ADR 0008", "[paper][outer][base_case]") {
  CHECK(std::string(kOuterBaseCaseScanReason) == expected_reason(Variant::kScan));
  CHECK(std::string(kOuterBaseCaseD1Reason) == expected_reason(Variant::kD1));
}

TEMPLATE_TEST_CASE("outer base case: hand-computed distances and tree", "[paper][outer][base_case]",
                   std::uint64_t, double) {
  using W = TestType;
  // e0 0-1 (4), e1 0-2 (1), e2 2-1 (2), e3 1-3 (5), e4 2-3 (8), e5 3-4 (3).
  // d = {0, 3, 1, 8, 11}; every shortest path is unique: 1 <- 2 (e2), 2 <- 0 (e1),
  // 3 <- 1 (e3), 4 <- 3 (e5).
  const auto g = make<W>(
      5, {{0, 1, W{4}}, {0, 2, W{1}}, {2, 1, W{2}}, {1, 3, W{5}}, {2, 3, W{8}}, {3, 4, W{3}}});
  for (const Variant v : kVariants) {
    const auto r = run(v, g, 0);
    CHECK(r.dist == std::vector<W>{W{0}, W{3}, W{1}, W{8}, W{11}});
    CHECK(r.parent == std::vector<Vertex>{kNoVertex, 2, 0, 1, 3});
    CHECK(r.parent_edge == std::vector<EdgeId>{kNoEdge, 2, 1, 3, 5});
    require_single_record(r, v, 5, 6);
    REQUIRE(sssp::verify_sssp(g, 0, r).ok);
  }
  // Same graph from vertex 4: d = {11, 8, 10, 3, 0}.
  for (const Variant v : kVariants) {
    const auto r = run(v, g, 4);
    CHECK(r.dist == std::vector<W>{W{11}, W{8}, W{10}, W{3}, W{0}});
    CHECK(r.parent == std::vector<Vertex>{2, 3, 1, 4, kNoVertex});
  }
}

TEMPLATE_TEST_CASE("outer base case: one vertex", "[paper][outer][base_case]", std::uint64_t,
                   double) {
  using W = TestType;
  const auto g = make<W>(1, {});
  for (const Variant v : kVariants) {
    const auto r = run(v, g, 0);
    CHECK(r.dist == std::vector<W>{W{0}});
    CHECK(r.parent == std::vector<Vertex>{kNoVertex});
    CHECK(r.parent_edge == std::vector<EdgeId>{kNoEdge});
    require_single_record(r, v, 1, 0);
  }
}

TEMPLATE_TEST_CASE("outer base case: unreachable vertices get infinity",
                   "[paper][outer][base_case]", std::uint64_t, double) {
  using W = TestType;
  // Components {0,1,2}, {3,4}, {5}.
  const auto g = make<W>(6, {{0, 1, W{2}}, {1, 2, W{2}}, {3, 4, W{1}}});
  const W inf = WeightTraits<W>::infinity();
  for (const Variant v : kVariants) {
    const auto r = run(v, g, 1);
    CHECK(r.dist == std::vector<W>{W{2}, W{0}, W{2}, inf, inf, inf});
    CHECK(r.parent == std::vector<Vertex>{1, kNoVertex, 1, kNoVertex, kNoVertex, kNoVertex});
    CHECK(r.parent_edge == std::vector<EdgeId>{0, kNoEdge, 1, kNoEdge, kNoEdge, kNoEdge});
    require_single_record(r, v, 6, 3);
    REQUIRE(sssp::verify_sssp(g, 1, r).ok);
  }
}

TEMPLATE_TEST_CASE("outer base case: zero weights are non-negative lengths",
                   "[paper][outer][base_case]", std::uint64_t, double) {
  using W = TestType;
  // Lemma 3.2 holds for non-negative lengths. 0 -0- 1 -3- 2 -0- 3, and 0 -5- 3.
  const auto g = make<W>(4, {{0, 1, W{0}}, {1, 2, W{3}}, {2, 3, W{0}}, {0, 3, W{5}}});
  for (const Variant v : kVariants) {
    const auto r = run(v, g, 0);
    CHECK(r.dist == std::vector<W>{W{0}, W{0}, W{3}, W{3}});
    CHECK(r.parent[1] == 0);
    CHECK(r.parent[2] == 1);
    CHECK(r.parent[3] == 2);
    REQUIRE(sssp::verify_sssp(g, 0, r).ok);
  }
}

TEMPLATE_TEST_CASE("outer base case: lightest edge per vertex pair (Alg 8 l.5)",
                   "[paper][outer][base_case]", std::uint64_t, double) {
  using W = TestType;
  // Pair {0,1}: e0 (7), e1 (3), e2 (3) -> e1 is (weight, id)-least. Loop e3 is dropped.
  // Pair {1,2}: e4 (2), e5 (1) -> e5. Kept: e1, e5; removed: e0, e2, e3, e4.
  const auto g = make<W>(
      3, {{0, 1, W{7}}, {0, 1, W{3}}, {1, 0, W{3}}, {1, 1, W{1}}, {1, 2, W{2}}, {2, 1, W{1}}});
  for (const Variant v : kVariants) {
    const auto r = run(v, g, 0);
    CHECK(r.dist == std::vector<W>{W{0}, W{3}, W{4}});
    CHECK(r.parent == std::vector<Vertex>{kNoVertex, 0, 1});
    CHECK(r.parent_edge == std::vector<EdgeId>{kNoEdge, 1, 5});
    CHECK(r.diag.counters.at("outer_base_pair_edges_kept") == 2);
    CHECK(r.diag.counters.at("outer_base_pair_edges_removed") == 4);
    // The record states the size of the graph handed to the base case (before l.5).
    require_single_record(r, v, 3, 6);
    REQUIRE(sssp::verify_sssp(g, 0, r).ok);
  }
  // From vertex 2 the lighter parallel edge with the larger id is chosen: 1 <- 2 via e5,
  // 0 <- 1 via e1.
  for (const Variant v : kVariants) {
    const auto r = run(v, g, 2);
    CHECK(r.dist == std::vector<W>{W{4}, W{1}, W{0}});
    CHECK(r.parent_edge == std::vector<EdgeId>{1, 5, kNoEdge});
  }
}

TEST_CASE("outer base case rejects an invalid source", "[paper][outer][base_case]") {
  const auto g = make<std::uint64_t>(3, {{0, 1, 1}});
  REQUIRE_THROWS_AS(outer_base_case(g, 3), std::out_of_range);
  REQUIRE_THROWS_AS(outer_base_case_d1(g, 3), std::out_of_range);
  const auto empty = make<double>(0, {});
  REQUIRE_THROWS_AS(outer_base_case(empty, 0), std::out_of_range);
  REQUIRE_THROWS_AS(outer_base_case_d1(empty, 0), std::out_of_range);
}

TEST_CASE("outer base case reports integer distance overflow", "[paper][outer][base_case]") {
  const std::uint64_t big = WeightTraits<std::uint64_t>::infinity() / 2 + 1;
  const auto g = make<std::uint64_t>(3, {{0, 1, big}, {1, 2, big}});
  REQUIRE_THROWS_AS(outer_base_case(g, 0), std::overflow_error);
  REQUIRE_THROWS_AS(outer_base_case_d1(g, 0), std::overflow_error);
}

#if SSSP_DIAGNOSTICS  // operation counts are compiled out when SSSP_DIAGNOSTICS=0
TEST_CASE("scan base case spends order N^2 comparisons, the D-1 heap order N log N",
          "[paper][outer][base_case][lemma3.2]") {
  // Star with centre 0 and k leaves at distinct weights: after the centre is settled all k leaves
  // are pending, so the scans of Lemma 3.2 cost (k-1) + (k-2) + ... + 0 = k(k-1)/2 comparisons.
  const std::uint64_t k = 2000;
  std::vector<Edge<std::uint64_t>> e;
  for (Vertex v = 1; v <= k; ++v) e.push_back({0, v, (v * 7919) % 1000 + 1});
  const auto g = make<std::uint64_t>(k + 1, std::move(e));
  const std::uint64_t n = k + 1;
  const std::uint64_t l = g.num_arcs();
  const auto scan = outer_base_case(g, 0);
  const auto heap = outer_base_case_d1(g, 0);
  REQUIRE(scan.dist == heap.dist);
  const std::uint64_t scan_cmp = scan.diag.ops.comparisons;
  const std::uint64_t heap_cmp = heap.diag.ops.comparisons;
  CHECK(scan_cmp >= k * (k - 1) / 2);
  CHECK(scan_cmp <= n * n + l + g.num_edges());  // Lemma 3.2: O(N^2 + l)
  // Binary heap: O(1 + log(N + 1)) per operation, N pushes and N pops (log2(2001) < 11).
  CHECK(heap_cmp <= 4 * n * 11 + l + g.num_edges());
  CHECK(heap_cmp * 10 < scan_cmp);
  CHECK(scan.diag.ops.vertex_scans == n);
  CHECK(heap.diag.ops.vertex_scans == n);
}
#endif

TEMPLATE_TEST_CASE("outer base case equals Bellman-Ford on every generator family at small N",
                   "[paper][outer][base_case][differential]", std::uint64_t, double) {
  using W = TestType;
  std::size_t runs = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    for (std::uint64_t seed = 1; seed <= 3; ++seed) {
      sssp::gen::GraphSpec spec;
      spec.type = std::string(info.name);
      spec.n = 120;
      spec.avg_degree = 3.0;
      spec.degree = 3;
      spec.multi = sssp::gen::MultiEdgePolicy::kKeep;
      spec.options.seed = seed;
      if constexpr (WeightTraits<W>::is_exact) {
        spec.options.weights.dist = sssp::gen::WeightDist::kUniformInt;
      } else {
        spec.options.weights.dist = sssp::gen::WeightDist::kUniformReal;
      }
      spec.options.weights.lo = 1.0;
      spec.options.weights.hi = 50.0;
      const auto g = sssp::to_csr(sssp::gen::generate<W>(spec));
      const auto bf = sssp::baseline::bellman_ford(g, 0);
      const auto tol = sssp::Tolerance::for_graph(g);
      for (const Variant v : kVariants) {
        INFO(info.name << " seed=" << seed << " variant=" << expected_reason(v));
        const auto r = run(v, g, 0);
        const auto rep = sssp::verify_sssp(g, 0, r, &bf.dist, tol);
        INFO(rep.summary());
        REQUIRE(rep.ok);
        if constexpr (WeightTraits<W>::is_exact) REQUIRE(r.dist == bf.dist);
        require_single_record(r, v, g.num_vertices(), g.num_edges());
        ++runs;
      }
    }
  }
  REQUIRE(runs == sssp::gen::generator_catalog().size() * 3 * 2);
}

TEST_CASE("scan base case at N = 2000 and D-1 base case at N = 90000 equal baseline Dijkstra",
          "[paper][outer][base_case][differential]") {
  for (const std::size_t side : {std::size_t{44}, std::size_t{300}}) {
    sssp::gen::GraphSpec spec;
    spec.type = "grid";
    spec.rows = side;
    spec.cols = side;
    spec.n = side * side;
    spec.options.seed = 7;
    spec.options.weights.dist = sssp::gen::WeightDist::kUniformInt;
    spec.options.weights.lo = 1.0;
    spec.options.weights.hi = 1000.0;
    const auto g = sssp::to_csr(sssp::gen::generate<std::uint64_t>(spec));
    const auto s = static_cast<Vertex>(g.num_vertices() / 3);
    const auto ref = sssp::baseline::dijkstra_binary_heap(g, s);
    INFO("N=" << g.num_vertices());
    const auto d1 = outer_base_case_d1(g, s);
    REQUIRE(d1.dist == ref.dist);
    REQUIRE(sssp::verify_sssp(g, s, d1, &ref.dist).ok);
    require_single_record(d1, Variant::kD1, g.num_vertices(), g.num_edges());
    if (g.num_vertices() <= 2000) {
      const auto scan = outer_base_case(g, s);
      REQUIRE(scan.dist == ref.dist);
      REQUIRE(sssp::verify_sssp(g, s, scan, &ref.dist).ok);
      require_single_record(scan, Variant::kScan, g.num_vertices(), g.num_edges());
    }
  }
}
