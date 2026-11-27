// Tests for the entry point of Algorithm 4 (LocalHierarchy, l.2): the size test k <= 2 L_k^4
// (Sec 2 l.284-290) with saturating integer arithmetic, the base-case branch with its
// FallbackRecord (project fallback policy; reason string per open question Q-10), and the
// parameterised threshold used by Mode B (open question Q-04: one implementation, Mode A fixes
// the paper's constants).

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/base_case.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "paper/model/params.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::Vertex;
using sssp::paper::kPaperLocalParams;
using sssp::paper::local_construction_applies;
using sssp::paper::local_hierarchy;
using sssp::paper::local_hierarchy_base_case;
using sssp::paper::LocalHierarchyParams;

namespace {

template <class W>
CsrGraph<W> path_graph(std::size_t n) {
  std::vector<Edge<W>> e;
  e.reserve(n);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    e.push_back({static_cast<Vertex>(i), static_cast<Vertex>(i + 1), static_cast<W>(1 + i % 7)});
  }
  return CsrGraph<W>::from_edges(n, std::move(e));
}

}  // namespace

TEST_CASE("local threshold uses the paper constants c = 2, e = 4 by default") {
  STATIC_REQUIRE(kPaperLocalParams.threshold_coeff == 2);
  STATIC_REQUIRE(kPaperLocalParams.threshold_exponent == 4);
  // L_N = 18 for 2^17 - 2 < N <= 2^18 - 2, and 2 * 18^4 = 209952 (Sec 2 l.284-286).
  CHECK_FALSE(local_construction_applies(209952, kPaperLocalParams));
  CHECK(local_construction_applies(209953, kPaperLocalParams));
  CHECK_FALSE(local_construction_applies(0, kPaperLocalParams));
  CHECK_FALSE(local_construction_applies(1, kPaperLocalParams));
  CHECK_FALSE(local_construction_applies(32, kPaperLocalParams));
  for (std::uint64_t n : {0ULL, 1ULL, 2ULL, 1000ULL, 209952ULL, 209953ULL, 1ULL << 40}) {
    CHECK(local_construction_applies(n, kPaperLocalParams) ==
          sssp::paper::model::local_nontrivial(n));
  }
}

TEST_CASE("local threshold never wraps around for large exponents or sizes") {
  constexpr auto kMax = std::numeric_limits<std::uint64_t>::max();
  // L_max = 65 and 2 * 65^4 is tiny compared to 2^64 - 1.
  CHECK(local_construction_applies(kMax, kPaperLocalParams));
  // 65^64 overflows every word; the capped evaluation must still say "not above".
  CHECK_FALSE(local_construction_applies(kMax, LocalHierarchyParams{1, 64}));
  CHECK_FALSE(local_construction_applies(kMax, LocalHierarchyParams{kMax, 1}));
  // 2^20 vertices: L = 21, 21^64 saturates; must be false, not a wrapped small number.
  CHECK_FALSE(local_construction_applies(1ULL << 20, LocalHierarchyParams{1, 64}));
}

TEST_CASE("local threshold is a parameter for Mode B (N > c * L_N^e)") {
  // c = 1, e = 1: N > L_N = ceil(log2(N + 2)). L_3 = 3, L_4 = 3.
  const LocalHierarchyParams p{1, 1};
  CHECK_FALSE(local_construction_applies(3, p));
  CHECK(local_construction_applies(4, p));
  // c = 0 gives N > 0.
  CHECK_FALSE(local_construction_applies(0, LocalHierarchyParams{0, 4}));
  CHECK(local_construction_applies(1, LocalHierarchyParams{0, 4}));
}

TEMPLATE_TEST_CASE("Alg 4 below the threshold runs the base case and records the fallback", "",
                   std::uint64_t, double) {
  using W = TestType;
  // A 5-cycle with a chord and a self-loop: m = 7 edges in the input.
  const std::vector<Edge<W>> edges{{0, 1, W{3}}, {1, 2, W{1}}, {2, 3, W{3}}, {3, 4, W{2}},
                                   {4, 0, W{5}}, {1, 3, W{4}}, {2, 2, W{9}}};
  const auto g = CsrGraph<W>::from_edges(5, edges);
  Diagnostics diag;
  const auto got = local_hierarchy(g, &diag);
  const auto want = local_hierarchy_base_case(g);

  REQUIRE(diag.fallbacks.size() == 1);
  CHECK(diag.fallbacks[0].reason == "local: k <= 2*L_k^4 (paper base case)");
  CHECK(diag.fallbacks[0].vertices == 5);
  CHECK(diag.fallbacks[0].edges == 7);

  // Same output as calling the base case directly.
  CHECK(got.native_entry == want.native_entry);
  CHECK(got.arc_edge == want.arc_edge);
  REQUIRE(got.hierarchy.nodes.size() == want.hierarchy.nodes.size());
  for (std::size_t i = 0; i < got.hierarchy.nodes.size(); ++i) {
    CHECK(got.hierarchy.nodes[i].children == want.hierarchy.nodes[i].children);
    CHECK(got.hierarchy.nodes[i].height == want.hierarchy.nodes[i].height);
    CHECK(got.hierarchy.nodes[i].links.size() == want.hierarchy.nodes[i].links.size());
  }
  CHECK(got.hierarchy.roots == want.hierarchy.roots);
  CHECK(got.hierarchy.c0 == want.hierarchy.c0);
}

TEST_CASE("Alg 4 records path counters and accepts a null Diagnostics") {
  const auto g = path_graph<std::uint64_t>(4);
  Diagnostics diag;
  static_cast<void>(local_hierarchy(g, &diag));
  static_cast<void>(local_hierarchy(g, &diag));
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["local_hierarchy_calls"] == 2);
  CHECK(diag.counters["local_hierarchy_base_case"] == 2);
  CHECK(diag.counters["local_hierarchy_radial"] == 0);
#endif
  CHECK(diag.fallbacks.size() == 2);
  CHECK_NOTHROW(local_hierarchy(g, nullptr));
  CHECK(local_hierarchy(g).hierarchy.num_states() == 4);
}

TEST_CASE("Alg 4 on a single vertex is a base case with one leaf") {
  const auto g = CsrGraph<std::uint64_t>::from_edges(1, {});
  Diagnostics diag;
  const auto r = local_hierarchy(g, &diag);
  CHECK(r.hierarchy.num_states() == 1);
  REQUIRE(diag.fallbacks.size() == 1);
  CHECK(diag.fallbacks[0].vertices == 1);
  CHECK(diag.fallbacks[0].edges == 0);
}

TEST_CASE("Alg 4 Mode B reason string names the threshold that was used") {
  const auto g = path_graph<std::uint64_t>(3);
  Diagnostics diag;
  static_cast<void>(local_hierarchy(g, &diag, LocalHierarchyParams{1, 1}));
  REQUIRE(diag.fallbacks.size() == 1);
  CHECK(diag.fallbacks[0].reason == "local: k <= 1*L_k^1 (paper base case)");
}

TEST_CASE("Alg 4 base case rejects zero weights (contracted by the caller)") {
  const auto g = CsrGraph<std::uint64_t>::from_edges(2, {{0, 1, 0}});
  CHECK_THROWS_AS(local_hierarchy(g), std::invalid_argument);
}

TEST_CASE("Alg 4 above the threshold does not take the base case") {
  // Mode A at the first size where the radial construction applies (k = 209953, b = L_k^4 =
  // 18^4), and Mode B with c = 1, e = 1. The radial branch runs Sec 5.1-5.4 and records no
  // fallback; it requires N >= 32 and h = b + 1 <= N (Sec 5.4 l.1206).
  SECTION("Mode A, k = 209953") {
    const auto g = path_graph<std::uint64_t>(209953);
    Diagnostics diag;
    const auto r = local_hierarchy(g, &diag);
    CHECK(diag.fallbacks.empty());
    CHECK(r.hierarchy.c0 == 2);
    REQUIRE(r.native_entry.size() == 209953);
    for (Vertex v = 0; v < 209953; v += 997) CHECK(r.hierarchy.projection[r.native_entry[v]] == v);
#if SSSP_DIAGNOSTICS
    CHECK(diag.counters["local_hierarchy_radial"] == 1);
    CHECK(diag.counters["local_hierarchy_base_case"] == 0);
    CHECK(diag.counters["local_assemble_calls"] == 1);
#endif
  }
  SECTION("Mode A, k = 209952 still takes the base case") {
    const auto g = path_graph<std::uint64_t>(209952);
    Diagnostics diag;
    const auto r = local_hierarchy(g, &diag);
    CHECK(r.hierarchy.num_states() == 209952);
    REQUIRE(diag.fallbacks.size() == 1);
    CHECK(diag.fallbacks[0].vertices == 209952);
    CHECK(diag.fallbacks[0].edges == 209951);
  }
  SECTION("Mode B, k = 4 with N > L_N violates N >= 32 of Sec 5.4") {
    const auto g = path_graph<double>(4);
    Diagnostics diag;
    CHECK_THROWS_AS(local_hierarchy(g, &diag, LocalHierarchyParams{1, 1}), std::invalid_argument);
    CHECK(diag.fallbacks.empty());
  }
  SECTION("Mode B, k = 40 with the paper's b = L_k^4 > k violates h <= N") {
    const auto g = path_graph<double>(40);
    CHECK_THROWS_AS(local_hierarchy(g, nullptr, LocalHierarchyParams{1, 1}), std::invalid_argument);
  }
}

TEMPLATE_TEST_CASE("Alg 4 radial branch with a separator override yields Def 4.1 (c0 = 2)", "",
                   std::uint64_t, double) {
  using W = TestType;
  // A 6 x 6 grid with varied weights, b = 3: several pivots, macros and corridors.
  std::vector<Edge<W>> e;
  for (Vertex v = 0; v < 36; ++v) {
    if (v % 6 != 5) e.push_back({v, v + 1, static_cast<W>(1 + (v * 7) % 13)});
    if (v + 6 < 36) e.push_back({v, v + 6, static_cast<W>(1 + (v * 5) % 11)});
  }
  const auto g = CsrGraph<W>::from_edges(36, e);
  Diagnostics diag;
  const LocalHierarchyParams p{1, 1, 3};
  CHECK_FALSE(p == kPaperLocalParams);
  const auto r = local_hierarchy(g, &diag, p);
  CHECK(diag.fallbacks.empty());
  CHECK(r.hierarchy.c0 == 2);
  const auto report = sssp::paper::validate_ordered_hierarchy(r.hierarchy, g);
  INFO(report.summary());
  CHECK(report.ok());
  REQUIRE(r.native_entry.size() == 36);
  for (Vertex v = 0; v < 36; ++v) CHECK(r.hierarchy.projection[r.native_entry[v]] == v);
  CHECK(r.arc_edge.size() == r.hierarchy.graph.arcs.size());
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["local_hierarchy_radial"] == 1);
  CHECK(diag.counters["owner_decomposition_calls"] == 1);
  CHECK(diag.counters["threshold_joins_calls"] == 1);
  CHECK(diag.counters["local_assemble_calls"] == 1);
  CHECK(diag.counters["local_hierarchy_states"] == r.hierarchy.num_states());
  CHECK(diag.counters["local_hierarchy_arcs"] == r.hierarchy.graph.arcs.size());
  CHECK(diag.counters["local_hierarchy_nodes"] == r.hierarchy.nodes.size());
#endif
  // Zero weights are rejected in the radial branch as well (Sec 5.1 l.962).
  auto z = e;
  z[0].weight = W{0};
  CHECK_THROWS_AS(local_hierarchy(CsrGraph<W>::from_edges(36, z), nullptr, p),
                  std::invalid_argument);
}

TEST_CASE("Alg 4 radial branch requires a connected input") {
  // Alg 4 assumes a connected graph (Sec 5.1 l.961); the base case tolerates more, but the
  // radial construction (MST, owner cells) does not.
  std::vector<Edge<std::uint64_t>> two;  // two paths of 20 vertices: N = 40 >= 32, b = 3
  for (Vertex v = 0; v + 1 < 40; ++v) {
    if (v != 19) two.push_back({v, v + 1, 1 + v % 3});
  }
  const auto g40 = CsrGraph<std::uint64_t>::from_edges(40, two);
  CHECK_THROWS_AS(local_hierarchy(g40, nullptr, LocalHierarchyParams{1, 1, 3}),
                  std::invalid_argument);
  const auto g = CsrGraph<std::uint64_t>::from_edges(5, {{0, 1, 1}, {2, 3, 1}, {3, 4, 2}});
  Diagnostics diag;
  const auto r = local_hierarchy(g, &diag);  // base case: one tree per component
  CHECK(r.hierarchy.roots.size() == 2);
  CHECK(diag.fallbacks.size() == 1);
}
