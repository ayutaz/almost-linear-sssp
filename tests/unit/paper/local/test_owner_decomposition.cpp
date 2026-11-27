// Tests for the owner decomposition of Sec 5.1 (Def 5.1, Lemmas 5.2-5.3) on hand-computed
// examples. The separator size b is a parameter (local: L_N^4; outer: its own b); the small b
// used here make the construction nontrivial on small graphs, as Mode B does (Q-03).
//
// Paper-unspecified choices tested here (open question Q-31): T is rooted at vertex 0 for the
// postorder marking; with no mark the root is the single pivot; among maximum corridor edges
// the smallest edge id is cut.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::EdgeId;
using sssp::Vertex;
using sssp::paper::kNoCorridor;
using sssp::paper::kNoPivot;
using sssp::paper::local_separator_size;
using sssp::paper::owner_decomposition;

namespace {

template <class W>
std::vector<Vertex> owner_vertices(const sssp::paper::OwnerDecomposition<W>& od) {
  std::vector<Vertex> out;
  for (const auto p : od.owner) out.push_back(od.pivots[p]);
  return out;
}

}  // namespace

TEST_CASE("local separator size is L_N^4") {
  CHECK(local_separator_size(1) == 16);           // L_1 = 2
  CHECK(local_separator_size(6) == 81);           // L_6 = 3
  CHECK(local_separator_size(209953) == 104976);  // L = 18
}

TEMPLATE_TEST_CASE("owner decomposition of a path with b = 3", "", std::uint64_t, double) {
  using W = TestType;
  // 0 -5- 1 -2- 2 -7- 3 -7- 4 -1- 5 -3- 6; edge i joins i and i+1.
  const std::vector<Edge<W>> edges{{0, 1, W{5}}, {1, 2, W{2}}, {2, 3, W{7}},
                                   {3, 4, W{7}}, {4, 5, W{1}}, {5, 6, W{3}}};
  const auto g = CsrGraph<W>::from_edges(7, edges);
  const auto od = owner_decomposition(g, 3);

  CHECK(od.separator == 3);
  CHECK(od.hop_bound == 4);
  CHECK(od.tree_edges == std::vector<EdgeId>{0, 1, 2, 3, 4, 5});
  // Postorder from root 0: 6 (1), 5 (2), 4 (3: mark), 3 (1), 2 (2), 1 (3: mark), 0 (1).
  CHECK(od.marked == std::vector<std::uint8_t>{0, 1, 0, 0, 1, 0, 0});
  CHECK(od.pivots == std::vector<Vertex>{1, 4});
  CHECK(od.pivot_index[1] == 0);
  CHECK(od.pivot_index[4] == 1);
  CHECK(od.pivot_index[0] == kNoPivot);
  CHECK(od.is_pivot(4));
  CHECK_FALSE(od.is_pivot(3));
  // Corridor 1-2-3-4 has maximum 7 on edges 2 and 3; the smaller id (2) is cut (Q-31).
  REQUIRE(od.corridors.size() == 1);
  CHECK(od.corridors[0].edge == 2);
  CHECK(od.corridors[0].weight == W{7});
  CHECK(od.cut == std::vector<std::uint8_t>{0, 0, 1, 0, 0, 0});
  CHECK(owner_vertices(od) == std::vector<Vertex>{1, 1, 1, 4, 4, 4, 4});
  CHECK(od.radius == std::vector<W>{W{5}, W{0}, W{2}, W{7}, W{0}, W{1}, W{3}});
  CHECK(od.owner_distance == std::vector<W>{W{5}, W{0}, W{2}, W{7}, W{0}, W{1}, W{4}});
  CHECK(od.owner_hops == std::vector<std::uint32_t>{1, 0, 1, 1, 0, 1, 2});
  CHECK(od.sigma(0, 1) == W{7});
  CHECK(od.sigma(1, 0) == W{7});
  CHECK(od.sigma(0, 0) == W{0});
  CHECK(od.corridor_max(0, 1) == 0);
  CHECK(od.corridor_max(1, 1) == kNoCorridor);
}

TEST_CASE("owner decomposition adds the branch vertex of the marks' Steiner tree") {
  using W = std::uint64_t;
  // Spider: centre 0 with legs 0-1-2, 0-3-4, 0-5-6.
  const std::vector<Edge<W>> edges{{0, 1, 4}, {1, 2, 1}, {0, 3, 6},
                                   {3, 4, 1}, {0, 5, 2}, {5, 6, 1}};
  const auto g = CsrGraph<W>::from_edges(7, edges);
  SECTION("b = 2: marks 1, 3, 5 and branch vertex 0") {
    const auto od = owner_decomposition(g, 2);
    CHECK(od.marked == std::vector<std::uint8_t>{0, 1, 0, 1, 0, 1, 0});
    CHECK(od.pivots == std::vector<Vertex>{0, 1, 3, 5});
    // Adjacent pivots: their common edge is cut.
    REQUIRE(od.corridors.size() == 3);
    CHECK(od.cut == std::vector<std::uint8_t>{1, 0, 1, 0, 1, 0});
    CHECK(owner_vertices(od) == std::vector<Vertex>{0, 1, 1, 3, 3, 5, 5});
    CHECK(od.radius == std::vector<W>{0, 0, 1, 0, 1, 0, 1});
    // T_A is the star 0-1 (4), 0-3 (6), 0-5 (2); sigma(1, 3) = 6, sigma(1, 5) = 4.
    const auto idx = [&](Vertex v) { return od.pivot_index[v]; };
    CHECK(od.sigma(idx(1), idx(3)) == 6);
    CHECK(od.sigma(idx(1), idx(5)) == 4);
    CHECK(od.sigma(idx(5), idx(0)) == 2);
    CHECK(od.corridors[od.corridor_max(idx(1), idx(3))].edge == 2);
  }
  SECTION("b = 3: only the centre is marked") {
    const auto od = owner_decomposition(g, 3);
    CHECK(od.pivots == std::vector<Vertex>{0});
    CHECK(od.corridors.empty());
    CHECK(od.radius == std::vector<W>{0, 4, 4, 6, 6, 2, 2});
    CHECK(od.owner_distance == std::vector<W>{0, 4, 5, 6, 7, 2, 3});
  }
}

TEST_CASE("owner decomposition with b = 1 makes every vertex a pivot") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(4, {{0, 1, 3}, {1, 2, 1}, {2, 3, 2}, {0, 3, 9}});
  const auto od = owner_decomposition(g, 1);
  CHECK(od.pivots == std::vector<Vertex>{0, 1, 2, 3});
  CHECK(od.corridors.size() == 3);
  CHECK(od.radius == std::vector<W>{0, 0, 0, 0});
  CHECK(od.sigma(0, 3) == 3);
}

TEST_CASE("owner decomposition with b > N has the root as its single pivot") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(4, {{2, 1, 3}, {1, 0, 1}, {3, 2, 2}});
  Diagnostics diag;
  const auto od = owner_decomposition(g, 5, &diag);
  CHECK(od.marked == std::vector<std::uint8_t>{1, 0, 0, 0});
  CHECK(od.pivots == std::vector<Vertex>{0});
  CHECK(od.radius == std::vector<W>{0, 1, 3, 3});
  CHECK(od.owner_distance == std::vector<W>{0, 1, 4, 6});
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["owner_decomposition_calls"] == 1);
  CHECK(diag.counters["owner_decomposition_pivots"] == 1);
  CHECK(diag.counters["owner_decomposition_unmarked"] == 1);
#endif
}

TEST_CASE("owner decomposition on a single vertex") {
  const auto g = CsrGraph<std::uint64_t>::from_edges(1, {{0, 0, 4}});
  for (const std::size_t b : {1, 2}) {
    const auto od = owner_decomposition(g, b);
    CHECK(od.pivots == std::vector<Vertex>{0});
    CHECK(od.corridors.empty());
    CHECK(od.reduced == std::vector<std::uint8_t>{0});  // a loop is heavier than its tree path
    CHECK(od.sigma(0, 0) == 0);
  }
}

TEST_CASE("reduced graph omits nontree edges heavier than their tree path") {
  using W = std::uint64_t;
  // MST {0-1 (1), 1-2 (1)}; D_T(0, 2) = 2.
  const std::vector<Edge<W>> edges{
      {0, 1, 1},  // 0: tree
      {1, 2, 1},  // 1: tree
      {0, 2, 5},  // 2: 5 > 2, deleted
      {0, 2, 2},  // 3: 2 = D_T, kept
      {0, 1, 1},  // 4: parallel, equal to the tree edge, kept
      {1, 0, 3},  // 5: parallel, heavier, deleted
      {2, 2, 1},  // 6: loop, deleted
  };
  const auto od = owner_decomposition(CsrGraph<W>::from_edges(3, edges), 2);
  CHECK(od.tree_edges == std::vector<EdgeId>{0, 1});
  CHECK(od.reduced == std::vector<std::uint8_t>{1, 1, 0, 1, 1, 0, 0});
}

TEST_CASE("owner decomposition rejects invalid input") {
  using W = std::uint64_t;
  const auto ok = CsrGraph<W>::from_edges(2, {{0, 1, 1}});
  CHECK_THROWS_AS(owner_decomposition(ok, 0), std::invalid_argument);
  CHECK_THROWS_AS(owner_decomposition(CsrGraph<W>::from_edges(0, {}), 1), std::invalid_argument);
  CHECK_THROWS_AS(owner_decomposition(CsrGraph<W>::from_edges(3, {{0, 1, 1}}), 1),
                  std::invalid_argument);
  CHECK_THROWS_AS(owner_decomposition(CsrGraph<W>::from_edges(2, {{0, 1, 0}}), 1),
                  std::invalid_argument);
  CHECK_THROWS_AS(owner_decomposition(CsrGraph<double>::from_edges(2, {{0, 1, 0.0}}), 1),
                  std::invalid_argument);
}
