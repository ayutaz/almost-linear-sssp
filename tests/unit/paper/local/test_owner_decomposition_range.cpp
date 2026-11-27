// Range of the arithmetic of the owner decomposition (Sec 5.1) at the top of the weight type.
//
// Overflow rule: std::overflow_error is reserved for a tentative path length from the source
// that reaches the infinity sentinel. Sec 5.1 has no source; its sums (root distances of the
// MST for the reduced-graph test of l.965-966, R_v = D_T(v, p(v)) of l.986-988) are auxiliary
// quantities and must not raise. Specified here:
//   * integer W: root distances are exact (TreeLength<W> is a 128-bit integer), the test
//     w <= D_T(u, v) is decided exactly, and an R_v that is not below the sentinel is infinity();
//   * double: nothing raises; a root distance or R_v beyond the range is +infinity;
//   * scale invariance: multiplying every weight by 2^s changes no decision of Sec 5.1 (the
//     construction only compares and adds lengths, Sec 2), also when the sums leave 64 bits.
// The component is tested directly: some of the small graphs below (e.g. two edges of 2^63) have
// no source from which every tentative length fits; Sec 5.1 must not raise on them either.
// Expected values are computed by hand or by brute-force tree walks in 128-bit integers written
// in this file; none is copied from a run of the implementation.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/owner/aux_metric.hpp"
#include "sssp/graph.hpp"
#include "unit/paper/local/local_test_support.hpp"

namespace {

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeId;
using sssp::Vertex;
using sssp::paper::owner_decomposition;
using sssp::paper::OwnerDecomposition;
namespace lt = sssp::test::local;

using U64 = std::uint64_t;
__extension__ typedef unsigned __int128 U128;

constexpr U64 kInf = sssp::WeightTraits<U64>::infinity();
constexpr U64 kQ = U64{1} << 62;

std::vector<Edge<U64>> uniform_path(std::size_t n, U64 len) {
  std::vector<Edge<U64>> edges;
  for (Vertex v = 0; v + 1 < n; ++v) edges.push_back({v, v + 1, len});
  return edges;
}

/// D_T(u, v) by walking the tree path, in 128 bits (fewer than 2^32 edges below 2^64 each).
U128 tree_distance(const std::vector<Edge<U64>>& edges, const lt::RootedTree<U64>& t, Vertex u,
                   Vertex v) {
  U128 sum = 0;
  for (const EdgeId id : lt::tree_path(t, u, v)) sum += edges[id].weight;
  return sum;
}

}  // namespace

TEST_CASE("the root distance of the owner decomposition is a 128-bit integer for uint64") {
  STATIC_REQUIRE(std::is_same_v<sssp::paper::TreeLength<U64>, U128>);
  STATIC_REQUIRE(std::is_same_v<sssp::paper::TreeLength<double>, double>);
  STATIC_REQUIRE(
      std::is_same_v<decltype(OwnerDecomposition<U64>{}.root_distance), std::vector<U128>>);
}

TEST_CASE("owner decomposition: root distance equal to the sentinel (2^63 + 2^63 - 1)") {
  const U64 a = U64{1} << 63;
  const U64 c = a - 1;
  const auto g = CsrGraph<U64>::from_edges(3, {{0, 1, a}, {1, 2, c}});
  OwnerDecomposition<U64> od;
  REQUIRE_NOTHROW(od = owner_decomposition(g, 2));
  CHECK(od.tree_edges == std::vector<EdgeId>{0, 1});
  REQUIRE(od.root_distance.size() == 3);
  CHECK(od.root_distance[0] == U128{0});
  CHECK(od.root_distance[1] == U128{a});
  CHECK(od.root_distance[2] == U128{kInf});  // exact: a + c = 2^64 - 1
  CHECK(od.marked == std::vector<std::uint8_t>{0, 1, 0});
  CHECK(od.pivots == std::vector<Vertex>{1});
  CHECK(od.corridors.empty());
  CHECK(od.radius == std::vector<U64>{a, 0, c});
  CHECK(od.owner_distance == std::vector<U64>{a, 0, c});
  CHECK(od.owner_hops == std::vector<std::uint32_t>{1, 0, 1});
}

TEST_CASE("owner decomposition of the path of four edges of 2^62 with b = 2") {
  const auto g = CsrGraph<U64>::from_edges(5, uniform_path(5, kQ));
  OwnerDecomposition<U64> od;
  REQUIRE_NOTHROW(od = owner_decomposition(g, 2));
  // Root distances 0, q, 2q, 3q, 4q = 2^64.
  REQUIRE(od.root_distance.size() == 5);
  for (Vertex v = 0; v < 5; ++v) CHECK(od.root_distance[v] == U128{v} * kQ);
  CHECK(od.root_distance[4] == U128{1} << 64);
  // Postorder from root 0: 4 (1), 3 (2: mark), 2 (1), 1 (2: mark), 0 (1).
  CHECK(od.marked == std::vector<std::uint8_t>{0, 1, 0, 1, 0});
  CHECK(od.pivots == std::vector<Vertex>{1, 3});
  // Corridor 1 - 2 - 3: both edges have length q, the smaller id (edge 1) is cut (Q-31).
  REQUIRE(od.corridors.size() == 1);
  CHECK(od.corridors[0].edge == 1);
  CHECK(od.corridors[0].weight == kQ);
  CHECK(od.cut == std::vector<std::uint8_t>{0, 1, 0, 0});
  CHECK(od.owner == std::vector<std::uint32_t>{0, 0, 1, 1, 1});
  CHECK(od.radius == std::vector<U64>{kQ, 0, kQ, 0, kQ});
  CHECK(od.owner_distance == std::vector<U64>{kQ, 0, kQ, 0, kQ});
  CHECK(od.owner_hops == std::vector<std::uint32_t>{1, 0, 1, 0, 1});
  CHECK(od.reduced == std::vector<std::uint8_t>{1, 1, 1, 1});
  CHECK(od.sigma(0, 1) == kQ);
}

TEST_CASE("owner decomposition: R_v that does not fit is the infinity sentinel") {
  // Path of four edges of 2^62 with b = 6 > N: no mark, the root 0 is the single pivot, so
  // R_v = D_T(0, v) = v * 2^62 and R_4 = 2^64 does not fit. r_v = b_T(0, v) = 2^62.
  const auto g = CsrGraph<U64>::from_edges(5, uniform_path(5, kQ));
  OwnerDecomposition<U64> od;
  REQUIRE_NOTHROW(od = owner_decomposition(g, 6));
  CHECK(od.pivots == std::vector<Vertex>{0});
  CHECK(od.radius == std::vector<U64>{0, kQ, kQ, kQ, kQ});
  CHECK(od.owner_distance == std::vector<U64>{0, kQ, 2 * kQ, 3 * kQ, kInf});
  CHECK(od.owner_hops == std::vector<std::uint32_t>{0, 1, 2, 3, 4});

  // 2^64 - 2 is the largest R_v that is returned exactly: 2^63 + (2^63 - 2).
  const U64 a = U64{1} << 63;
  const auto fits =
      owner_decomposition(CsrGraph<U64>::from_edges(3, {{0, 1, a}, {1, 2, a - 2}}), 4);
  CHECK(fits.pivots == std::vector<Vertex>{0});
  CHECK(fits.owner_distance == std::vector<U64>{0, a, kInf - 1});
  // One more unit reaches the sentinel.
  const auto sat = owner_decomposition(CsrGraph<U64>::from_edges(3, {{0, 1, a}, {1, 2, a - 1}}), 4);
  CHECK(sat.owner_distance == std::vector<U64>{0, a, kInf});
  // Saturation is absorbing along the cell: 2^63, 2^63, 1.
  const auto past =
      owner_decomposition(CsrGraph<U64>::from_edges(4, {{0, 1, a}, {1, 2, a}, {2, 3, U64{1}}}), 5);
  CHECK(past.owner_distance == std::vector<U64>{0, a, kInf, kInf});
  CHECK(past.radius == std::vector<U64>{0, a, a, a});
}

TEST_CASE("reduced graph: w <= D_T(u, v) is exact when the root distances exceed 64 bits") {
  // T: stem 0 - 1 - 2 and two legs 2 - 3 - 4 and 2 - 5 - 6, every edge q = 2^62 (edges 0..5).
  // Root distances: 0, q, 2q, 3q, 4q = 2^64, 3q, 4q = 2^64.
  //   edge 6 = (3, 6): D_T = 3q (lca 2), w = 3q        -> kept (equality)
  //   edge 7 = (4, 5): D_T = 3q (lca 2), w = 3q + 1    -> dropped
  //   edge 8 = (4, 6): D_T = 4q = 2^64,  w = 2^64 - 2  -> kept (every finite w is below D_T)
  //   edge 9 = (1, 4): D_T = 3q (lca 1), w = 3q + 1    -> dropped
  //   edge 10 = (0, 6): D_T = 4q (lca 0), w = 3q + 1   -> kept
  const std::vector<Edge<U64>> edges{{0, 1, kQ},         {1, 2, kQ},         {2, 3, kQ},
                                     {3, 4, kQ},         {2, 5, kQ},         {5, 6, kQ},
                                     {3, 6, 3 * kQ},     {4, 5, 3 * kQ + 1}, {4, 6, kInf - 1},
                                     {1, 4, 3 * kQ + 1}, {0, 6, 3 * kQ + 1}};
  const auto g = CsrGraph<U64>::from_edges(7, edges);
  for (const std::size_t b : {1U, 2U, 3U, 8U}) {
    INFO("b = " << b);
    OwnerDecomposition<U64> od;
    REQUIRE_NOTHROW(od = owner_decomposition(g, b));
    CHECK(od.tree_edges == std::vector<EdgeId>{0, 1, 2, 3, 4, 5});
    CHECK(od.reduced == std::vector<std::uint8_t>{1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1});
    REQUIRE(od.root_distance.size() == 7);
    CHECK(od.root_distance[4] == U128{1} << 64);
    CHECK(od.root_distance[5] == U128{3} * kQ);
    CHECK(od.root_distance[6] == U128{1} << 64);
  }
}

TEST_CASE("owner decomposition is invariant under scaling the weights by 2^53") {
  // Base weights are at most 1000 < 2^10, so the scaled weights are below 2^63 and valid, while
  // tree distances of the scaled graph leave 64 bits on all but the smallest graphs.
  constexpr unsigned kShift = 53;
  std::size_t wide_root_distances = 0;
  std::size_t saturated = 0;
  std::size_t dropped = 0;
  std::uint64_t seed = 9000;
  for (const auto family : {lt::Family::kPath, lt::Family::kSpider, lt::Family::kRandomTree,
                            lt::Family::kSparse, lt::Family::kGrid, lt::Family::kComplete}) {
    for (const auto weights : {lt::Weights::kTies, lt::Weights::kUniform}) {
      for (const std::size_t n : {2U, 5U, 9U, 17U, 33U}) {
        ++seed;
        const auto base = lt::make_graph<U64>(family, weights, n, seed);
        auto scaled = base;
        for (auto& e : scaled) e.weight <<= kShift;
        const auto gb = CsrGraph<U64>::from_edges(n, base);
        const auto gs = CsrGraph<U64>::from_edges(n, scaled);
        for (const std::size_t b : std::set<std::size_t>{1, 2, 3, 4, n, n + 1}) {
          INFO(lt::name(family) << " " << lt::name(weights) << " n = " << n << " seed = " << seed
                                << " b = " << b);
          const auto ob = owner_decomposition(gb, b);
          OwnerDecomposition<U64> os;
          REQUIRE_NOTHROW(os = owner_decomposition(gs, b));

          // Decisions are the same at both scales.
          CHECK(os.tree_edges == ob.tree_edges);
          CHECK(os.tree_parent == ob.tree_parent);
          CHECK(os.tree_parent_edge == ob.tree_parent_edge);
          CHECK(os.reduced == ob.reduced);
          CHECK(os.marked == ob.marked);
          CHECK(os.pivots == ob.pivots);
          CHECK(os.pivot_index == ob.pivot_index);
          CHECK(os.cut == ob.cut);
          CHECK(os.owner == ob.owner);
          CHECK(os.owner_hops == ob.owner_hops);
          CHECK(os.merge_parent == ob.merge_parent);
          CHECK(os.merge_corridor == ob.merge_corridor);
          REQUIRE(os.corridors.size() == ob.corridors.size());
          for (std::size_t i = 0; i < ob.corridors.size(); ++i) {
            CHECK(os.corridors[i].lower == ob.corridors[i].lower);
            CHECK(os.corridors[i].upper == ob.corridors[i].upper);
            CHECK(os.corridors[i].edge == ob.corridors[i].edge);
            CHECK(os.corridors[i].weight == ob.corridors[i].weight << kShift);
          }

          // Lengths scale; the root distance exactly, R_v up to the sentinel.
          const auto tree = lt::root_tree(n, scaled, os.tree_edges, os.root);
          REQUIRE(os.root_distance.size() == n);
          REQUIRE(ob.root_distance.size() == n);
          for (Vertex v = 0; v < n; ++v) {
            CHECK(os.radius[v] == ob.radius[v] << kShift);
            const U128 rd = U128{ob.root_distance[v]} << kShift;
            CHECK(os.root_distance[v] == rd);
            CHECK(os.root_distance[v] == tree_distance(scaled, tree, os.root, v));
            if (rd >= U128{kInf}) ++wide_root_distances;
            const U128 big_r = U128{ob.owner_distance[v]} << kShift;
            CHECK(big_r == tree_distance(scaled, tree, v, os.pivots[os.owner[v]]));
            if (big_r >= U128{kInf}) {
              ++saturated;
              CHECK(os.owner_distance[v] == kInf);
            } else {
              CHECK(U128{os.owner_distance[v]} == big_r);
            }
          }

          // The reduced graph against brute-force tree walks (l.965-966).
          std::vector<std::uint8_t> in_tree(scaled.size(), 0);
          for (const EdgeId id : os.tree_edges) in_tree[id] = 1;
          for (std::size_t i = 0; i < scaled.size(); ++i) {
            const auto& e = scaled[i];
            const bool want =
                in_tree[i] != 0 || U128{e.weight} <= tree_distance(scaled, tree, e.u, e.v);
            CHECK((os.reduced[i] != 0) == want);
            if (!want) ++dropped;
          }
        }
      }
    }
  }
  // The property was exercised beyond 64 bits.
  CHECK(wide_root_distances > 100);
  CHECK(saturated > 100);
  CHECK(dropped > 100);
}

TEST_CASE("owner decomposition of double weights does not raise when a tree sum overflows") {
  constexpr double kHuge = 1e308;
  constexpr double kDInf = std::numeric_limits<double>::infinity();
  const auto g = CsrGraph<double>::from_edges(3, {{0, 1, kHuge}, {1, 2, kHuge}});
  // From vertex 1 every distance is 1e308; the root distance of vertex 2 is 2e308 > DBL_MAX.
  OwnerDecomposition<double> od;
  REQUIRE_NOTHROW(od = owner_decomposition(g, 2));
  CHECK(od.pivots == std::vector<Vertex>{1});
  CHECK(od.root_distance == std::vector<double>{0.0, kHuge, kDInf});
  CHECK(od.reduced == std::vector<std::uint8_t>{1, 1});
  CHECK(od.radius == std::vector<double>{kHuge, 0.0, kHuge});
  CHECK(od.owner_distance == std::vector<double>{kHuge, 0.0, kHuge});

  // b = 3: the root is the single pivot and R_2 = 2e308 is +infinity.
  REQUIRE_NOTHROW(od = owner_decomposition(g, 3));
  CHECK(od.pivots == std::vector<Vertex>{0});
  CHECK(od.radius == std::vector<double>{0.0, kHuge, kHuge});
  CHECK(od.owner_distance == std::vector<double>{0.0, kHuge, kDInf});
}

TEST_CASE("reduced graph of double weights near the top of the range") {
  // Path 0 - 1 - 2 - 3 with edges of 0.5e308: root distances 0.5e308, 1e308, 1.5e308.
  //   edge 3 = (1, 2): D_T = 0.5e308, w = 0.9e308 -> dropped; w + 2 rd(1) = 1.9e308 overflows
  //                    while rd(1) + rd(2) = 1.5e308 does not
  //   edge 4 = (0, 3): D_T = 1.5e308, w = 1.25e308 -> kept
  // (When both sides of the test overflow the edge is kept, see owner_decomposition.hpp; that
  // conservative outcome is not a statement of the paper and is not asserted here.)
  constexpr double kA = 0.5e308;
  const std::vector<Edge<double>> edges{
      {0, 1, kA}, {1, 2, kA}, {2, 3, kA}, {1, 2, 0.9e308}, {0, 3, 1.25e308}};
  const auto g = CsrGraph<double>::from_edges(4, edges);
  OwnerDecomposition<double> od;
  REQUIRE_NOTHROW(od = owner_decomposition(g, 2));
  CHECK(od.tree_edges == std::vector<EdgeId>{0, 1, 2});
  CHECK(od.root_distance == std::vector<double>{0.0, kA, 1e308, 1.5e308});
  CHECK(od.reduced == std::vector<std::uint8_t>{1, 1, 1, 0, 1});
}

TEST_CASE("a saturated R_v never enters a projected length silently (Sec 7.1)") {
  // Path 0 - 1 - ... - 9 with nine edges of 2^62 and b = 5. Postorder from root 0: vertex 5 and
  // vertex 0 are marked; the corridor 0 - ... - 5 has equal edges, so edge 0 is cut (Q-31). The
  // cell of pivot 5 is {1, ..., 9}: R_1 = R_9 = 4 * 2^62 = 2^64 does not fit.
  const auto g = CsrGraph<U64>::from_edges(10, uniform_path(10, kQ));
  OwnerDecomposition<U64> od;
  REQUIRE_NOTHROW(od = owner_decomposition(g, 5));
  CHECK(od.pivots == std::vector<Vertex>{0, 5});
  CHECK(od.cut == std::vector<std::uint8_t>{1, 0, 0, 0, 0, 0, 0, 0, 0});
  CHECK(od.owner_distance ==
        std::vector<U64>{0, kInf, 3 * kQ, 2 * kQ, kQ, 0, kQ, 2 * kQ, 3 * kQ, kInf});
  // Edge 0 joins the two cells: A_e = R_0 + w + R_1 is not representable, and the consumer
  // that needs it as a length raises instead of using the sentinel as a number.
  CHECK_THROWS_AS(sssp::paper::owner::projected_interowner_edges(g, od), std::overflow_error);
}
