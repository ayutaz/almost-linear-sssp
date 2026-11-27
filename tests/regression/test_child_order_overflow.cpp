// Regression 20261127-child-order-overflow (record: cases/20261127-child-order-overflow.json).
//
// Valid uint64 input: the path 0 - 1 - 2 - 3 with three edges of length q = ceil(2^64 / 5).
// Every distance fits (at most 3q < 2^64 - 1; exact reference
// tests/python/reference/exact_sssp.py). Algorithm 4 (Mode A, base case) returns one node of
// height q with three links of length q; its mass h + 3q = 4q fits as well. Sec 4.2
// (order_children) then threw std::overflow_error("integer length overflow").
//
// Cause: the forest distance of Lemma 4.5 was formed as the signed pair
// (rd(a) + rd(b), 2 rd(lca)) of root distances of the link forest F, and the pair of a segment
// was accumulated before one subtraction. Moving from state 2 to state 3 adds
// rd(2) + rd(3) = 2q + 3q = 5q >= 2^64 to the positive component, although the distance itself
// is q. Root distances and these sums are auxiliary quantities, not path lengths: by the
// overflow rule of include/sssp/paper.hpp they must not raise.
//
// Expected: Evaluate returns least emissions equal to the exact distances, for every source.
//
// Found by: P3-EVAL-ARITH-FIX (search for further overflows after saturating the Sec 4.1 mass):
// the path of four edges of length 2^62 from vertex 2, and lt::clustered(209953, 107) with
// bridges of 2^52 or more. Those two also overflow the mass S(C), which Evaluate rejects until
// D-23 (P3-EVAL-ARITH-FIX) is merged; this case is minimised so that only Sec 4.2 overflows:
// four vertices are needed (with three, every sum of Sec 4.2 is at most the mass), and
// 5q >= 2^64 gives the least weight. The four-edge path is covered at the order_children level
// below, on the normalised hierarchy that normalize() returns for it.

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "sssp/graph.hpp"

namespace {

using W = std::uint64_t;
using sssp::paper::NodeId;
using sssp::paper::StateId;

constexpr W kInf = sssp::WeightTraits<W>::infinity();

std::vector<W> least_emissions(const sssp::CsrGraph<W>& g, sssp::Vertex src) {
  const auto lh = sssp::paper::local_hierarchy(g);
  sssp::paper::SeededOrderedHierarchy<W> s{
      lh.hierarchy, {{src, W{0}}}, {{lh.native_entry[src], W{0}}}};
  const auto valid = sssp::paper::validate_seeded_ordered_hierarchy(s, g);
  INFO(valid.summary());
  REQUIRE(valid.ok());
  const auto r = sssp::paper::evaluate(s);
  std::vector<W> best(g.num_vertices(), kInf);
  for (const auto& [x, v] : r.emissions) {
    const auto p = lh.hierarchy.projection[x];
    if (v < best[p]) best[p] = v;
  }
  return best;
}

sssp::CsrGraph<W> path(std::size_t n, W len) {
  std::vector<sssp::Edge<W>> edges;
  for (sssp::Vertex v = 0; v + 1 < n; ++v) edges.push_back({v, v + 1, len});
  return sssp::CsrGraph<W>::from_edges(n, edges);
}

}  // namespace

TEST_CASE(
    "20261127-child-order-overflow: Alg 4 + Evaluate on a path of three edges of "
    "ceil(2^64 / 5)",
    "[regression][paper]") {
  constexpr W q = 3689348814741910324ULL;  // 5 q = 2^64 + 4
  const auto g = path(4, q);
  // Exact reference (exact_sssp.py): the distances on a path are multiples of q.
  std::vector<W> got;
  REQUIRE_NOTHROW(got = least_emissions(g, 0));
  CHECK(got == std::vector<W>{0, q, 2 * q, 3 * q});
  REQUIRE_NOTHROW(got = least_emissions(g, 1));
  CHECK(got == std::vector<W>{q, 0, q, 2 * q});
  REQUIRE_NOTHROW(got = least_emissions(g, 2));
  CHECK(got == std::vector<W>{2 * q, q, 0, q});
  REQUIRE_NOTHROW(got = least_emissions(g, 3));
  CHECK(got == std::vector<W>{3 * q, 2 * q, q, 0});
}

TEST_CASE(
    "20261127-child-order-overflow: the same path with shorter edges (passes before the "
    "fix)",
    "[regression][paper]") {
  constexpr W q = W{1} << 61;  // 5 q < 2^64
  CHECK(least_emissions(path(4, q), 0) == std::vector<W>{0, q, 2 * q, 3 * q});
}

TEST_CASE("20261127-child-order-overflow: order_children on the path of four edges of 2^62",
          "[regression][paper]") {
  // The reported case. Its mass 5 * 2^62 does not fit, so normalize() is run on the same path
  // with unit edges and every length of the result is then multiplied by 2^62: Sec 4.1 only
  // compares and adds lengths, so its decisions are the same at both scales, and order_children
  // does not read the mass (set to the D-23 saturation value here).
  constexpr W q = W{1} << 62;
  const auto lh = sssp::paper::local_hierarchy(path(5, W{1}));
  auto nh = sssp::paper::normalize(lh.hierarchy);
  for (auto& a : nh.forest.graph.arcs) a.length *= q;
  for (auto& nd : nh.forest.nodes) {
    nd.height *= q;
    for (auto& lk : nd.links) lk.length *= q;
  }
  for (NodeId v = 0; v < nh.num_nodes(); ++v) {
    nh.unit[v] *= q;
    nh.original_height[v] *= q;
    if (!nh.is_leaf(v)) nh.mass[v] = kInf;
  }
  REQUIRE(nh.forest.roots.size() == 1);
  const NodeId x = nh.forest.roots[0];
  REQUIRE(nh.forest.nodes[x].children.size() == 5);
  REQUIRE(nh.forest.nodes[x].links.size() == 4);
  REQUIRE(nh.unit[x] == q);

  // The link forest is the path itself: the root distance of its far end is 4 * 2^62 = 2^64.
  sssp::paper::LinkForest<W> lf;
  REQUIRE_NOTHROW(lf = sssp::paper::build_link_forest(nh));
  __extension__ typedef unsigned __int128 U128;
  U128 farthest = 0;
  for (StateId s = 0; s < 5; ++s) {
    if (lf.root_distance[s] > farthest) farthest = lf.root_distance[s];
  }
  CHECK(farthest == U128{4} * q);

  // No light arcs (every arc has length nu(X)): the five leaves are DFS roots in F preorder and
  // each gap is the F distance between consecutive leaves of the path, 2^62 (l.747).
  REQUIRE_NOTHROW(sssp::paper::order_children(nh));
  const auto& ch = nh.forest.nodes[x].children;
  for (std::size_t k = 0; k + 1 < ch.size(); ++k) CHECK(nh.gap_next[ch[k]] == q);
  CHECK(nh.gap_next[ch.back()] == W{0});
}
