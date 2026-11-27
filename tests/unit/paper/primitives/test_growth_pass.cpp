// Tests for GrowthPass (Sec 3.2, Algorithm 1, Lemma 3.4) and the Kruskal reference.
#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <vector>

#include "paper/primitives/mst.hpp"
#include "sssp/graph.hpp"
#include "unit/paper/primitives/mst_test_support.hpp"

using sssp::Edge;
using sssp::EdgeId;
using sssp::Vertex;
using sssp::paper::growth_pass;
using sssp::paper::kUnboundedGrowth;
using sssp::paper::mst_kruskal_ref;
using sssp::paper::MstEdgeKey;
using sssp::paper::QuotientGraph;
using sssp::paper::singleton_quotient;

namespace {

template <class W>
QuotientGraph<W> quotient_of(std::size_t n, const std::vector<Edge<W>>& edges) {
  return singleton_quotient<W>(n, edges);
}

template <class W>
W w(std::uint64_t r) {
  return mst_test::make_weight<W>(r);
}

}  // namespace

TEMPLATE_TEST_CASE("MstEdgeKey orders by length, then input index", "[paper][mst]", std::uint64_t,
                   double) {
  using K = MstEdgeKey<TestType>;
  REQUIRE(K{w<TestType>(1), 9} < K{w<TestType>(2), 0});
  REQUIRE(K{w<TestType>(2), 3} < K{w<TestType>(2), 4});
  REQUIRE_FALSE(K{w<TestType>(2), 4} < K{w<TestType>(2), 3});
  REQUIRE_FALSE(K{w<TestType>(2), 4} < K{w<TestType>(2), 4});
}

TEMPLATE_TEST_CASE("singleton_quotient drops loops and keeps original indices", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  const std::vector<Edge<T>> edges{{0, 1, w<T>(3)}, {2, 2, w<T>(1)}, {2, 1, w<T>(4)}};
  const auto q = quotient_of<T>(3, edges);
  REQUIRE(q.num_vertices == 3);
  REQUIRE(q.edges.size() == 2);
  REQUIRE(q.edges[0].id == 0);
  REQUIRE(q.edges[1].id == 2);
  REQUIRE(q.edges[1].a == 2);
  REQUIRE(q.edges[1].b == 1);
  REQUIRE(q.edges[1].weight == w<T>(4));
  const std::vector<Edge<T>> bad{{0, 3, w<T>(1)}};
  REQUIRE_THROWS_AS(quotient_of<T>(3, bad), std::out_of_range);
}

TEMPLATE_TEST_CASE("growth_pass with B = infinity grows one tree over a connected path",
                   "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  // 0 -1- 1 -5- 2 -2- 3, given out of order.
  const std::vector<Edge<T>> edges{{2, 3, w<T>(2)}, {0, 1, w<T>(1)}, {1, 2, w<T>(5)}};
  const auto r = growth_pass(quotient_of<T>(4, edges), kUnboundedGrowth);
  REQUIRE(r.selected == std::vector<EdgeId>{1, 2, 0});  // selection order of the growth from 0
  REQUIRE(r.quotient.num_vertices == 1);
  REQUIRE(r.quotient.edges.empty());
  REQUIRE(r.component_of == std::vector<Vertex>{0, 0, 0, 0});
  REQUIRE(r.stats.growths == 1);
  REQUIRE(r.stats.stopped_empty == 1);
  REQUIRE(r.stats.stopped_large == 0);
  REQUIRE(r.stats.stopped_attached == 0);
  REQUIRE(r.stats.extractions == 3);
  REQUIRE(r.stats.vertices_in == 4);
  REQUIRE(r.stats.edges_in == 3);
  REQUIRE(r.stats.vertices_out == 1);
  REQUIRE(r.stats.edges_out == 0);
  REQUIRE(r.stats.limit == kUnboundedGrowth);
}

TEMPLATE_TEST_CASE("growth_pass stops a growth whose heap exceeds B (Alg 1 l.9)", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  // Star centred at 0 with 4 leaves: the growth from 0 sees 4 > 2 outside vertices and stops
  // without extracting; each leaf then extracts its edge to 0 and attaches (Alg 1 l.12-13).
  const std::vector<Edge<T>> edges{
      {0, 1, w<T>(1)}, {0, 2, w<T>(2)}, {0, 3, w<T>(3)}, {0, 4, w<T>(4)}};
  const auto r = growth_pass(quotient_of<T>(5, edges), 2);
  REQUIRE(r.selected == std::vector<EdgeId>{0, 1, 2, 3});
  REQUIRE(r.stats.growths == 5);
  REQUIRE(r.stats.stopped_large == 1);
  REQUIRE(r.stats.stopped_attached == 4);
  REQUIRE(r.stats.stopped_empty == 0);
  REQUIRE(r.stats.extractions == 4);
  REQUIRE(r.stats.max_extract_heap_size == 1);
  REQUIRE(r.quotient.num_vertices == 1);
  REQUIRE(r.stats.limit == 2);
}

TEMPLATE_TEST_CASE("growth_pass leaves the quotient by the components of F", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  // Two stars {0;1,2,3} and {4;5,6,7} joined by the heavy edge 3-7 (id 6). With B = 2 both
  // centres stop on a large heap, every leaf attaches to its centre, and the bridge survives as
  // the single quotient edge.
  const std::vector<Edge<T>> edges{{0, 1, w<T>(1)}, {0, 2, w<T>(1)}, {0, 3, w<T>(1)},
                                   {4, 5, w<T>(1)}, {4, 6, w<T>(1)}, {4, 7, w<T>(1)},
                                   {3, 7, w<T>(10)}};
  const auto r = growth_pass(quotient_of<T>(8, edges), 2);
  REQUIRE(std::set<EdgeId>(r.selected.begin(), r.selected.end()) ==
          std::set<EdgeId>{0, 1, 2, 3, 4, 5});
  REQUIRE(r.selected.size() == 6);
  REQUIRE(r.component_of == std::vector<Vertex>{0, 0, 0, 0, 1, 1, 1, 1});
  REQUIRE(r.quotient.num_vertices == 2);
  REQUIRE(r.quotient.edges.size() == 1);
  REQUIRE(r.quotient.edges[0].id == 6);
  REQUIRE(r.quotient.edges[0].weight == w<T>(10));
  REQUIRE(std::min(r.quotient.edges[0].a, r.quotient.edges[0].b) == 0);
  REQUIRE(std::max(r.quotient.edges[0].a, r.quotient.edges[0].b) == 1);
  REQUIRE(r.stats.growths == 8);
  REQUIRE(r.stats.stopped_large == 2);
  REQUIRE(r.stats.stopped_attached == 6);
  // Lemma 3.4: at most 1 + 2 l'/B quotient vertices for a connected input.
  REQUIRE(r.stats.vertices_out * 2 <= 2 + 2 * r.stats.edges_in);
}

TEMPLATE_TEST_CASE("growth_pass keeps the lightest parallel edge by (length, index)",
                   "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  const std::vector<Edge<T>> edges{
      {0, 1, w<T>(5)}, {1, 0, w<T>(3)}, {0, 1, w<T>(3)}, {0, 0, w<T>(0)}};
  const auto r = growth_pass(quotient_of<T>(2, edges), 1);
  REQUIRE(r.selected == std::vector<EdgeId>{1});  // weight 3, smaller index than edge 2
  REQUIRE(r.quotient.num_vertices == 1);
  REQUIRE(r.quotient.edges.empty());
}

TEMPLATE_TEST_CASE("growth_pass heap size counts distinct outside vertices (A11)", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  // Vertex 0 has 4 boundary edges but only 2 distinct outside vertices: with B = 2 the growth
  // must continue (2 is not > 2) and cover all three vertices.
  const std::vector<Edge<T>> edges{
      {0, 1, w<T>(1)}, {0, 1, w<T>(2)}, {1, 0, w<T>(3)}, {0, 2, w<T>(4)}};
  const auto r = growth_pass(quotient_of<T>(3, edges), 2);
  REQUIRE(r.selected == std::vector<EdgeId>{0, 3});
  REQUIRE(r.stats.growths == 1);
  REQUIRE(r.stats.stopped_empty == 1);
  REQUIRE(r.stats.max_extract_heap_size == 2);
  REQUIRE(r.quotient.num_vertices == 1);
}

TEMPLATE_TEST_CASE("growth_pass handles a single vertex and an edgeless quotient", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  const std::vector<Edge<T>> loop_only{{0, 0, w<T>(7)}};
  const auto one = growth_pass(quotient_of<T>(1, loop_only), 1);
  REQUIRE(one.selected.empty());
  REQUIRE(one.quotient.num_vertices == 1);
  REQUIRE(one.stats.growths == 1);
  const auto three = growth_pass(quotient_of<T>(3, std::vector<Edge<T>>{}), kUnboundedGrowth);
  REQUIRE(three.selected.empty());
  REQUIRE(three.quotient.num_vertices == 3);
  REQUIRE(three.component_of == std::vector<Vertex>{0, 1, 2});
  REQUIRE(three.stats.stopped_empty == 3);
  const auto none = growth_pass(quotient_of<T>(0, std::vector<Edge<T>>{}), 1);
  REQUIRE(none.selected.empty());
  REQUIRE(none.quotient.num_vertices == 0);
}

TEMPLATE_TEST_CASE("growth_pass rejects B = 0 and malformed quotients", "[paper][mst]",
                   std::uint64_t, double) {
  using T = TestType;
  const auto q = quotient_of<T>(2, std::vector<Edge<T>>{{0, 1, w<T>(1)}});
  REQUIRE_THROWS_AS(growth_pass(q, 0), std::invalid_argument);
  QuotientGraph<T> bad;
  bad.num_vertices = 2;
  bad.edges.push_back({0, 2, w<T>(1), 0});
  REQUIRE_THROWS_AS(growth_pass(bad, 1), std::out_of_range);
}

TEMPLATE_TEST_CASE("mst_kruskal_ref returns the minimum spanning forest under (length, index)",
                   "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  // Triangle with a tie at weight 2 (edge 1 beats edge 2 by index), an isolated vertex 3, a
  // separate edge 4-5, and a loop.
  const std::vector<Edge<T>> edges{
      {0, 1, w<T>(1)}, {1, 2, w<T>(2)}, {0, 2, w<T>(2)}, {4, 5, w<T>(9)}, {5, 5, w<T>(0)}};
  REQUIRE(mst_kruskal_ref<T>(6, edges) == std::vector<EdgeId>{0, 1, 3});
  REQUIRE(mst_kruskal_ref<T>(1, std::vector<Edge<T>>{}).empty());
}

TEMPLATE_TEST_CASE("growth_pass selects only MST edges, each once, within the Lemma 3.4 bounds",
                   "[paper][mst]", std::uint64_t, double) {
  using T = TestType;
  for (std::uint64_t seed = 1; seed <= 40; ++seed) {
    const std::size_t n = 1 + seed % 23;
    const bool connected = seed % 4 != 0;
    const auto edges = mst_test::random_edges<T>(n, 2 * n, seed, seed % 3 == 0 ? 2 : 50, connected);
    const auto ref = mst_kruskal_ref<T>(n, edges);
    const std::set<EdgeId> mst(ref.begin(), ref.end());
    const std::size_t comps = mst_test::count_components(n, edges);
    for (const std::size_t limit :
         {std::size_t{1}, std::size_t{2}, std::size_t{3}, kUnboundedGrowth}) {
      const auto q = quotient_of<T>(n, edges);
      const auto r = growth_pass(q, limit);
      const std::set<EdgeId> sel(r.selected.begin(), r.selected.end());
      REQUIRE(sel.size() == r.selected.size());  // no edge selected twice
      for (const EdgeId id : sel) REQUIRE(mst.count(id) == 1);
      REQUIRE(r.quotient.num_vertices == n - r.selected.size());
      REQUIRE(r.stats.vertices_out == r.quotient.num_vertices);
      // Lemma 3.4 (per connected component of the input): at most c + 2 l'/B vertices.
      if (limit != kUnboundedGrowth) {
        REQUIRE(r.stats.vertices_out * limit <= comps * limit + 2 * r.stats.edges_in);
        REQUIRE(r.stats.max_extract_heap_size <= std::min(n, limit));
      } else {
        REQUIRE(r.stats.vertices_out == comps);  // B = +inf covers each component
      }
      REQUIRE(r.stats.extractions == r.selected.size());
      // The output quotient holds exactly the input edges joining different components.
      REQUIRE(r.component_of.size() == n);
      std::size_t crossing = 0;
      for (const auto& e : q.edges) {
        if (r.component_of[e.a] != r.component_of[e.b]) ++crossing;
      }
      REQUIRE(r.quotient.edges.size() == crossing);
      for (const auto& e : r.quotient.edges) {
        REQUIRE(e.a != e.b);
        REQUIRE(e.a < r.quotient.num_vertices);
        REQUIRE(e.b < r.quotient.num_vertices);
        const auto& orig = edges[e.id];
        REQUIRE(e.weight == orig.weight);
        REQUIRE(((r.component_of[orig.u] == e.a && r.component_of[orig.v] == e.b) ||
                 (r.component_of[orig.u] == e.b && r.component_of[orig.v] == e.a)));
      }
      // Every component of F is one output vertex: F's endpoints agree on component_of.
      for (const EdgeId id : r.selected) {
        REQUIRE(r.component_of[edges[id].u] == r.component_of[edges[id].v]);
      }
    }
  }
}
