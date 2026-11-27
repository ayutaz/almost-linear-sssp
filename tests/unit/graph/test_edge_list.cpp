// Unit tests for the edge-list utilities (P2-GRAPH-EDGELIST): hand examples whose expected
// values follow from the documented contract in include/sssp/edge_list.hpp.
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "sssp/edge_list.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"

using sssp::Edge;
using sssp::EdgeId;
using sssp::EdgeList;
using sssp::kNoVertex;
using sssp::LoopPolicy;
using sssp::Vertex;

namespace {

template <class W>
void require_edge(const Edge<W>& e, Vertex u, Vertex v, W w) {
  REQUIRE(e.u == u);
  REQUIRE(e.v == v);
  REQUIRE(e.weight == w);
}

/// Parallel edges in both orientations, a tie, two loops on one vertex, an isolated vertex 4.
template <class W>
EdgeList<W> multigraph() {
  return EdgeList<W>{5,
                     {{0, 1, W{5}},
                      {1, 0, W{2}},
                      {1, 2, W{4}},
                      {0, 1, W{2}},
                      {2, 2, W{7}},
                      {2, 2, W{3}},
                      {3, 1, W{1}}}};
}

}  // namespace

// ------------------------------------------------------------------ keep_lightest_per_pair

TEMPLATE_TEST_CASE("keep_lightest_per_pair keeps the (weight, index)-least edge of each pair",
                   "[graph][edge_list]", std::uint64_t, double) {
  using W = TestType;
  const auto g = multigraph<W>();
  const auto r = sssp::keep_lightest_per_pair(g);
  REQUIRE(r.graph.num_vertices == 5);
  // {0,1}: ids 0 (5), 1 (2), 3 (2) -> id 1 (tie with 3 broken by the smaller id), orientation of
  // id 1 kept. {1,2}: id 2. {2,2}: ids 4 (7), 5 (3) -> id 5. {1,3}: id 6.
  REQUIRE(r.original_edge == std::vector<EdgeId>{1, 2, 5, 6});
  REQUIRE(r.graph.edges.size() == 4);
  require_edge<W>(r.graph.edges[0], 1, 0, W{2});
  require_edge<W>(r.graph.edges[1], 1, 2, W{4});
  require_edge<W>(r.graph.edges[2], 2, 2, W{3});
  require_edge<W>(r.graph.edges[3], 3, 1, W{1});
}

TEMPLATE_TEST_CASE("keep_lightest_per_pair with LoopPolicy::kDrop removes every self-loop",
                   "[graph][edge_list]", std::uint64_t, double) {
  using W = TestType;
  const auto r = sssp::keep_lightest_per_pair(multigraph<W>(), LoopPolicy::kDrop);
  REQUIRE(r.graph.num_vertices == 5);
  REQUIRE(r.original_edge == std::vector<EdgeId>{1, 2, 6});
  require_edge<W>(r.graph.edges[0], 1, 0, W{2});
  require_edge<W>(r.graph.edges[1], 1, 2, W{4});
  require_edge<W>(r.graph.edges[2], 3, 1, W{1});
}

TEMPLATE_TEST_CASE("keep_lightest_per_pair leaves a simple graph and empty graphs unchanged",
                   "[graph][edge_list]", std::uint64_t, double) {
  using W = TestType;
  const EdgeList<W> simple{4, {{2, 3, W{9}}, {0, 1, W{1}}, {1, 2, W{1}}}};
  const auto r = sssp::keep_lightest_per_pair(simple);
  REQUIRE(r.original_edge == std::vector<EdgeId>{0, 1, 2});
  require_edge<W>(r.graph.edges[0], 2, 3, W{9});
  require_edge<W>(r.graph.edges[1], 0, 1, W{1});
  require_edge<W>(r.graph.edges[2], 1, 2, W{1});

  const auto empty = sssp::keep_lightest_per_pair(EdgeList<W>{});
  REQUIRE(empty.graph.num_vertices == 0);
  REQUIRE(empty.graph.edges.empty());
  REQUIRE(empty.original_edge.empty());

  const auto no_edges = sssp::keep_lightest_per_pair(EdgeList<W>{3, {}});
  REQUIRE(no_edges.graph.num_vertices == 3);
  REQUIRE(no_edges.graph.edges.empty());
}

TEST_CASE("keep_lightest_per_pair compares uint64 weights exactly near the maximum",
          "[graph][edge_list]") {
  constexpr std::uint64_t kMax = std::numeric_limits<std::uint64_t>::max();
  const EdgeList<std::uint64_t> g{2, {{0, 1, kMax}, {1, 0, kMax - 1}, {0, 1, kMax - 1}}};
  const auto r = sssp::keep_lightest_per_pair(g);
  REQUIRE(r.original_edge == std::vector<EdgeId>{1});
  require_edge<std::uint64_t>(r.graph.edges[0], 1, 0, kMax - 1);
}

TEST_CASE("keep_lightest_per_pair treats -0.0 and +0.0 as equal and keeps the first",
          "[graph][edge_list]") {
  const EdgeList<double> g{2, {{0, 1, -0.0}, {1, 0, 0.0}, {0, 1, 0x1p-1074}}};
  const auto r = sssp::keep_lightest_per_pair(g);
  REQUIRE(r.original_edge == std::vector<EdgeId>{0});
  REQUIRE(r.graph.edges[0].weight == 0.0);
  REQUIRE(std::signbit(r.graph.edges[0].weight));  // the edge is copied, not rebuilt
}

TEST_CASE("keep_lightest_per_pair validates its input", "[graph][edge_list]") {
  REQUIRE_THROWS_AS(sssp::keep_lightest_per_pair(EdgeList<std::uint64_t>{2, {{0, 2, 1}}}),
                    std::out_of_range);
  REQUIRE_THROWS_AS(sssp::keep_lightest_per_pair(EdgeList<double>{2, {{0, 1, -1.0}}}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(sssp::keep_lightest_per_pair(
                        EdgeList<double>{2, {{0, 1, std::numeric_limits<double>::quiet_NaN()}}}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(sssp::keep_lightest_per_pair(
                        EdgeList<double>{2, {{0, 1, std::numeric_limits<double>::infinity()}}}),
                    std::invalid_argument);
}

// ------------------------------------------------------------------ relabel_vertices

TEMPLATE_TEST_CASE("relabel_vertices applies a permutation and keeps edge order and weights",
                   "[graph][edge_list]", std::uint64_t, double) {
  using W = TestType;
  const EdgeList<W> g{3, {{0, 1, W{1}}, {1, 2, W{2}}, {2, 2, W{3}}, {2, 0, W{4}}}};
  const std::vector<Vertex> new_id{2, 0, 1};
  const auto r = sssp::relabel_vertices(g, new_id, 3);
  REQUIRE(r.graph.num_vertices == 3);
  REQUIRE(r.original_edge == std::vector<EdgeId>{0, 1, 2, 3});
  require_edge<W>(r.graph.edges[0], 2, 0, W{1});
  require_edge<W>(r.graph.edges[1], 0, 1, W{2});
  require_edge<W>(r.graph.edges[2], 1, 1, W{3});
  require_edge<W>(r.graph.edges[3], 1, 2, W{4});
}

TEMPLATE_TEST_CASE("relabel_vertices drops kNoVertex vertices with their edges",
                   "[graph][edge_list]", std::uint64_t, double) {
  using W = TestType;
  const EdgeList<W> g{4, {{0, 1, W{1}}, {1, 2, W{2}}, {2, 3, W{3}}, {3, 0, W{4}}, {3, 3, W{5}}}};
  // Keep 0, 2, 3 as 4, 0, 1 in a graph of 6 vertices (new ids 2, 3, 5 stay isolated).
  const std::vector<Vertex> new_id{4, kNoVertex, 0, 1};
  const auto r = sssp::relabel_vertices(g, new_id, 6);
  REQUIRE(r.graph.num_vertices == 6);
  REQUIRE(r.original_edge == std::vector<EdgeId>{2, 3, 4});
  require_edge<W>(r.graph.edges[0], 0, 1, W{3});
  require_edge<W>(r.graph.edges[1], 1, 4, W{4});
  require_edge<W>(r.graph.edges[2], 1, 1, W{5});
}

TEST_CASE("relabel_vertices rejects a bad vertex map", "[graph][edge_list]") {
  const EdgeList<std::uint64_t> g{3, {{0, 1, 1}, {1, 2, 2}}};
  // wrong length
  REQUIRE_THROWS_AS(sssp::relabel_vertices(g, std::vector<Vertex>{0, 1}, 3), std::invalid_argument);
  // id out of range of the new vertex count
  REQUIRE_THROWS_AS(sssp::relabel_vertices(g, std::vector<Vertex>{0, 1, 3}, 3), std::out_of_range);
  // two vertices mapped to the same id (contraction is not a relabelling)
  REQUIRE_THROWS_AS(sssp::relabel_vertices(g, std::vector<Vertex>{0, 1, 0}, 3),
                    std::invalid_argument);
  // too many new vertices for Vertex ids
  REQUIRE_THROWS_AS(
      sssp::relabel_vertices(g, std::vector<Vertex>{0, 1, 2}, static_cast<std::size_t>(kNoVertex)),
      std::length_error);
  // the input edge list itself is validated
  REQUIRE_THROWS_AS(
      sssp::relabel_vertices(EdgeList<std::uint64_t>{2, {{0, 5, 1}}}, std::vector<Vertex>{0, 1}, 2),
      std::out_of_range);
}

// ------------------------------------------------------------------ connected_components

TEMPLATE_TEST_CASE("connected_components numbers components by their smallest vertex",
                   "[graph][edge_list]", std::uint64_t, double) {
  using W = TestType;
  // Components: {0,3,5}, {1}, {2,4}, {6} (loop only), {7}.
  const EdgeList<W> g{
      8, {{5, 3, W{1}}, {4, 2, W{0}}, {3, 0, W{2}}, {6, 6, W{1}}, {5, 0, W{3}}, {2, 4, W{1}}}};
  const auto c = sssp::connected_components(g);
  REQUIRE(c.count == 5);
  REQUIRE(c.label == std::vector<std::uint32_t>{0, 1, 2, 0, 2, 0, 3, 4});
  REQUIRE(c.size == std::vector<std::size_t>{3, 1, 2, 1, 1});
}

TEST_CASE("connected_components of empty graphs", "[graph][edge_list]") {
  const auto none = sssp::connected_components(EdgeList<double>{});
  REQUIRE(none.count == 0);
  REQUIRE(none.label.empty());
  REQUIRE(none.size.empty());
  const auto isolated = sssp::connected_components(EdgeList<double>{3, {}});
  REQUIRE(isolated.count == 3);
  REQUIRE(isolated.label == std::vector<std::uint32_t>{0, 1, 2});
  REQUIRE(isolated.size == std::vector<std::size_t>{1, 1, 1});
  REQUIRE_THROWS_AS(sssp::connected_components(EdgeList<double>{1, {{0, 1, 1.0}}}),
                    std::out_of_range);
}

// ------------------------------------------------------------------ component subgraphs

TEMPLATE_TEST_CASE("extract_component returns the induced subgraph with vertex and edge maps",
                   "[graph][edge_list]", std::uint64_t, double) {
  using W = TestType;
  const EdgeList<W> g{
      8, {{5, 3, W{1}}, {4, 2, W{0}}, {3, 0, W{2}}, {6, 6, W{1}}, {5, 0, W{3}}, {2, 4, W{1}}}};
  const auto c = sssp::connected_components(g);

  const auto first = sssp::extract_component(g, c, 0);  // {0,3,5}
  REQUIRE(first.graph.num_vertices == 3);
  REQUIRE(first.original_vertex == std::vector<Vertex>{0, 3, 5});
  REQUIRE(first.original_edge == std::vector<EdgeId>{0, 2, 4});
  require_edge<W>(first.graph.edges[0], 2, 1, W{1});
  require_edge<W>(first.graph.edges[1], 1, 0, W{2});
  require_edge<W>(first.graph.edges[2], 2, 0, W{3});

  const auto loop = sssp::extract_component(g, c, 3);  // {6} with its loop
  REQUIRE(loop.original_vertex == std::vector<Vertex>{6});
  REQUIRE(loop.original_edge == std::vector<EdgeId>{3});
  require_edge<W>(loop.graph.edges[0], 0, 0, W{1});

  REQUIRE_THROWS_AS(sssp::extract_component(g, c, 5), std::out_of_range);
  sssp::ComponentLabels wrong = c;
  wrong.label.pop_back();
  REQUIRE_THROWS_AS(sssp::extract_component(g, wrong, 0), std::invalid_argument);
}

TEMPLATE_TEST_CASE("component_of returns the component that contains a vertex",
                   "[graph][edge_list]", std::uint64_t, double) {
  using W = TestType;
  const EdgeList<W> g{
      8, {{5, 3, W{1}}, {4, 2, W{0}}, {3, 0, W{2}}, {6, 6, W{1}}, {5, 0, W{3}}, {2, 4, W{1}}}};
  const auto r = sssp::component_of(g, 4);
  REQUIRE(r.original_vertex == std::vector<Vertex>{2, 4});
  REQUIRE(r.original_edge == std::vector<EdgeId>{1, 5});
  require_edge<W>(r.graph.edges[0], 1, 0, W{0});
  require_edge<W>(r.graph.edges[1], 0, 1, W{1});
  REQUIRE_THROWS_AS(sssp::component_of(g, 8), std::out_of_range);
}

TEMPLATE_TEST_CASE("largest_component picks the most vertices, ties by the smallest label",
                   "[graph][edge_list]", std::uint64_t, double) {
  using W = TestType;
  // {0,1} and {2,3} both have 2 vertices; the first one (smallest vertex 0) wins.
  const EdgeList<W> tie{5, {{3, 2, W{7}}, {1, 0, W{1}}, {1, 0, W{2}}}};
  const auto t = sssp::largest_component(tie);
  REQUIRE(t.original_vertex == std::vector<Vertex>{0, 1});
  REQUIRE(t.original_edge == std::vector<EdgeId>{1, 2});
  require_edge<W>(t.graph.edges[0], 1, 0, W{1});
  require_edge<W>(t.graph.edges[1], 1, 0, W{2});

  // {1,3,4} is strictly larger than {0,2}.
  const EdgeList<W> g{5, {{0, 2, W{1}}, {4, 1, W{2}}, {3, 4, W{3}}}};
  const auto r = sssp::largest_component(g);
  REQUIRE(r.graph.num_vertices == 3);
  REQUIRE(r.original_vertex == std::vector<Vertex>{1, 3, 4});
  REQUIRE(r.original_edge == std::vector<EdgeId>{1, 2});
  require_edge<W>(r.graph.edges[0], 2, 0, W{2});
  require_edge<W>(r.graph.edges[1], 1, 2, W{3});

  const auto empty = sssp::largest_component(EdgeList<W>{});
  REQUIRE(empty.graph.num_vertices == 0);
  REQUIRE(empty.graph.edges.empty());
  REQUIRE(empty.original_vertex.empty());
}
