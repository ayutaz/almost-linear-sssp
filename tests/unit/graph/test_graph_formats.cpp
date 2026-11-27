#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"

using sssp::AdjacencyList;
using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeList;

namespace {

template <class W>
EdgeList<W> sample_edge_list() {
  // multi-edge (0,1) twice, a self-loop on 2, an isolated vertex 4
  return EdgeList<W>{5, {{0, 1, W{3}}, {1, 2, W{4}}, {0, 1, W{1}}, {2, 2, W{7}}, {3, 1, W{2}}}};
}

template <class W>
bool same_edges(std::span<const Edge<W>> a, std::span<const Edge<W>> b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i].u != b[i].u || a[i].v != b[i].v || a[i].weight != b[i].weight) return false;
  }
  return true;
}

}  // namespace

TEST_CASE("EdgeList -> CSR -> EdgeList preserves vertices and edge order", "[graph][formats]") {
  const auto el = sample_edge_list<std::uint64_t>();
  const CsrGraph<std::uint64_t> g = sssp::to_csr(el);
  REQUIRE(g.num_vertices() == 5);
  REQUIRE(g.num_edges() == 5);
  const auto back = sssp::to_edge_list(g);
  REQUIRE(back.num_vertices == 5);
  REQUIRE(same_edges<std::uint64_t>(back.edges, el.edges));
}

TEST_CASE("to_csr validates the edge list", "[graph][formats][validation]") {
  REQUIRE_THROWS_AS(sssp::to_csr(EdgeList<double>{2, {{0, 2, 1.0}}}), std::out_of_range);
  REQUIRE_THROWS_AS(sssp::to_csr(EdgeList<double>{2, {{0, 1, -1.0}}}), std::invalid_argument);
}

TEST_CASE("AdjacencyList built from CSR has identical arcs in identical order",
          "[graph][formats]") {
  const auto g = sssp::to_csr(sample_edge_list<double>());
  const AdjacencyList<double> adj = sssp::to_adjacency_list(g);
  REQUIRE(adj.num_vertices() == g.num_vertices());
  REQUIRE(adj.num_edges() == g.num_edges());
  REQUIRE(adj.num_arcs() == g.num_arcs());
  for (sssp::Vertex v = 0; v < g.num_vertices(); ++v) {
    const auto a = adj.arcs(v);
    const auto c = g.arcs(v);
    REQUIRE(a.size() == c.size());
    REQUIRE(adj.degree(v) == g.degree(v));
    for (std::size_t i = 0; i < a.size(); ++i) {
      REQUIRE(a[i].to == c[i].to);
      REQUIRE(a[i].weight == c[i].weight);
      REQUIRE(a[i].edge == c[i].edge);
    }
  }
  REQUIRE(adj.degree(4) == 0);
  REQUIRE(adj.degree(2) == 2);  // edge {1,2} plus the self-loop as one arc
}

TEST_CASE("AdjacencyList built from EdgeList matches the CSR arc order", "[graph][formats]") {
  const auto el = sample_edge_list<std::uint64_t>();
  const auto direct = sssp::to_adjacency_list(el);
  const auto via_csr = sssp::to_adjacency_list(sssp::to_csr(el));
  REQUIRE(direct.num_vertices() == via_csr.num_vertices());
  for (sssp::Vertex v = 0; v < direct.num_vertices(); ++v) {
    const auto a = direct.arcs(v);
    const auto b = via_csr.arcs(v);
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      REQUIRE(a[i].to == b[i].to);
      REQUIRE(a[i].edge == b[i].edge);
    }
  }
}

TEST_CASE("AdjacencyList round-trips through CSR and EdgeList", "[graph][formats]") {
  const auto el = sample_edge_list<std::uint64_t>();
  const auto adj = sssp::to_adjacency_list(el);
  REQUIRE(same_edges<std::uint64_t>(sssp::to_edge_list(adj).edges, el.edges));
  const auto g = sssp::to_csr(adj);
  REQUIRE(same_edges<std::uint64_t>(g.edges(), el.edges));
  REQUIRE(g.num_arcs() == adj.num_arcs());
}

TEST_CASE("empty graphs convert", "[graph][formats]") {
  const EdgeList<double> none{0, {}};
  REQUIRE(sssp::to_csr(none).num_vertices() == 0);
  REQUIRE(sssp::to_adjacency_list(none).num_vertices() == 0);
  const EdgeList<double> isolated{3, {}};
  const auto adj = sssp::to_adjacency_list(sssp::to_csr(isolated));
  REQUIRE(adj.num_vertices() == 3);
  REQUIRE(adj.num_arcs() == 0);
}

TEST_CASE("timed conversion reports a non-negative duration and the result",
          "[graph][formats][timing]") {
  const auto el = sample_edge_list<std::uint64_t>();
  const auto timed = sssp::timed([&] { return sssp::to_csr(el); });
  REQUIRE(timed.seconds >= 0.0);
  REQUIRE(timed.value.num_edges() == 5);
}
