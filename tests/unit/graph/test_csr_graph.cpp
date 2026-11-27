#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "sssp/graph.hpp"

using sssp::CsrGraph;
using sssp::Edge;

TEST_CASE("CSR stores both orientations of every undirected edge", "[graph][csr]") {
  const std::vector<Edge<std::uint64_t>> edges{{0, 1, 5}, {1, 2, 7}};
  const auto g = CsrGraph<std::uint64_t>::from_edges(3, edges);
  REQUIRE(g.num_vertices() == 3);
  REQUIRE(g.num_edges() == 2);
  REQUIRE(g.num_arcs() == 4);
  REQUIRE(g.degree(1) == 2);
  std::uint64_t sum = 0;
  for (const auto& a : g.arcs(1)) {
    sum += a.weight;
    REQUIRE(g.edge(a.edge).weight == a.weight);
  }
  REQUIRE(sum == 12);
}

TEST_CASE("CSR keeps multi-edges, self-loops and isolated vertices", "[graph][csr]") {
  const std::vector<Edge<double>> edges{{0, 1, 1.0}, {0, 1, 2.0}, {2, 2, 3.0}};
  const auto g = CsrGraph<double>::from_edges(4, edges);
  REQUIRE(g.num_edges() == 3);
  REQUIRE(g.degree(0) == 2);
  REQUIRE(g.degree(2) == 1);  // a self-loop is one arc (u -> u)
  REQUIRE(g.degree(3) == 0);
}

TEST_CASE("CSR rejects invalid input with clear errors", "[graph][csr][validation]") {
  using W = double;
  REQUIRE_THROWS_AS(CsrGraph<W>::from_edges(2, {{0, 2, 1.0}}), std::out_of_range);
  REQUIRE_THROWS_AS(CsrGraph<W>::from_edges(2, {{0, 1, -1.0}}), std::invalid_argument);
  REQUIRE_THROWS_AS(CsrGraph<W>::from_edges(2, {{0, 1, std::nan("")}}), std::invalid_argument);
  REQUIRE_THROWS_AS(CsrGraph<W>::from_edges(2, {{0, 1, std::numeric_limits<W>::infinity()}}),
                    std::invalid_argument);
  REQUIRE_NOTHROW(CsrGraph<W>::from_edges(2, {{0, 1, 0.0}}));  // zero weight is allowed
}

TEST_CASE("weight traits define infinity per weight type", "[graph][weights]") {
  REQUIRE(sssp::WeightTraits<std::uint64_t>::infinity() ==
          std::numeric_limits<std::uint64_t>::max());
  REQUIRE(std::isinf(sssp::WeightTraits<double>::infinity()));
  REQUIRE(sssp::WeightTraits<std::uint64_t>::is_exact);
  REQUIRE_FALSE(sssp::WeightTraits<double>::is_exact);
}

TEST_CASE("integer distance addition detects overflow", "[graph][weights]") {
  using T = sssp::WeightTraits<std::uint64_t>;
  REQUIRE(T::checked_add(2, 3) == 5);
  REQUIRE_THROWS_AS(T::checked_add(std::numeric_limits<std::uint64_t>::max() - 1, 2),
                    std::overflow_error);
}
