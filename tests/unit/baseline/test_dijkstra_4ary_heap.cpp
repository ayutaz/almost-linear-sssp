#include <catch2/catch_test_macros.hpp>
#include <cstdint>

#include "sssp/baseline.hpp"
#include "unit/baseline/dijkstra_test_support.hpp"

using namespace sssp::test;

namespace {
struct Algo {
  template <class W>
  sssp::SsspResult<W> operator()(const sssp::CsrGraph<W>& g, sssp::Vertex s) const {
    return sssp::baseline::dijkstra_4ary_heap(g, s);
  }
};
}  // namespace

TEST_CASE("4-ary-heap Dijkstra: hand-computed classic graph", "[dijkstra][4ary]") {
  expect_classic_graph<std::uint64_t>(Algo{});
  expect_classic_graph<double>(Algo{});
}

TEST_CASE("4-ary-heap Dijkstra: zero weights, multi-edges, self-loops, components",
          "[dijkstra][4ary]") {
  expect_special_graph<std::uint64_t>(Algo{});
  expect_special_graph<double>(Algo{});
  expect_zero_weight_cycle<std::uint64_t>(Algo{});
  expect_zero_weight_cycle<double>(Algo{});
}

TEST_CASE("4-ary-heap Dijkstra: single vertex", "[dijkstra][4ary]") {
  expect_single_vertex<std::uint64_t>(Algo{});
  expect_single_vertex<double>(Algo{});
}

TEST_CASE("4-ary-heap Dijkstra: exact dyadic doubles", "[dijkstra][4ary]") {
  expect_dyadic_doubles(Algo{});
}

TEST_CASE("4-ary-heap Dijkstra: invalid source and integer overflow", "[dijkstra][4ary]") {
  expect_invalid_source_rejected<std::uint64_t>(Algo{});
  expect_invalid_source_rejected<double>(Algo{});
  expect_integer_overflow_detected(Algo{});
}

TEST_CASE("4-ary-heap Dijkstra: matches Floyd-Warshall on random graphs",
          "[dijkstra][4ary][random]") {
  expect_matches_floyd_warshall<std::uint64_t>(Algo{}, int_weights(20), true);
  expect_matches_floyd_warshall<std::uint64_t>(Algo{}, int_weights(std::uint64_t{1} << 40), true);
  expect_matches_floyd_warshall<double>(Algo{}, dyadic_weights(), true);
  expect_matches_floyd_warshall<double>(Algo{}, real_weights(), false);
}

#if SSSP_DIAGNOSTICS
TEST_CASE("4-ary-heap Dijkstra: operation counters", "[dijkstra][4ary][diagnostics]") {
  expect_classic_counters<std::uint64_t>(Algo{}, true);
  expect_classic_counters<double>(Algo{}, true);
}
#endif
