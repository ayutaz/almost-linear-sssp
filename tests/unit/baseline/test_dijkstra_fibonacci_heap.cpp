#include <catch2/catch_test_macros.hpp>
#include <cstdint>

#include "sssp/baseline.hpp"
#include "unit/baseline/dijkstra_test_support.hpp"

using namespace sssp::test;

namespace {
struct Algo {
  template <class W>
  sssp::SsspResult<W> operator()(const sssp::CsrGraph<W>& g, sssp::Vertex s) const {
    return sssp::baseline::dijkstra_fibonacci_heap(g, s);
  }
};
}  // namespace

TEST_CASE("Fibonacci-heap Dijkstra: hand-computed classic graph", "[dijkstra][fibonacci]") {
  expect_classic_graph<std::uint64_t>(Algo{});
  expect_classic_graph<double>(Algo{});
}

TEST_CASE("Fibonacci-heap Dijkstra: zero weights, multi-edges, self-loops, components",
          "[dijkstra][fibonacci]") {
  expect_special_graph<std::uint64_t>(Algo{});
  expect_special_graph<double>(Algo{});
  expect_zero_weight_cycle<std::uint64_t>(Algo{});
  expect_zero_weight_cycle<double>(Algo{});
}

TEST_CASE("Fibonacci-heap Dijkstra: single vertex", "[dijkstra][fibonacci]") {
  expect_single_vertex<std::uint64_t>(Algo{});
  expect_single_vertex<double>(Algo{});
}

TEST_CASE("Fibonacci-heap Dijkstra: exact dyadic doubles", "[dijkstra][fibonacci]") {
  expect_dyadic_doubles(Algo{});
}

TEST_CASE("Fibonacci-heap Dijkstra: invalid source and integer overflow", "[dijkstra][fibonacci]") {
  expect_invalid_source_rejected<std::uint64_t>(Algo{});
  expect_invalid_source_rejected<double>(Algo{});
  expect_integer_overflow_detected(Algo{});
}

TEST_CASE("Fibonacci-heap Dijkstra: matches Floyd-Warshall on random graphs",
          "[dijkstra][fibonacci][random]") {
  expect_matches_floyd_warshall<std::uint64_t>(Algo{}, int_weights(20), true);
  expect_matches_floyd_warshall<std::uint64_t>(Algo{}, int_weights(std::uint64_t{1} << 40), true);
  expect_matches_floyd_warshall<double>(Algo{}, dyadic_weights(), true);
  expect_matches_floyd_warshall<double>(Algo{}, real_weights(), false);
}

#if SSSP_DIAGNOSTICS
TEST_CASE("Fibonacci-heap Dijkstra: operation counters", "[dijkstra][fibonacci][diagnostics]") {
  expect_classic_counters<std::uint64_t>(Algo{}, true);
  expect_classic_counters<double>(Algo{}, true);
}
#endif
