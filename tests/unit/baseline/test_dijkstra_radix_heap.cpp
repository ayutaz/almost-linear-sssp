#include <catch2/catch_test_macros.hpp>
#include <cstdint>

#include "sssp/baseline.hpp"
#include "unit/baseline/dijkstra_test_support.hpp"

using namespace sssp::test;

namespace {
// The radix heap needs monotone integer keys, so only std::uint64_t weights exist.
struct Algo {
  sssp::SsspResult<std::uint64_t> operator()(const sssp::CsrGraph<std::uint64_t>& g,
                                             sssp::Vertex s) const {
    return sssp::baseline::dijkstra_radix_heap(g, s);
  }
};
}  // namespace

TEST_CASE("radix-heap Dijkstra: hand-computed classic graph", "[dijkstra][radix]") {
  expect_classic_graph<std::uint64_t>(Algo{});
}

TEST_CASE("radix-heap Dijkstra: zero weights, multi-edges, self-loops, components",
          "[dijkstra][radix]") {
  expect_special_graph<std::uint64_t>(Algo{});
  expect_zero_weight_cycle<std::uint64_t>(Algo{});
}

TEST_CASE("radix-heap Dijkstra: single vertex", "[dijkstra][radix]") {
  expect_single_vertex<std::uint64_t>(Algo{});
}

TEST_CASE("radix-heap Dijkstra: invalid source and integer overflow", "[dijkstra][radix]") {
  expect_invalid_source_rejected<std::uint64_t>(Algo{});
  expect_integer_overflow_detected(Algo{});
}

TEST_CASE("radix-heap Dijkstra: matches Floyd-Warshall on random graphs",
          "[dijkstra][radix][random]") {
  expect_matches_floyd_warshall<std::uint64_t>(Algo{}, int_weights(20), true);
  expect_matches_floyd_warshall<std::uint64_t>(Algo{}, int_weights(std::uint64_t{1} << 40), true);
  // Near the top of the 64-bit range (buckets 60..64 are exercised).
  expect_matches_floyd_warshall<std::uint64_t>(Algo{}, int_weights(std::uint64_t{1} << 55), true);
}

#if SSSP_DIAGNOSTICS
TEST_CASE("radix-heap Dijkstra: operation counters", "[dijkstra][radix][diagnostics]") {
  expect_classic_counters<std::uint64_t>(Algo{}, true);
}
#endif
