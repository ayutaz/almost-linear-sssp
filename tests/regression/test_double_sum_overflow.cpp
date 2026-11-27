// Regression 20261127-double-sum-overflow (record: cases/20261127-double-sum-overflow.json).
//
// Valid double input (finite, non-negative weights): path 0 -1e308- 1 -1e308- 2. The true
// distance of vertex 2 is 2e308 > DBL_MAX, so fl(1e308 + 1e308) = +inf, the "unreachable"
// sentinel. bellman_ford reports this as std::overflow_error ("exceeds the representable range").
// Every double Dijkstra used to return dist[2] = inf with no parent for a REACHABLE vertex: a
// result that verify_sssp rejects (reachable_but_infinite) and that silently claimed vertex 2
// was unreachable. Cause: WeightTraits<double>::checked_add (include/sssp/graph.hpp) was a plain
// addition, and the strict test nd < dist[v] dropped the infinite candidate without a trace.
// Expected behaviour chosen here: the same as bellman_ford and as the integer policy in
// docs/design/baselines.md ("the result is never silently wrong") -> throw std::overflow_error.
//
// FIXED (task P2-DIFF-DOUBLE-OVERFLOW): WeightTraits<double>::checked_add throws
// std::overflow_error when the rounded sum of two finite distances is not finite, as the uint64
// version throws when the sum would reach infinity(). Broader coverage (every baseline against
// an independent rounded oracle): tests/differential/test_double_overflow_baselines.cpp.

#include <catch2/catch_test_macros.hpp>
#include <limits>
#include <stdexcept>

#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"

namespace {

sssp::CsrGraph<double> overflowing_path() {
  return sssp::CsrGraph<double>::from_edges(3, {{0, 1, 1e308}, {1, 2, 1e308}});
}

}  // namespace

TEST_CASE("20261127-double-sum-overflow: bellman_ford reports the unrepresentable distance",
          "[regression][baselines]") {
  // Reference behaviour (passed before the fix): the oracle refuses to return a wrong result.
  CHECK_THROWS_AS(sssp::baseline::bellman_ford(overflowing_path(), 0), std::overflow_error);
}

TEST_CASE(
    "20261127-double-sum-overflow: Dijkstra baselines must not return inf for a reachable "
    "vertex",
    "[regression][baselines]") {
  const auto g = overflowing_path();
  CHECK_THROWS_AS(sssp::baseline::dijkstra_binary_heap(g, 0), std::overflow_error);
  CHECK_THROWS_AS(sssp::baseline::dijkstra_4ary_heap(g, 0), std::overflow_error);
  CHECK_THROWS_AS(sssp::baseline::dijkstra_fibonacci_heap(g, 0), std::overflow_error);
}

TEST_CASE("20261127-double-sum-overflow: double checked_add throws iff the rounded sum is inf",
          "[regression][baselines]") {
  using T = sssp::WeightTraits<double>;
  constexpr double kMax = std::numeric_limits<double>::max();
  CHECK_THROWS_AS(T::checked_add(1e308, 1e308), std::overflow_error);
  CHECK_THROWS_AS(T::checked_add(kMax, kMax), std::overflow_error);
  // ulp(DBL_MAX) = 2^971: adding half an ulp rounds to +inf, a quarter ulp is absorbed.
  CHECK_THROWS_AS(T::checked_add(kMax, 0x1p970), std::overflow_error);
  CHECK(T::checked_add(kMax, 0x1p969) == kMax);
  CHECK(T::checked_add(kMax / 2, kMax / 2) == kMax);
  CHECK(T::checked_add(kMax, 0.0) == kMax);
  CHECK(T::checked_add(0.5, 0.25) == 0.75);
}
