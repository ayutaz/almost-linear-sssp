// Regression 20261127-double-sum-overflow (record: cases/20261127-double-sum-overflow.json).
//
// Valid double input (finite, non-negative weights): path 0 -1e308- 1 -1e308- 2. The true
// distance of vertex 2 is 2e308 > DBL_MAX, so fl(1e308 + 1e308) = +inf, the "unreachable"
// sentinel. bellman_ford reports this as std::overflow_error ("exceeds the representable range").
// Every double Dijkstra instead returns dist[2] = inf with no parent for a REACHABLE vertex: a
// result that verify_sssp rejects (reachable_but_infinite) and that silently claims vertex 2 is
// unreachable. Cause: WeightTraits<double>::checked_add (include/sssp/graph.hpp) is a plain
// addition, and the strict test nd < dist[v] drops the infinite candidate without a trace.
// Expected behaviour chosen here: the same as bellman_ford and as the integer policy in
// docs/design/baselines.md ("the result is never silently wrong") -> throw std::overflow_error.
//
// KNOWN BUG, NOT YET FIXED: tagged [!shouldfail] so the suite stays green; Catch2 reports the test
// as FAILED once the baselines are fixed, and the fix commit must remove the tag (RED is recorded
// in docs/tdd/tasks/P2-DIFF-DOUBLE-OVERFLOW.json).

#include <catch2/catch_test_macros.hpp>
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
  // Reference behaviour (passes today): the oracle refuses to return a wrong result.
  CHECK_THROWS_AS(sssp::baseline::bellman_ford(overflowing_path(), 0), std::overflow_error);
}

TEST_CASE(
    "20261127-double-sum-overflow: Dijkstra baselines must not return inf for a reachable "
    "vertex",
    "[regression][baselines][!shouldfail]") {
  const auto g = overflowing_path();
  CHECK_THROWS_AS(sssp::baseline::dijkstra_binary_heap(g, 0), std::overflow_error);
  CHECK_THROWS_AS(sssp::baseline::dijkstra_4ary_heap(g, 0), std::overflow_error);
  CHECK_THROWS_AS(sssp::baseline::dijkstra_fibonacci_heap(g, 0), std::overflow_error);
}
