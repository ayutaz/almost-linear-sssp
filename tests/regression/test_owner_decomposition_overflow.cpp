// Regression 20261127-owner-decomposition-overflow
// (record: cases/20261127-owner-decomposition-overflow.json).
//
// Valid uint64 input: the path 0 - 1 - 2 - 3 with edges of length 2^63 - 1, 2^62, 2^62. From
// vertex 1 the distances are 2^63 - 1, 0, 2^62, 2^63 (exact reference
// tests/python/reference/exact_sssp.py), and every tentative length d[u] + w of an edge is at
// most 2^64 - 2, below the infinity sentinel: the baselines solve it. Sec 5.1
// (owner_decomposition) threw std::overflow_error("integer distance overflow") for every
// separator size b.
//
// Cause: the reduced-graph test w <= D_T(u, v) of l.965-966 is decided from the distances of
// the MST from its root (vertex 0, Lemma 3.9). They were accumulated with
// WeightTraits::checked_add, and D_T(0, 3) = 2^64 - 1 is the infinity sentinel. A root distance
// is an auxiliary quantity, not a tentative path length from the source: by the overflow rule it
// must not raise. The same holds for R_v = D_T(v, p(v)) when the cell is deep (b = 5 below: the
// single pivot is vertex 0, so R_3 = 2^64 - 1).
//
// Expected: the decomposition is returned; root distances are exact (128-bit), and an R_v that
// is not below the sentinel is the infinity sentinel.
//
// Found by: review of the overflow rule after P3-CHILDORDER-OVERFLOW; reported on the path of
// four edges of 2^62 with b = 2 (kept below). Minimised to four vertices and to the least sum:
// with three vertices 0 - 1 - 2 (weights a, c) a root distance reaches the sentinel only if
// a + c >= 2^64 - 1, and then one of the tentative lengths 2a, 2c (source 1) or a + c (source 0
// or 2) reaches it too, so a baseline may raise as well; with four vertices the sum
// (2^63 - 1) + 2^62 + 2^62 is exactly the sentinel.

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/verify.hpp"

namespace {

using W = std::uint64_t;
using sssp::Vertex;
using sssp::paper::owner_decomposition;
using sssp::paper::OwnerDecomposition;

constexpr W kInf = sssp::WeightTraits<W>::infinity();
constexpr W kQ = W{1} << 62;
constexpr W kA = (W{1} << 63) - 1;

sssp::CsrGraph<W> minimal_path() {
  return sssp::CsrGraph<W>::from_edges(4, {{0, 1, kA}, {1, 2, kQ}, {2, 3, kQ}});
}

sssp::CsrGraph<W> uniform_path(std::size_t n, W len) {
  std::vector<sssp::Edge<W>> edges;
  for (Vertex v = 0; v + 1 < n; ++v) edges.push_back({v, v + 1, len});
  return sssp::CsrGraph<W>::from_edges(n, edges);
}

}  // namespace

TEST_CASE("20261127-owner-decomposition-overflow: every distance of the input fits",
          "[regression][paper]") {
  // The premise (passed before the fix): expected_dist of the record, from exact_sssp.py.
  const auto g = minimal_path();
  const std::vector<W> expected{kA, 0, kQ, 2 * kQ};
  const auto r = sssp::baseline::dijkstra_binary_heap(g, 1);
  CHECK(r.dist == expected);
  const auto report = sssp::verify_sssp(g, 1, r, &expected);
  INFO(report.summary());
  CHECK(report.ok);
}

TEST_CASE("20261127-owner-decomposition-overflow: path of 2^63 - 1, 2^62, 2^62 with b = 2",
          "[regression][paper]") {
  const auto g = minimal_path();
  OwnerDecomposition<W> od;
  REQUIRE_NOTHROW(od = owner_decomposition(g, 2));
  // Postorder from root 0: 3 (1), 2 (2: mark), 1 (1), 0 (2: mark). Corridor 0 - 1 - 2 is cut at
  // its maximum, edge 0.
  CHECK(od.pivots == std::vector<Vertex>{0, 2});
  CHECK(od.cut == std::vector<std::uint8_t>{1, 0, 0});
  CHECK(od.owner == std::vector<std::uint32_t>{0, 1, 1, 1});
  CHECK(od.reduced == std::vector<std::uint8_t>{1, 1, 1});
  CHECK(od.radius == std::vector<W>{0, kQ, 0, kQ});
  CHECK(od.owner_distance == std::vector<W>{0, kQ, 0, kQ});
}

TEST_CASE("20261127-owner-decomposition-overflow: the same path with b = 5 (R_3 does not fit)",
          "[regression][paper]") {
  const auto g = minimal_path();
  OwnerDecomposition<W> od;
  REQUIRE_NOTHROW(od = owner_decomposition(g, 5));
  // b > N: no mark, the root is the single pivot (l.971, Q-31).
  CHECK(od.pivots == std::vector<Vertex>{0});
  CHECK(od.radius == std::vector<W>{0, kA, kA, kA});
  // R_3 = (2^63 - 1) + 2^62 + 2^62 = 2^64 - 1 is not below the sentinel.
  CHECK(od.owner_distance == std::vector<W>{0, kA, kA + kQ, kInf});
}

TEST_CASE("20261127-owner-decomposition-overflow: path of four edges of 2^62 with b = 2",
          "[regression][paper]") {
  // The reported case; from vertex 2 the distances are 2^63, 2^62, 0, 2^62, 2^63 (exact_sssp.py).
  const auto g = uniform_path(5, kQ);
  const std::vector<W> expected{2 * kQ, kQ, 0, kQ, 2 * kQ};
  CHECK(sssp::baseline::dijkstra_binary_heap(g, 2).dist == expected);

  OwnerDecomposition<W> od;
  REQUIRE_NOTHROW(od = owner_decomposition(g, 2));
  CHECK(od.pivots == std::vector<Vertex>{1, 3});
  CHECK(od.owner == std::vector<std::uint32_t>{0, 0, 1, 1, 1});
  CHECK(od.reduced == std::vector<std::uint8_t>{1, 1, 1, 1});
  CHECK(od.owner_distance == std::vector<W>{kQ, 0, kQ, 0, kQ});
}

TEST_CASE(
    "20261127-owner-decomposition-overflow: path of four edges of 2^61 (passes before the fix)",
    "[regression][paper]") {
  constexpr W q = W{1} << 61;
  const auto od = owner_decomposition(uniform_path(5, q), 2);
  CHECK(od.pivots == std::vector<Vertex>{1, 3});
  CHECK(od.owner_distance == std::vector<W>{q, 0, q, 0, q});
}
