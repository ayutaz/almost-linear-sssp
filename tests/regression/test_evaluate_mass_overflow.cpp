// Regression 20261127-evaluate-mass-overflow (record: cases/20261127-evaluate-mass-overflow.json)
// and 20261127-evaluate-relax-overflow (record: cases/20261127-evaluate-relax-overflow.json).
//
// Valid uint64 input: one edge {0, 1} of length 2^63. Every distance fits (d(0, 1) = 2^63 <
// 2^64 - 1, exact reference tests/python/reference/exact_sssp.py). Algorithm 4 (Mode A, base
// case k = 2 <= 2 L_k^4) returns a hierarchy that passes validate_seeded_ordered_hierarchy, but
// Evaluate (Alg 3) throws std::overflow_error("integer length overflow").
// Suspected step: the mass S(C) = h(C) + lengths of C's links + S of internal children of
// Sec 4.1 (l.686-688, src/paper/hierarchy/normalize.hpp), here 2^63 + 2^63 = 2^64, or the expiry
// c = D0 + S(X) of Alg 3 l.9. These are bookkeeping bounds, not distances: Lemma 4.6 promises
// delta at protected states, and delta fits the type. At the paper's size (P4-LOCAL-LARGE) the
// same failure appears on lt::clustered(209953, 107) (bridges of 2^54 .. 2^56, MST total about
// 2^70) for source 0.
//
// Expected behaviour (from Lemma 4.6 and the overflow policy of docs/design/numerics.md §6 /
// include/sssp/paper.hpp: an integer overflow is an error only if a DISTANCE would reach the
// infinity sentinel): Evaluate returns, for each source, least emissions equal to the exact
// distances. An input whose true distance does not fit still throws std::overflow_error (never
// a silent "unreachable").
//
// Found by: tests/differential/paper/test_local_hierarchy_large_vs_dijkstra.cpp (P4-LOCAL-LARGE),
// minimised by hand to one edge; 2^63 - 1 passes, 2^63 fails.
//
// 20261127-evaluate-relax-overflow (found while searching the same class, P3-EVAL-ARITH-FIX):
// the relaxation value emitted(tail) + l of Alg 2 l.11 overflows on a walk that is not a
// shortest path. Edges 0 -1- 1, 1 -(2^64 - 2)- 2, 1 -1- 3, 3 -1- 2, source 0: when vertex 1
// emits at 1, the arc 1 -> 2 gives 1 + 2^64 - 2 = 2^64 - 1 (the sentinel) while vertex 2 still
// has no offer; its distance is 3 (via 3). Evaluate threw std::overflow_error. Expected: the
// exact distances {0, 1, 3, 2} (Bellman-Ford's rule in numerics.md §6: an overflowing candidate
// is skipped; only an unrepresentable distance is an error).

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "sssp/graph.hpp"

namespace {

using W = std::uint64_t;

std::vector<W> least_emissions(const sssp::CsrGraph<W>& g, sssp::Vertex src) {
  const auto lh = sssp::paper::local_hierarchy(g);
  sssp::paper::SeededOrderedHierarchy<W> s{
      lh.hierarchy, {{src, W{0}}}, {{lh.native_entry[src], W{0}}}};
  const auto valid = sssp::paper::validate_seeded_ordered_hierarchy(s, g);
  INFO(valid.summary());
  REQUIRE(valid.ok());
  const auto r = sssp::paper::evaluate(s);
  std::vector<W> best(g.num_vertices(), sssp::WeightTraits<W>::infinity());
  for (const auto& [x, v] : r.emissions) {
    const auto p = lh.hierarchy.projection[x];
    if (v < best[p]) best[p] = v;
  }
  return best;
}

}  // namespace

TEST_CASE("20261127-evaluate-mass-overflow: Alg 4 + Evaluate on one edge of length 2^63",
          "[regression][paper]") {
  constexpr W kLen = W{1} << 63;
  const auto g = sssp::CsrGraph<W>::from_edges(2, {{0, 1, kLen}});
  // Exact reference (exact_sssp.py): source 0 -> {0, 2^63}; source 1 -> {2^63, 0}.
  std::vector<W> got0;
  REQUIRE_NOTHROW(got0 = least_emissions(g, 0));
  CHECK(got0 == std::vector<W>{0, kLen});
  std::vector<W> got1;
  REQUIRE_NOTHROW(got1 = least_emissions(g, 1));
  CHECK(got1 == std::vector<W>{kLen, 0});
}

TEST_CASE("20261127-evaluate-mass-overflow: one edge of length 2^63 - 1 (passes before the fix)",
          "[regression][paper]") {
  constexpr W kLen = (W{1} << 63) - 1;
  const auto g = sssp::CsrGraph<W>::from_edges(2, {{0, 1, kLen}});
  CHECK(least_emissions(g, 0) == std::vector<W>{0, kLen});
}

TEST_CASE("20261127-evaluate-mass-overflow: a mass of 2^64 with slots up to the sentinel",
          "[regression][paper]") {
  // Path 1 -q- 0 -q- 2 -q- 3 with q = 2^62: one base-case node of height q with three links,
  // S = 4q = 2^64. From vertex 1 the distances are {q, 0, 2q, 3q}; the root's slots of width q
  // reach [3q, 4q), whose end is the sentinel + 1.
  constexpr W q = W{1} << 62;
  const auto g = sssp::CsrGraph<W>::from_edges(4, {{0, 1, q}, {0, 2, q}, {2, 3, q}});
  std::vector<W> got;
  REQUIRE_NOTHROW(got = least_emissions(g, 1));
  CHECK(got == std::vector<W>{q, 0, 2 * q, 3 * q});
}

TEST_CASE("20261127-evaluate-mass-overflow: a distance of 2^64 still throws std::overflow_error",
          "[regression][paper]") {
  // 1 -2^63- 0 -2^63- 2 from vertex 1: d(2) = 2^64 does not fit. The masses overflow as well;
  // whatever saturates, the unrepresentable distance must not come out as "unreachable".
  constexpr W kLen = W{1} << 63;
  const auto g = sssp::CsrGraph<W>::from_edges(3, {{0, 1, kLen}, {0, 2, kLen}});
  CHECK_THROWS_AS(least_emissions(g, 1), std::overflow_error);
}

TEST_CASE("20261127-evaluate-relax-overflow: an overflowing relaxation on a non-shortest walk",
          "[regression][paper]") {
  constexpr W kBig = sssp::WeightTraits<W>::infinity() - 1;  // 2^64 - 2, a valid length
  const auto g = sssp::CsrGraph<W>::from_edges(4, {{0, 1, 1}, {1, 2, kBig}, {1, 3, 1}, {3, 2, 1}});
  // Exact reference: {0, 1, 3, 2}.
  std::vector<W> got;
  REQUIRE_NOTHROW(got = least_emissions(g, 0));
  CHECK(got == std::vector<W>{0, 1, 3, 2});
}
