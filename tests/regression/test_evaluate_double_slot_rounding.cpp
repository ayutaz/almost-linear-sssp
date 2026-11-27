// Regression 20261127-evaluate-double-slot-rounding
// (record: cases/20261127-evaluate-double-slot-rounding.json).
//
// Valid double input: the path 0 -(2^53)- 1 -(1)- 2, source 0. Exact distances (exact_sssp.py)
// are {0, 2^53, 2^53 + 1}; 2^53 + 1 is not a double. Algorithm 4 (Mode A, base case k = 3)
// returns a hierarchy that passes validate_seeded_ordered_hierarchy, but Evaluate (Alg 3) threw
// std::logic_error("evaluate: child slot behind the portion").
// Cause: Alg 2 l.9-10 / [Interpretation Q-27] in src/paper/hierarchy/evaluate.cpp
// (start_child_slot). The child {1, 2} (height 1, unit nu = 1) opens at D0 = 2^53, where
// D0 + nu rounds back to D0, so its first slot [D0, D0 + nu) is empty in double and no longer
// reaches into the parent's portion; the invariant "the current slot always reaches into the
// portion" does not survive rounding (D-19 covers only the origin rule b - j nu).
// The same exception appears at the paper's size (P4-LOCAL-LARGE) on light clusters joined by
// bridges of 2^200 .. 2^202, k = 209953, source 0. Thresholds found while minimising: the
// weights (2^52, 1) and (2^53, 2) pass, (2^52, 0.5) and (2^53, 1) fail (a sum that is not
// representable).
//
// Expected behaviour (lead's decision, option (ii) of P3-EVAL-ARITH-FIX): a slot whose end
// rounds back to its start cannot be evaluated by Alg 2 / Alg 3 at all. The slot sequence makes
// no progress (every later slot is the same empty one), and widening it would let a heavy arc
// (length >= nu) of the slot land at an earlier rank inside the same portion, which the
// correctness argument of Lemma 4.6 (case 1: "the crossing costs less than v - u <= nu(C)")
// excludes. Such an input is outside the double scope of Q-06 / Q-55 (option A: correctness is
// established where the arithmetic is exact, i.e. integer and dyadic lengths whose sums stay
// representable), and Evaluate must say so with std::domain_error, not with an internal
// std::logic_error and never with a silent wrong value.
//
// Where no slot rounds away, a non-representable distance is still evaluated and must meet the
// project's double tolerance (docs/design/numerics.md, correctness-review skill §3): each least
// emission d satisfies |d - d*| <= k eps d* with k the hop count of the exact shortest path and
// eps = 2^-53. Example: the same path from source 2, where the unit-1 child opens at 0 and the
// leaf 0 is emitted in the root's slot [2^53, 2^54). Inputs whose sums are representable keep
// their exact results (last test case).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "sssp/graph.hpp"

namespace {

using W = double;
using Catch::Matchers::ContainsSubstring;

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

TEST_CASE("20261127-evaluate-double-slot-rounding: Alg 4 + Evaluate on 0 -2^53- 1 -1- 2",
          "[regression][paper]") {
  const auto g = sssp::CsrGraph<W>::from_edges(3, {{0, 1, 0x1p53}, {1, 2, 1.0}});
  CHECK_THROWS_AS(least_emissions(g, 0), std::domain_error);
  CHECK_THROWS_WITH(least_emissions(g, 0), ContainsSubstring("rounds back to its start"));
  // The same with (2^52, 0.5): 2^52 + 0.5 is not a double either.
  const auto h = sssp::CsrGraph<W>::from_edges(3, {{0, 1, 0x1p52}, {1, 2, 0.5}});
  CHECK_THROWS_AS(least_emissions(h, 0), std::domain_error);
  CHECK_THROWS_WITH(least_emissions(h, 0), ContainsSubstring("rounds back to its start"));
}

TEST_CASE("20261127-evaluate-double-slot-rounding: the same path from source 2 (no empty slot)",
          "[regression][paper]") {
  const auto g = sssp::CsrGraph<W>::from_edges(3, {{0, 1, 0x1p53}, {1, 2, 1.0}});
  std::vector<W> got;
  REQUIRE_NOTHROW(got = least_emissions(g, 2));
  REQUIRE(got.size() == 3);
  // Exact d* = {2^53 + 1, 1, 0}, hops {2, 1, 0}. d* of vertex 0 is formed in long double
  // (64-bit significand on x86, otherwise double: then the bound is checked against 2^53, one
  // unit below d*, which only tightens it by 2^-53 relative).
  constexpr double kEps = 0x1p-53;
  const long double d0 = 0x1p53L + 1.0L;
  CHECK(std::fabs(static_cast<long double>(got[0]) - d0) <= 2.0L * kEps * d0);
  CHECK(got[1] == 1.0);
  CHECK(got[2] == 0.0);
}

TEST_CASE("20261127-evaluate-double-slot-rounding: representable sums (pass before the fix)",
          "[regression][paper]") {
  const auto a = sssp::CsrGraph<W>::from_edges(3, {{0, 1, 0x1p52}, {1, 2, 1.0}});
  CHECK(least_emissions(a, 0) == std::vector<W>{0.0, 0x1p52, 0x1p52 + 1.0});
  const auto b = sssp::CsrGraph<W>::from_edges(3, {{0, 1, 0x1p53}, {1, 2, 2.0}});
  CHECK(least_emissions(b, 0) == std::vector<W>{0.0, 0x1p53, 0x1p53 + 2.0});
}
