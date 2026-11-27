// Property test of the Mode C local-hierarchy SSSP (task P5-PRACTICAL-LOCAL) outside its exact
// domain: "a result is never silently wrong".
//
// local_sssp.hpp ("Numerics") promises exact results where every sum is exact, and otherwise
// EITHER a result that is a valid shortest-path solution OR an exception; it never returns
// labels it could not certify. This file checks that promise on the two kinds of input the
// differential test leaves out:
//   * double weights whose sums round (uniform_real, log_uniform, exponential). A returned
//     result must pass verify_sssp against binary-heap Dijkstra under
//     Tolerance::for_path_hops(2(n-1)) (docs/design/numerics.md Sec 3.2, the rule the benchmark
//     runner applies to non-relaxation algorithms); a refusal must be std::invalid_argument
//     (Lemma 8.5: a label not attained by a tight arc) or std::logic_error.
//   * uint64 weights around 2^62 .. 2^63, where masses and scaled thresholds can overflow
//     although the distances are representable (P3-EVAL-ARITH-FIX). A returned result must
//     equal Bellman-Ford exactly; a refusal must be std::overflow_error.
// How often each outcome occurs is reported, not asserted: it is a property of the components,
// not of this entry point.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "experimental/practical/local_sssp.hpp"
#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/verify.hpp"

using sssp::CsrGraph;
using sssp::Vertex;
namespace pr = sssp::practical;

namespace {

std::vector<pr::PracticalOptions> hierarchy_options() {
  std::vector<pr::PracticalOptions> v;
  for (const auto a : {pr::Assembly::kExactHeight, pr::Assembly::kLambdaGrid}) {
    pr::PracticalOptions o;
    o.assembly = a;
    o.direct_search_max_vertices = 0;
    v.push_back(o);
  }
  return v;
}

}  // namespace

TEST_CASE("practical-local on rounding doubles: a valid result or a refusal, never a wrong one",
          "[property][practical][local]") {
  using sssp::gen::WeightDist;
  std::size_t returned = 0;
  std::size_t refused = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    for (const auto dist :
         {WeightDist::kUniformReal, WeightDist::kLogUniform, WeightDist::kExponential}) {
      for (std::uint64_t seed = 1; seed <= 4; ++seed) {
        sssp::gen::GraphSpec spec;
        spec.type = std::string(info.name);
        spec.n = 40;
        spec.avg_degree = 3.0;
        spec.options.seed = seed;
        spec.options.weights.dist = dist;
        spec.options.weights.lo = dist == WeightDist::kLogUniform ? 1e-6 : 0.0;
        spec.options.weights.hi = dist == WeightDist::kLogUniform ? 1e6 : 10.0;
        const auto g = sssp::to_csr(sssp::gen::generate<double>(spec));
        const std::size_t n = g.num_vertices();
        const auto dij = sssp::baseline::dijkstra_binary_heap(g, 0);
        const auto tol = sssp::Tolerance::for_path_hops(2 * (n - 1));
        for (const auto& o : hierarchy_options()) {
          INFO(info.name << " " << sssp::gen::to_string(dist) << " seed=" << seed
                         << " assembly=" << static_cast<int>(o.assembly));
          try {
            const auto r = pr::local_sssp(g, 0, o);
            const auto rep = sssp::verify_sssp(g, 0, r, &dij.dist, tol);
            INFO(rep.summary());
            REQUIRE(rep.ok);
            ++returned;
          } catch (const std::invalid_argument&) {
            ++refused;
          } catch (const std::logic_error&) {
            ++refused;
          }
        }
      }
    }
  }
  INFO("returned=" << returned << " refused=" << refused);
  CHECK(returned + refused == sssp::gen::generator_catalog().size() * 3 * 4 * 2);
  CHECK(returned > 0);
}

TEST_CASE("practical-local on uint64 weights near 2^63: exact or std::overflow_error",
          "[property][practical][local]") {
  using U64 = std::uint64_t;
  const U64 p62 = U64{1} << 62;
  const U64 p63 = U64{1} << 63;
  struct Case {
    const char* name;
    std::size_t n;
    std::vector<sssp::Edge<U64>> edges;
  };
  const std::vector<Case> cases{
      {"one edge 2^63", 2, {{0, 1, p63}}},
      {"one edge 2^62", 2, {{0, 1, p62}}},
      {"path 2^62, 2^62", 3, {{0, 1, p62}, {1, 2, p62}}},
      {"path 1, 2^63 - 2", 3, {{0, 1, 1}, {1, 2, p63 - 2}}},
      {"triangle 2^62", 3, {{0, 1, p62}, {1, 2, p62}, {0, 2, p62 + 1}}},
      {"light path with heavy chord", 4, {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {0, 3, p63}}},
      {"star 2^61", 5, {{0, 1, p62 / 2}, {0, 2, p62 / 2}, {0, 3, p62 / 2}, {0, 4, p62 / 2}}},
  };
  std::size_t returned = 0;
  std::size_t refused = 0;
  for (const auto& c : cases) {
    const auto g = CsrGraph<U64>::from_edges(c.n, c.edges);
    for (Vertex s = 0; s < c.n; ++s) {
      std::vector<U64> want;
      try {
        want = sssp::baseline::bellman_ford(g, s).dist;
      } catch (const std::overflow_error&) {
        continue;  // a distance is not representable: no reference to compare with
      }
      for (const auto& o : hierarchy_options()) {
        INFO(c.name << " s=" << s << " assembly=" << static_cast<int>(o.assembly));
        try {
          const auto r = pr::local_sssp(g, s, o);
          REQUIRE(r.dist == want);
          REQUIRE(sssp::verify_sssp(g, s, r).ok);
          ++returned;
        } catch (const std::overflow_error&) {
          ++refused;
        }
      }
    }
  }
  INFO("returned=" << returned << " refused=" << refused);
  CHECK(returned + refused > 0);
}
