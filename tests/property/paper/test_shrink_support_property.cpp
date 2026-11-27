// The shrinking driver of shrink_support.hpp (P4-PROPERTY: "shrinking reproduces failures with
// a seed") [P4-PROPERTY-LEMMAS].
//
// A property test is only useful when its failures can be replayed. These tests plant failing
// properties and check that the driver (1) stops at the first failing seed, (2) regenerates the
// same original input from that seed, (3) shrinks it deterministically (two runs give the same
// minimal case), (4) to a case that still fails and is 1-minimal (deleting any one edge or
// lowering any one weight makes the property hold), and (5) names the seed and both edge lists
// in its report. The second planted property is a deliberately false bound on the real
// Estimate output (Z <= 2 d_S), so shrinking is exercised on a library component.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/pivots/estimate.hpp"
#include "property/paper/shrink_support.hpp"
#include "sssp/graph.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

using sssp::Edge;
using sssp::Vertex;
namespace sh = sssp::test::shrink;
namespace pv = sssp::paper::pivots;

namespace {

using W = std::uint64_t;

sh::Case<W> random_case(std::uint64_t seed) {
  sssp::test::SplitMix64 g(seed);
  sh::Case<W> c;
  c.n = 3 + g.below(10);
  c.b = 1 + g.below(4);
  const std::size_t m = g.below(3 * c.n);
  for (std::size_t i = 0; i < m; ++i) {
    const auto u = static_cast<Vertex>(g.below(c.n));
    auto v = static_cast<Vertex>(g.below(c.n));
    if (u == v) v = static_cast<Vertex>((v + 1) % c.n);
    c.edges.push_back({u, v, 1 + g.below(9)});
  }
  return c;
}

// Planted property: "the edges of weight <= 3 form a forest". Fails on any light cycle
// (including two parallel light edges).
std::optional<std::string> light_cycle(const sh::Case<W>& c) {
  std::vector<std::size_t> p(c.n);
  for (std::size_t i = 0; i < c.n; ++i) p[i] = i;
  const auto find = [&](std::size_t x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  };
  for (const auto& e : c.edges) {
    if (e.weight > 3) continue;
    const std::size_t a = find(e.u);
    const std::size_t b = find(e.v);
    if (a == b) return "light cycle closed by " + std::to_string(e.u) + "-" + std::to_string(e.v);
    p[a] = b;
  }
  return std::nullopt;
}

// Planted false bound on the library's Estimate: Z <= 2 d_S for every connected demand.
std::optional<std::string> z_at_most_twice(const sh::Case<W>& c) {
  if (c.n < 2) return std::nullopt;
  for (const auto& e : c.edges) {
    if (e.u == e.v || e.weight == 0) return std::nullopt;
  }
  std::vector<pv::Demand> q;
  for (Vertex x = 0; x < c.n; ++x) {
    for (Vertex y = x + 1; y < c.n; ++y) q.push_back({x, y});
  }
  sssp::paper::model::WordRng rng(5);
  const auto r = pv::estimate<W>(c.n, c.edges, q, rng);
  const auto d = sssp::test::floyd_warshall<W>(c.n, c.edges);
  for (std::size_t i = 0; i < q.size(); ++i) {
    if (r.index[i] == pv::kDisconnected) continue;
    const auto& z = r.table[r.index[i]];
    const unsigned __int128 zz = ((static_cast<unsigned __int128>(4) << z.scale) + 1) * z.anchor;
    if (zz > 2 * static_cast<unsigned __int128>(d[q[i].x][q[i].y])) {
      return "Z > 2 d_S at demand " + std::to_string(q[i].x) + "-" + std::to_string(q[i].y);
    }
  }
  return std::nullopt;
}

template <class Violation>
void check_driver(Violation violation, const std::vector<W>& weights, std::size_t max_edges) {
  const auto a = sh::run_seeded<W>(1, 200, random_case, violation, weights);
  INFO(a.report());
  REQUIRE(a.failing_seed.has_value());
  // (1) the first failing seed: every earlier seed passes.
  for (std::uint64_t s = 1; s < *a.failing_seed; ++s) CHECK_FALSE(violation(random_case(s)));
  // (2) the seed regenerates the original input.
  CHECK(random_case(*a.failing_seed) == a.original);
  // (3) deterministic shrinking.
  const auto b = sh::run_seeded<W>(1, 200, random_case, violation, weights);
  CHECK(b.failing_seed == a.failing_seed);
  CHECK(b.minimal == a.minimal);
  sh::Case<W> replay = random_case(*a.failing_seed);
  sh::shrink_case(replay, violation, weights);
  CHECK(replay == a.minimal);
  // (4) still failing and 1-minimal.
  CHECK(violation(a.minimal).has_value());
  CHECK(a.minimal.edges.size() <= max_edges);
  CHECK(a.minimal.edges.size() <= a.original.edges.size());
  for (std::size_t i = 0; i < a.minimal.edges.size(); ++i) {
    auto fewer = a.minimal;
    fewer.edges.erase(fewer.edges.begin() + static_cast<std::ptrdiff_t>(i));
    CHECK_FALSE(violation(fewer));
    for (const W w : weights) {
      if (!(w < a.minimal.edges[i].weight)) break;
      auto lighter = a.minimal;
      lighter.edges[i].weight = w;
      CHECK_FALSE(violation(lighter));
    }
  }
  // (5) the report names the seed and the minimal edge list.
  const std::string rep = a.report();
  CHECK(rep.find("seed " + std::to_string(*a.failing_seed)) != std::string::npos);
  CHECK(rep.find(sh::describe(a.minimal)) != std::string::npos);
}

}  // namespace

TEST_CASE("shrinking: a planted light cycle shrinks to a reproducible 1-minimal cycle",
          "[property][shrink]") {
  check_driver(light_cycle, {1, 2, 3, 4, 5, 6, 7, 8, 9}, 3);
  // The minimal light cycle is two parallel edges or a triangle of weight-1 edges on exactly
  // its own vertices.
  const auto a = sh::run_seeded<W>(1, 200, random_case, light_cycle, {1, 2, 3, 4, 5, 6, 7, 8, 9});
  CHECK((a.minimal.n == 2 || a.minimal.n == 3));
  CHECK(a.minimal.edges.size() == a.minimal.n);
  for (const auto& e : a.minimal.edges) CHECK(e.weight == 1);
  CHECK(a.minimal.b == 1);
}

TEST_CASE("shrinking: a false bound on Estimate shrinks to a reproducible minimal sample",
          "[property][shrink]") {
  check_driver(z_at_most_twice, {1, 2, 3, 4, 5, 6, 7, 8, 9}, 1);
  // Z >= 5a > 2 d_S already on one edge, so the minimum is a single edge on two pivots.
  const auto a = sh::run_seeded<W>(1, 200, random_case, z_at_most_twice, {1, 2, 3});
  CHECK(a.minimal.n == 2);
  REQUIRE(a.minimal.edges.size() == 1);
  CHECK(a.minimal.edges[0].weight == 1);
}

TEST_CASE("shrinking: no failure leaves the outcome empty", "[property][shrink]") {
  const auto a = sh::run_seeded<W>(1, 50, random_case,
                                   [](const sh::Case<W>&) { return std::optional<std::string>(); });
  CHECK_FALSE(a.failing_seed.has_value());
  CHECK(a.cases == 50);
  CHECK(a.report() == "no failure in 50 cases");
}
