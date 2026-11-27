// Definition 6.1, Sec 6.1, Algorithm 5 (NameRoutes) of arXiv 2611.34911v1 on hand-computed
// examples. Constants are derived from Alg 5 l.2-3; certificates are checked against the
// Floyd-Warshall reference on G_P (tests/support).
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/graph.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

using sssp::Edge;
using sssp::Vertex;
using sssp::paper::model::WordRng;
using sssp::paper::pivots::check_route_names;
using sssp::paper::pivots::median_prune;
using sssp::paper::pivots::name_routes;
using sssp::paper::pivots::name_routes_constants;
using sssp::paper::pivots::NameRoutesParams;
using sssp::paper::pivots::partial_shuffle;
using sssp::paper::pivots::PivotRoutingInstance;

namespace {

using W = std::uint64_t;

std::vector<std::vector<W>> gp_dist(const PivotRoutingInstance<W>& g) {
  return sssp::test::floyd_warshall<W>(g.t, g.edges);
}

PivotRoutingInstance<W> complete(std::size_t t, std::uint64_t lo, std::uint64_t hi,
                                 std::uint64_t seed) {
  sssp::test::SplitMix64 r(seed);
  PivotRoutingInstance<W> g;
  g.t = t;
  for (Vertex x = 0; x < t; ++x) {
    for (Vertex y = x + 1; y < t; ++y) g.edges.push_back({x, y, r.between(lo, hi)});
  }
  return g;
}

std::uint64_t counter(const sssp::Diagnostics& d, const std::string& k) {
  const auto it = d.counters.find(k);
  return it == d.counters.end() ? 0 : it->second;
}

}  // namespace

TEST_CASE("name_routes constants follow Alg 5 l.2-3", "[name_routes]") {
  // t = 3: u = ceil(log2 5) = 3, g = 15, 2048 * 15 * 9 = 276480 in (2^18, 2^19].
  auto c = name_routes_constants(3);
  CHECK(c.u == 3);
  CHECK(c.g == 15);
  CHECK(c.p_log2 == 19);
  // t = 80: u = ceil(log2 82) = 7, g = 31, 2048 * 31 * 49 = 3110912 in (2^21, 2^22].
  c = name_routes_constants(80);
  CHECK(c.u == 7);
  CHECK(c.g == 31);
  CHECK(c.p_log2 == 22);
  // t = 2: u = 2, g = 11, 2048 * 11 * 4 = 90112 in (2^16, 2^17].
  c = name_routes_constants(2);
  CHECK(c.u == 2);
  CHECK(c.g == 11);
  CHECK(c.p_log2 == 17);
}

TEST_CASE("median_prune finds the least feasible key (Sec 6.1)", "[name_routes]") {
  const std::vector<W> keys{5, 1, 9, 3, 7};
  std::uint64_t tests = 0;
  auto r = median_prune<W>(keys, [](W k) { return k >= 4; }, nullptr, &tests);
  REQUIRE(r.has_value());
  CHECK(*r == 5);
  CHECK(tests <= 3);  // each test halves the remaining keys: <= floor(log2 5) + 1
  CHECK_FALSE(median_prune<W>(keys, [](W k) { return k >= 100; }).has_value());
  r = median_prune<W>(keys, [](W) { return true; });
  REQUIRE(r.has_value());
  CHECK(*r == 1);
  // Ties are discarded on both sides of the median.
  const std::vector<W> ties{2, 2, 2, 2};
  tests = 0;
  r = median_prune<W>(ties, [](W k) { return k >= 2; }, nullptr, &tests);
  REQUIRE(r.has_value());
  CHECK(*r == 2);
  CHECK(tests == 1);
  CHECK_FALSE(median_prune<W>(std::span<const W>{}, [](W) { return true; }).has_value());
  // Against brute force on random keys.
  sssp::test::SplitMix64 g(99);
  for (int rep = 0; rep < 50; ++rep) {
    std::vector<W> k(1 + g.below(1000));
    for (auto& x : k) x = g.below(500);
    const W thr = g.below(520);
    W expect = std::numeric_limits<W>::max();
    for (const W x : k) {
      if (x >= thr) expect = std::min(expect, x);
    }
    tests = 0;
    const auto got = median_prune<W>(k, [&](W x) { return x >= thr; }, nullptr, &tests);
    if (expect == std::numeric_limits<W>::max()) {
      CHECK_FALSE(got.has_value());
    } else {
      REQUIRE(got.has_value());
      CHECK(*got == expect);
    }
    std::uint64_t bound = 1;
    while ((std::uint64_t{1} << bound) <= k.size()) ++bound;  // floor(log2 n) + 1
    CHECK(tests <= bound);
  }
}

TEST_CASE("partial_shuffle draws a uniform k-subset (Alg 5 l.6)", "[name_routes]") {
  WordRng rng(2024);
  std::map<std::pair<int, int>, int> count;
  constexpr int kTrials = 20000;
  for (int rep = 0; rep < kTrials; ++rep) {
    std::vector<int> a{0, 1, 2, 3, 4};
    partial_shuffle(std::span<int>(a), 2, rng);
    std::vector<int> s = a;
    std::sort(s.begin(), s.end());
    REQUIRE(s == std::vector<int>({0, 1, 2, 3, 4}));  // a permutation of the input
    ++count[{std::min(a[0], a[1]), std::max(a[0], a[1])}];
  }
  REQUIRE(count.size() == 10);
  for (const auto& [pair, c] : count) {
    // Expected 2000 per subset; the band is about 7 standard deviations wide.
    CHECK(c > 1700);
    CHECK(c < 2300);
  }
  std::vector<int> b{1, 2};
  CHECK_THROWS_AS(partial_shuffle(std::span<int>(b), 3, rng), std::invalid_argument);
}

TEST_CASE("name_routes rejects inputs outside Definition 6.1", "[name_routes]") {
  WordRng rng(1);
  PivotRoutingInstance<W> loop{3, {{1, 1, 4}}};
  CHECK_THROWS_AS(name_routes(loop, rng), std::invalid_argument);
  PivotRoutingInstance<W> zero{3, {{0, 1, 0}}};
  CHECK_THROWS_AS(name_routes(zero, rng), std::invalid_argument);
  PivotRoutingInstance<W> out{3, {{0, 3, 1}}};
  CHECK_THROWS_AS(name_routes(out, rng), std::out_of_range);
  PivotRoutingInstance<W> ok{3, {{0, 1, 1}}};
  NameRoutesParams bad;
  bad.sample_factor = 4;
  bad.residual_factor = 2;  // a sample must fit into the residual
  CHECK_THROWS_AS(name_routes(ok, rng, bad), std::invalid_argument);
  NameRoutesParams cap0;
  cap0.retry_cap = 0;
  CHECK_THROWS_AS(name_routes(ok, rng, cap0), std::invalid_argument);
}

TEST_CASE("name_routes with one pivot or without pivot edges is empty", "[name_routes]") {
  WordRng rng(3);
  // t = 1: the loopless instance has no edges (l.1422).
  const PivotRoutingInstance<W> one{1, {}};
  auto r = name_routes(one, rng);
  CHECK(r.dictionary.empty());
  CHECK(r.name.empty());
  CHECK(r.p_log2 == name_routes_constants(1).p_log2);
  // Every physical edge radial: no pivot-routing edge at all.
  const PivotRoutingInstance<W> none{5, {}};
  r = name_routes(none, rng);
  CHECK(r.dictionary.empty());
  CHECK(r.name.empty());
  CHECK(counter(r.diag, "name_routes.trials") == 0);
}

TEST_CASE("name_routes keeps the lightest parallel edge and coalesces ties", "[name_routes]") {
  // |U| = 2 <= 32t: no sampling; every representative certifies its own weight (l.15) and the
  // parallel edges reuse their representative's certificate (Sec 6.1).
  const PivotRoutingInstance<W> g{3, {{0, 1, 5}, {1, 0, 3}, {1, 2, 7}, {2, 1, 7}}};
  WordRng rng(5);
  const auto r = name_routes(g, rng);
  CHECK(r.dictionary == std::vector<W>({3, 7}));
  CHECK(r.name == std::vector<std::uint32_t>({0, 0, 1, 1}));
  CHECK(counter(r.diag, "name_routes.representatives") == 2);
  CHECK(counter(r.diag, "name_routes.trials") == 0);
  CHECK(counter(r.diag, "name_routes.terminal_edges") == 2);
  CHECK(check_route_names(g, r, gp_dist(g)).empty());
  // Equal weights on different pairs share one dictionary entry (l.16).
  const PivotRoutingInstance<W> h{4, {{0, 1, 4}, {2, 3, 4}, {1, 2, 9}}};
  const auto s = name_routes(h, rng);
  CHECK(s.dictionary == std::vector<W>({4, 9}));
  CHECK(s.name == std::vector<std::uint32_t>({0, 0, 1}));
}

TEST_CASE("name_routes samples with Mode A constants on a dense instance", "[name_routes]") {
  // t = 70, complete: 2415 representatives > 32t = 2240, so Alg 5 l.5-14 runs.
  const auto g = complete(70, 1, 1000, 77);
  const auto dist = gp_dist(g);
  WordRng rng(123);
  const auto r = name_routes(g, rng);
  const auto v = check_route_names(g, r, dist);
  INFO((v.empty() ? std::string() : v.front()));
  CHECK(v.empty());
  const auto trials = counter(r.diag, "name_routes.trials");
  CHECK(trials >= 1);
  CHECK(counter(r.diag, "name_routes.accepted_rounds") >= 1);
  CHECK(trials == counter(r.diag, "name_routes.accepted_rounds") +
                      counter(r.diag, "name_routes.rejected_trials"));
  CHECK(counter(r.diag, "estimate.calls") == trials);
  CHECK(counter(r.diag, "name_routes.max_trials_per_round") >= 1);
  CHECK(counter(r.diag, "name_routes.max_trials_per_round") <= 20);
  CHECK(counter(r.diag, "name_routes.terminal_edges") <= 32 * 70);
  CHECK(r.p_log2 == name_routes_constants(70).p_log2);
  CHECK(r.diag.fallbacks.empty());
}

TEST_CASE("name_routes with Mode B factors runs several rounds and stays valid", "[name_routes]") {
  // sample 1t, residual 2t (Mode B): the loop runs on t = 9, complete (36 > 18).
  NameRoutesParams p;
  p.sample_factor = 1;
  p.residual_factor = 2;
  for (std::uint64_t seed = 1; seed <= 10; ++seed) {
    const auto g = complete(9, 1, 50, seed);
    WordRng rng(seed);
    const auto r = name_routes(g, rng, p);
    CHECK(check_route_names(g, r, gp_dist(g)).empty());
    CHECK(counter(r.diag, "name_routes.trials") >= 1);
  }
}

TEST_CASE("name_routes retry cap (deviation D-10) is off by default and recorded when hit",
          "[name_routes]") {
  // P = 1 (Mode B) and unit weights: every bound is Z >= 5 > P * 1, so no edge ever qualifies
  // and every trial is rejected. With a cap of 1 the loop ends after one rejected trial and the
  // residual certifies its own weight.
  const auto g = complete(8, 1, 1, 0);
  NameRoutesParams p;
  p.sample_factor = 1;
  p.residual_factor = 2;
  p.p_log2 = 0;
  p.retry_cap = 1;
  WordRng rng(8);
  const auto r = name_routes(g, rng, p);
  CHECK(r.p_log2 == 0);
  CHECK(r.dictionary == std::vector<W>({1}));
  CHECK(check_route_names(g, r, gp_dist(g)).empty());
  CHECK(counter(r.diag, "name_routes.retry_cap_hits") == 1);
  CHECK(counter(r.diag, "name_routes.rejected_trials") == 1);
  CHECK(counter(r.diag, "name_routes.accepted_rounds") == 0);
  REQUIRE(r.diag.fallbacks.size() == 1);
  CHECK(r.diag.fallbacks[0].reason.find("D-10") != std::string::npos);
  CHECK(r.diag.fallbacks[0].edges == 28);
  // A cap that is never reached leaves no record.
  NameRoutesParams q;
  q.sample_factor = 1;
  q.residual_factor = 2;
  q.retry_cap = 64;
  const auto h = complete(9, 1, 50, 3);
  const auto s = name_routes(h, rng, q);
  CHECK(counter(s.diag, "name_routes.retry_cap_hits") == 0);
  CHECK(s.diag.fallbacks.empty());
}

TEST_CASE("name_routes is deterministic for a fixed seed", "[name_routes]") {
  const auto g = complete(70, 1, 1000, 5);
  WordRng a(9);
  WordRng b(9);
  const auto ra = name_routes(g, a);
  const auto rb = name_routes(g, b);
  CHECK(ra.dictionary == rb.dictionary);
  CHECK(ra.name == rb.name);
  CHECK(ra.diag.counters == rb.diag.counters);
}

TEST_CASE("check_route_names reports invalid certificates", "[name_routes]") {
  const PivotRoutingInstance<W> g{3, {{0, 1, 5}, {1, 2, 7}, {0, 2, 100}}};
  const auto dist = gp_dist(g);
  WordRng rng(1);
  auto r = name_routes(g, rng);
  REQUIRE(check_route_names(g, r, dist).empty());
  // lambda > a_e.
  auto bad = r;
  bad.name[0] = bad.name[2];  // edge 0 (a = 5) named 100
  CHECK_FALSE(check_route_names(g, bad, dist).empty());
  // dist > P lambda: P = 1, edge 2 (a = 100, dist 12) named 5.
  auto bad2 = r;
  bad2.p_log2 = 0;
  bad2.name[2] = bad2.name[0];
  bad2.name[1] = bad2.name[0];  // 5 <= 7 but dist(1,2) = 7 > 5
  CHECK_FALSE(check_route_names(g, bad2, dist).empty());
  // Unsorted dictionary or a name outside it.
  auto bad3 = r;
  std::reverse(bad3.dictionary.begin(), bad3.dictionary.end());
  CHECK_FALSE(check_route_names(g, bad3, dist).empty());
  auto bad4 = r;
  bad4.name[1] = 99;
  CHECK_FALSE(check_route_names(g, bad4, dist).empty());
}
