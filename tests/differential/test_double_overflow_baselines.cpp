// Differential test of the double-weight overflow behaviour of every baseline
// (task P2-DIFF-DOUBLE-OVERFLOW, regression 20261127-double-sum-overflow).
//
// Policy (docs/design/baselines.md, numerics.md Sec 6): a result is never silently wrong. A
// finite distance plus a finite weight whose rounded sum is +inf (the "unreachable" sentinel)
// raises std::overflow_error, with the same split as for uint64:
//   * Bellman-Ford throws iff some reachable vertex has rounded distance d^(v) = +inf;
//   * every Dijkstra (binary, 4-ary, Fibonacci; the radix heap is uint64-only) throws iff some
//     reachable u with finite d^(u) has an arc of weight w with fl(d^(u) + w) = +inf, even if
//     that arc is not on a shortest path (e.g. a huge edge back to the source).
// Hence Bellman-Ford throwing implies every Dijkstra throws, never the converse.
//
// d^(v) = min over s-v paths of the left-to-right rounded path sum (numerics.md Sec 2) is
// computed below by an independent relaxation over the edge list with plain IEEE-754 addition,
// in which +inf saturates; reachability is computed separately (a reachable vertex can have
// d^(v) = +inf, which is exactly the case under test).

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/verify.hpp"
#include "support/prng.hpp"

namespace {

using sssp::CsrGraph;
using sssp::Edge;
using sssp::SsspResult;
using sssp::Vertex;

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kMax = std::numeric_limits<double>::max();

struct Algorithm {
  std::string name;
  std::function<SsspResult<double>(const CsrGraph<double>&, Vertex)> run;
};

std::vector<Algorithm> dijkstras() {
  return {
      {"dijkstra_binary_heap", [](const CsrGraph<double>& g,
                                  Vertex s) { return sssp::baseline::dijkstra_binary_heap(g, s); }},
      {"dijkstra_4ary_heap", [](const CsrGraph<double>& g,
                                Vertex s) { return sssp::baseline::dijkstra_4ary_heap(g, s); }},
      {"dijkstra_fibonacci_heap",
       [](const CsrGraph<double>& g, Vertex s) {
         return sssp::baseline::dijkstra_fibonacci_heap(g, s);
       }},
  };
}

bool throws_overflow(const std::function<SsspResult<double>()>& f) {
  try {
    static_cast<void>(f());
  } catch (const std::overflow_error&) {
    return true;
  }
  return false;
}

// ------------------------------------------------------------------ independent double oracle

struct Expectation {
  std::vector<double> rounded;  ///< d^(v), +inf if unreachable or not representable
  bool bellman_ford_throws = false;
  bool dijkstra_throws = false;
};

Expectation expect(std::size_t n, const std::vector<Edge<double>>& edges, Vertex s) {
  Expectation x;
  x.rounded.assign(n, kInf);
  std::vector<char> reachable(n, 0);
  x.rounded[s] = 0.0;
  reachable[s] = 1;
  for (bool changed = true; changed;) {
    changed = false;
    for (const auto& e : edges) {
      for (const auto [a, b] : {std::pair{e.u, e.v}, std::pair{e.v, e.u}}) {
        if (reachable[a] == 0) continue;
        if (reachable[b] == 0) {
          reachable[b] = 1;
          changed = true;
        }
        const double c = x.rounded[a] + e.weight;  // +inf saturates: inf + w = inf
        if (c < x.rounded[b]) {
          x.rounded[b] = c;
          changed = true;
        }
      }
    }
  }
  for (std::size_t v = 0; v < n; ++v) {
    if (reachable[v] != 0 && x.rounded[v] == kInf) x.bellman_ford_throws = true;
  }
  for (const auto& e : edges) {
    for (const Vertex u : {e.u, e.v}) {
      if (x.rounded[u] != kInf && x.rounded[u] + e.weight == kInf) x.dijkstra_throws = true;
    }
  }
  return x;
}

std::string describe(std::size_t n, const std::vector<Edge<double>>& edges, Vertex s) {
  std::ostringstream out;
  out.precision(17);  // %.17g: the exact double
  out << "n=" << n << " source=" << s << " edges:";
  for (const auto& e : edges) out << " {" << e.u << "," << e.v << "," << e.weight << "}";
  return out.str();
}

}  // namespace

TEST_CASE("double overflow: every baseline throws on an unrepresentable distance",
          "[differential][baselines][overflow]") {
  // d(2) = 2e308 > DBL_MAX (regression 20261127-double-sum-overflow).
  const auto g = CsrGraph<double>::from_edges(3, {{0, 1, 1e308}, {1, 2, 1e308}});
  CHECK_THROWS_AS(sssp::baseline::bellman_ford(g, 0), std::overflow_error);
  for (const auto& alg : dijkstras()) {
    INFO(alg.name);
    CHECK_THROWS_AS(alg.run(g, 0), std::overflow_error);
  }
}

TEST_CASE("double overflow: Dijkstra checks every scanned arc, Bellman-Ford only distances",
          "[differential][baselines][overflow]") {
  // d = {0, 1e308} is representable, but scanning vertex 1 offers 1e308 + 1e308 = +inf back to
  // the source. Same split as the uint64 policy (baselines.md "Integer weights and overflow").
  const auto g = CsrGraph<double>::from_edges(2, {{0, 1, 1e308}});
  const auto bf = sssp::baseline::bellman_ford(g, 0);
  CHECK(bf.dist == std::vector<double>{0.0, 1e308});
  for (const auto& alg : dijkstras()) {
    INFO(alg.name);
    CHECK_THROWS_AS(alg.run(g, 0), std::overflow_error);
  }
}

TEST_CASE("double overflow: large but finite sums are not reported",
          "[differential][baselines][overflow]") {
  // Every offered sum is at most 3/4 DBL_MAX: no baseline may throw, distances are exact.
  const double q = kMax / 4;  // exact: division by a power of two
  const auto g = CsrGraph<double>::from_edges(3, {{0, 1, q}, {1, 2, q}});
  const std::vector<double> want{0.0, q, q + q};
  CHECK(sssp::baseline::bellman_ford(g, 0).dist == want);
  for (const auto& alg : dijkstras()) {
    INFO(alg.name);
    CHECK(alg.run(g, 0).dist == want);
  }
}

TEST_CASE("double overflow: baselines match the independent rounded oracle on huge weights",
          "[differential][baselines][overflow]") {
  // Weights clustered near the top of the double range so that both outcomes occur often;
  // ulp(DBL_MAX) = 2^971: DBL_MAX + 2^969 is absorbed (= DBL_MAX), DBL_MAX + 2^970 rounds to +inf.
  const std::vector<double> pool{0.0,      1.0,      0x1p969,  0x1p970,  kMax / 8,   kMax / 4,
                                 kMax / 3, kMax / 2, 0x1p1022, 0x1p1023, 0x1.8p1023, kMax};
  std::uint64_t instances = 0;
  std::uint64_t bf_throws = 0;
  std::uint64_t dijkstra_only_throws = 0;
  std::uint64_t clean = 0;
  for (const std::uint64_t seed : {1ULL, 7ULL, 42ULL, 2026ULL, 0xC0FFEEULL}) {
    sssp::test::SplitMix64 rng(seed);
    for (int round = 0; round < 300; ++round) {
      const auto n = static_cast<std::size_t>(rng.between(1, 7));
      const auto m = static_cast<std::size_t>(rng.between(0, 10));
      std::vector<Edge<double>> edges;
      for (std::size_t i = 0; i < m; ++i) {
        const auto u = static_cast<Vertex>(rng.below(n));
        const auto v = static_cast<Vertex>(rng.below(n));
        const double w = rng.below(4) == 0 ? rng.unit() * kMax : pool[rng.below(pool.size())];
        edges.push_back({u, v, w});
      }
      const auto g = CsrGraph<double>::from_edges(n, edges);
      for (Vertex s = 0; s < n; ++s) {
        ++instances;
        const auto ex = expect(n, edges, s);
        const std::string where = describe(n, edges, s);
        REQUIRE((!ex.bellman_ford_throws || ex.dijkstra_throws));  // oracle self-consistency

        const bool bf_threw = throws_overflow([&] { return sssp::baseline::bellman_ford(g, s); });
        if (bf_threw != ex.bellman_ford_throws) {
          FAIL("bellman_ford threw=" << bf_threw << " expected=" << ex.bellman_ford_throws << "\n"
                                     << where);
        }
        if (!bf_threw) {
          const auto bf = sssp::baseline::bellman_ford(g, s);
          if (bf.dist != ex.rounded) FAIL("bellman_ford dist differs from the oracle\n" << where);
        }
        for (const auto& alg : dijkstras()) {
          SsspResult<double> r;
          bool threw = false;
          try {
            r = alg.run(g, s);
          } catch (const std::overflow_error&) {
            threw = true;
          }
          if (threw != ex.dijkstra_throws) {
            FAIL(alg.name << " threw=" << threw << " expected=" << ex.dijkstra_throws << "\n"
                          << where);
          }
          if (threw) continue;
          if (r.dist != ex.rounded) FAIL(alg.name << " dist differs from the oracle\n" << where);
          const auto rep = sssp::verify_sssp(g, s, r);
          if (!rep.ok) FAIL(alg.name << " fails verify_sssp\n" << where << "\n" << rep.summary());
        }
        if (ex.bellman_ford_throws) {
          ++bf_throws;
        } else if (ex.dijkstra_throws) {
          ++dijkstra_only_throws;
        } else {
          ++clean;
        }
      }
    }
  }
  // The sample must exercise all three outcomes, or the test proves nothing.
  CAPTURE(instances, bf_throws, dijkstra_only_throws, clean);
  CHECK(bf_throws > 100);
  CHECK(dijkstra_only_throws > 100);
  CHECK(clean > 100);
}
