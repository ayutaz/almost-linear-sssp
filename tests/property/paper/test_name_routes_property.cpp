// Definition 6.1, Sec 6.1, Algorithm 5 (NameRoutes) and Lemmas 6.2, 6.3, 6.8 of arXiv
// 2611.34911v1 on random integer instances over many seeds. Certificates are compared with exact
// references written independently of the library: Floyd-Warshall on G_P (tests/support) and
// the baseline Bellman-Ford oracle. Randomized claims (acceptance probability >= 1/2, Lemma 6.3)
// are measured, not assumed.
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::model::WordRng;
using sssp::paper::pivots::check_route_names;
using sssp::paper::pivots::name_routes;
using sssp::paper::pivots::name_routes_constants;
using sssp::paper::pivots::NameRoutesParams;
using sssp::paper::pivots::partial_shuffle;
using sssp::paper::pivots::PivotRoutingInstance;

namespace {

using W = std::uint64_t;
constexpr W kInf = WeightTraits<W>::infinity();

std::uint64_t counter(const sssp::Diagnostics& d, const std::string& k) {
  const auto it = d.counters.find(k);
  return it == d.counters.end() ? 0 : it->second;
}

// Loopless multigraph with parallel edges, weights over several decades, possibly disconnected.
PivotRoutingInstance<W> random_instance(std::uint64_t seed, std::size_t t_lo, std::size_t t_hi,
                                        double density) {
  sssp::test::SplitMix64 g(seed);
  PivotRoutingInstance<W> inst;
  inst.t = t_lo + g.below(t_hi - t_lo + 1);
  const auto pairs = static_cast<double>(inst.t * (inst.t - 1) / 2);
  const auto m = static_cast<std::size_t>(density * pairs) + g.below(inst.t + 1);
  for (std::size_t i = 0; i < m; ++i) {
    const auto u = static_cast<Vertex>(g.below(inst.t));
    auto v = static_cast<Vertex>(g.below(inst.t));
    if (u == v) v = static_cast<Vertex>((v + 1) % inst.t);
    W w = 1 + g.below(9);
    for (std::uint64_t k = g.below(4); k > 0; --k) w *= 23;
    inst.edges.push_back({u, v, w});
  }
  return inst;
}

PivotRoutingInstance<W> dense_instance(std::uint64_t seed, std::size_t t) {
  sssp::test::SplitMix64 g(seed);
  PivotRoutingInstance<W> inst;
  inst.t = t;
  for (Vertex x = 0; x < t; ++x) {
    for (Vertex y = x + 1; y < t; ++y) inst.edges.push_back({x, y, 1 + g.below(1000)});
  }
  // Some parallel edges, heavier or lighter than the first copy.
  for (std::size_t i = 0; i < t; ++i) {
    const auto u = static_cast<Vertex>(g.below(t));
    const auto v = static_cast<Vertex>((u + 1 + g.below(t - 1)) % t);
    inst.edges.push_back({v, u, 1 + g.below(2000)});
  }
  return inst;
}

// All-pairs distances of G_P by the baseline Bellman-Ford oracle (one run per source).
std::vector<std::vector<W>> bellman_ford_apsp(const PivotRoutingInstance<W>& inst) {
  const auto g = CsrGraph<W>::from_edges(inst.t, inst.edges);
  std::vector<std::vector<W>> d(inst.t);
  for (Vertex s = 0; s < inst.t; ++s) d[s] = sssp::baseline::bellman_ford(g, s).dist;
  return d;
}

}  // namespace

TEST_CASE("name_routes: certificates on random multigraphs, Mode B factors, many seeds",
          "[name_routes][property]") {
  // Small factors (sample t, residual 2t) make the sampling loop run on small t; correctness of
  // the certificates does not depend on the factors (Lemma 6.8(1)).
  NameRoutesParams p;
  p.sample_factor = 1;
  p.residual_factor = 2;
  for (std::uint64_t seed = 1; seed <= 120; ++seed) {
    const auto inst = random_instance(seed, 2, 14, 1.5);
    const auto dist = sssp::test::floyd_warshall<W>(inst.t, inst.edges);
    WordRng rng(seed * 104729);
    const auto r = name_routes(inst, rng, p);
    const auto v = check_route_names(inst, r, dist);
    INFO("seed " << seed << " t " << inst.t << " M " << inst.edges.size()
                 << (v.empty() ? "" : ": " + v.front()));
    REQUIRE(v.empty());
    CHECK(counter(r.diag, "name_routes.trials") ==
          counter(r.diag, "name_routes.accepted_rounds") +
              counter(r.diag, "name_routes.rejected_trials"));
    CHECK(r.dictionary.size() <= counter(r.diag, "name_routes.dictionary_records"));
    CHECK(r.diag.fallbacks.empty());
  }
}

TEST_CASE("name_routes: Mode A certificates against Bellman-Ford, measured acceptance",
          "[name_routes][property]") {
  std::uint64_t trials = 0;
  std::uint64_t accepted = 0;
  std::uint64_t worst_round = 0;
  for (std::uint64_t seed = 1; seed <= 6; ++seed) {
    const auto inst = dense_instance(seed, 66 + 2 * seed);
    const auto dist = bellman_ford_apsp(inst);
    WordRng rng(seed);
    const auto r = name_routes(inst, rng);
    const auto v = check_route_names(inst, r, dist);
    INFO("seed " << seed << (v.empty() ? "" : ": " + v.front()));
    REQUIRE(v.empty());
    REQUIRE(counter(r.diag, "name_routes.trials") >= 1);  // the sampling loop ran
    trials += counter(r.diag, "name_routes.trials");
    accepted += counter(r.diag, "name_routes.accepted_rounds");
    worst_round = std::max(worst_round, counter(r.diag, "name_routes.max_trials_per_round"));
    CHECK(counter(r.diag, "name_routes.terminal_edges") <= 32 * inst.t);
    CHECK(r.p_log2 == name_routes_constants(inst.t).p_log2);
    CHECK(r.diag.fallbacks.empty());
  }
  // Lemma 6.8(2): each trial is accepted with probability >= 1/2 (measured over all trials).
  INFO("accepted " << accepted << " of " << trials << " trials; worst round " << worst_round);
  CHECK(2 * accepted >= trials);
  CHECK(worst_round <= 20);
}

TEST_CASE("Lemma 6.2: at most 2t critical edges (brute force)", "[name_routes][property]") {
  // Critical: after deleting e = xy, no x-y path of <= g edges each no heavier than e.
  for (std::uint64_t seed = 1; seed <= 150; ++seed) {
    sssp::test::SplitMix64 r(seed);
    const std::size_t t = 2 + r.below(11);
    const std::size_t m = r.below(4 * t + 1);
    std::vector<Edge<W>> e;
    for (std::size_t i = 0; i < m; ++i) {
      const auto u = static_cast<Vertex>(r.below(t));
      auto v = static_cast<Vertex>(r.below(t));
      if (u == v) v = static_cast<Vertex>((v + 1) % t);
      e.push_back({u, v, 1 + r.below(6)});
    }
    const unsigned g = name_routes_constants(t).g;
    std::size_t critical = 0;
    for (std::size_t k = 0; k < m; ++k) {
      // Breadth-first layers from e[k].u using edges != k with weight <= w_k, up to g hops.
      std::vector<std::size_t> hop(t, SIZE_MAX);
      std::vector<Vertex> frontier{e[k].u};
      hop[e[k].u] = 0;
      for (unsigned d = 0; d < g && !frontier.empty(); ++d) {
        std::vector<Vertex> next;
        for (const Vertex x : frontier) {
          for (std::size_t i = 0; i < m; ++i) {
            if (i == k || e[i].weight > e[k].weight) continue;
            Vertex y = sssp::kNoVertex;
            if (e[i].u == x) y = e[i].v;
            if (e[i].v == x) y = e[i].u;
            if (y == sssp::kNoVertex || hop[y] != SIZE_MAX) continue;
            hop[y] = d + 1;
            next.push_back(y);
          }
        }
        frontier = std::move(next);
      }
      if (hop[e[k].v] == SIZE_MAX) ++critical;
    }
    INFO("seed " << seed << " t " << t << " m " << m);
    CHECK(critical <= 2 * t);
  }
}

TEST_CASE("Lemma 6.3: uniform 16t-samples leave few edges unwitnessed (measured)",
          "[name_routes][property]") {
  // Residual U with h > 32t edges, S a uniform 16t-subset: E[#{e : d_S(x,y) > g a_e}] <= h/8.
  const std::size_t t = 70;
  const unsigned g = name_routes_constants(t).g;
  double fraction_sum = 0;
  const int kSeeds = 12;
  for (int seed = 1; seed <= kSeeds; ++seed) {
    const auto inst = dense_instance(static_cast<std::uint64_t>(seed) + 500, t);
    // Simple residual: the complete graph's first t(t-1)/2 edges.
    std::vector<Edge<W>> u(inst.edges.begin(), inst.edges.begin() + t * (t - 1) / 2);
    REQUIRE(u.size() > 32 * t);
    WordRng rng(static_cast<std::uint64_t>(seed));
    std::vector<std::uint32_t> pos(u.size());
    for (std::size_t i = 0; i < pos.size(); ++i) pos[i] = static_cast<std::uint32_t>(i);
    partial_shuffle(std::span<std::uint32_t>(pos), 16 * t, rng);
    std::vector<Edge<W>> s;
    for (std::size_t i = 0; i < 16 * t; ++i) s.push_back(u[pos[i]]);
    const auto ds = sssp::test::floyd_warshall<W>(t, s);
    std::size_t bad = 0;
    for (const auto& e : u) {
      if (ds[e.u][e.v] == kInf || ds[e.u][e.v] > g * e.weight) ++bad;
    }
    fraction_sum += static_cast<double>(bad) / static_cast<double>(u.size());
  }
  const double mean = fraction_sum / kSeeds;
  INFO("mean unwitnessed fraction " << mean);
  CHECK(mean <= 0.125);
}
