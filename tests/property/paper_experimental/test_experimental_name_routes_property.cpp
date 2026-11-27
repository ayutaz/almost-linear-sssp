// Mode B NameRoutes (src/experimental/paper_experimental/name_routes.hpp), task
// P5-EXP-NAMEROUTES: Sec 6 claims under the deviations D-9 (H), D-10 (retry cap) and D-17
// (factors, P) on random instances and many seeds (arXiv 2611.34911v1).
//
// Oracles are independent of the library: Floyd-Warshall on G_P (tests/support) and the
// baseline Bellman-Ford. Checked claims:
//   * Lemma 6.8(1): every edge has a certificate (lambda_e <= a_e, dist_GP <= P lambda_e) on
//     every seed, for every parameter value (it does not depend on H, P, the factors or a cap).
//   * Lemma 6.8(2): a trial is accepted with probability >= 1/2. Its proof uses Lemma 6.7(2),
//     whose proof needs only H >= 4tT^2 (Q-42), so the claim is measured under D-9 with the
//     paper's P and factors, next to Mode A on the same instances.
//   * D-10: a cap hit is recorded once, with the residual size as its edge count, and is the
//     only way the residual can exceed the loop limit.
//   * Work (SSSP_DIAGNOSTICS only): counters against bounds from Alg 5 and the proof of
//     Lemma 6.8(4) (accepted rounds halve U; one record per qualifying group, at most one group
//     per Estimate table entry; terminal edges <= the loop limit unless capped).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "experimental/paper_experimental/name_routes.hpp"
#include "experimental/paper_experimental/params.hpp"
#include "paper/model/random.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

namespace pe = sssp::paper_experimental;
namespace pv = sssp::paper::pivots;
using sssp::CsrGraph;
using sssp::Vertex;
using sssp::paper::model::WordRng;

namespace {

using W = std::uint64_t;

std::uint64_t counter(const sssp::Diagnostics& d, const std::string& k) {
  const auto it = d.counters.find(k);
  return it == d.counters.end() ? 0 : it->second;
}

// Loopless multigraph with parallel edges, weights over several decades, possibly disconnected.
pv::PivotRoutingInstance<W> random_instance(std::uint64_t seed, std::size_t t_lo, std::size_t t_hi,
                                            double density) {
  sssp::test::SplitMix64 g(seed);
  pv::PivotRoutingInstance<W> inst;
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

pv::PivotRoutingInstance<W> dense_instance(std::uint64_t seed, std::size_t t) {
  sssp::test::SplitMix64 g(seed);
  pv::PivotRoutingInstance<W> inst;
  inst.t = t;
  for (Vertex x = 0; x < t; ++x) {
    for (Vertex y = x + 1; y < t; ++y) inst.edges.push_back({x, y, 1 + g.below(1000)});
  }
  for (std::size_t i = 0; i < t; ++i) {
    const auto u = static_cast<Vertex>(g.below(t));
    const auto v = static_cast<Vertex>((u + 1 + g.below(t - 1)) % t);
    inst.edges.push_back({v, u, 1 + g.below(2000)});
  }
  return inst;
}

std::vector<std::vector<W>> bellman_ford_apsp(const pv::PivotRoutingInstance<W>& inst) {
  const auto g = CsrGraph<W>::from_edges(inst.t, inst.edges);
  std::vector<std::vector<W>> d(inst.t);
  for (Vertex s = 0; s < inst.t; ++s) d[s] = sssp::baseline::bellman_ford(g, s).dist;
  return d;
}

/// Distinct endpoint pairs (independent count of Alg 5 l.4 representatives).
std::uint64_t distinct_pairs(const pv::PivotRoutingInstance<W>& inst) {
  std::vector<std::uint64_t> keys;
  for (const auto& e : inst.edges) {
    keys.push_back(std::uint64_t{std::min(e.u, e.v)} * inst.t + std::max(e.u, e.v));
  }
  std::sort(keys.begin(), keys.end());
  return static_cast<std::uint64_t>(std::unique(keys.begin(), keys.end()) - keys.begin());
}

/// Checks the D-10 bookkeeping and (with diagnostics) the work bounds for one run.
void check_run(const pv::PivotRoutingInstance<W>& inst, const pe::ExperimentalRouteNames<W>& r) {
  const auto& d = r.names.diag;
  const std::uint64_t reps = distinct_pairs(inst);
  const std::uint64_t limit = r.setting.residual_limit;
  // Counters written unconditionally by the shared component.
  const std::uint64_t trials = counter(d, "name_routes.trials");
  const std::uint64_t accepted = counter(d, "name_routes.accepted_rounds");
  const std::uint64_t rejected = counter(d, "name_routes.rejected_trials");
  const std::uint64_t terminal = counter(d, "name_routes.terminal_edges");
  CHECK(trials == accepted + rejected);
  CHECK(r.names.p_log2 == r.setting.p_log2);
  CHECK(terminal <= reps);  // Alg 5 l.15: terminal edges are residual representatives
  if (r.retry_cap_reached) {
    REQUIRE(r.setting.retry_cap.has_value());
    REQUIRE(d.fallbacks.size() == 1);
    CHECK(d.fallbacks[0].vertices == inst.t);
    CHECK(d.fallbacks[0].edges == terminal);
    CHECK(r.fallback_edges == terminal);
    CHECK(terminal > limit);  // the cap stopped the loop with |U| above the limit
    CHECK(counter(d, "name_routes.retry_cap_hits") == 1);
    CHECK(counter(d, "name_routes.max_consecutive_rejections") == *r.setting.retry_cap);
  } else {
    CHECK(d.fallbacks.empty());
    CHECK(r.fallback_edges == 0);
    CHECK(terminal <= limit);  // Alg 5 l.5: the loop ends only once |U| <= limit
    CHECK(counter(d, "name_routes.retry_cap_hits") == 0);
  }
  if (r.setting.retry_cap.has_value()) {
    CHECK(counter(d, "name_routes.max_consecutive_rejections") <= *r.setting.retry_cap);
  }
#if SSSP_DIAGNOSTICS
  CHECK(counter(d, "name_routes.representatives") == reps);
  // Alg 5 l.13-14: an accepted round removes at least half of U, and a round runs only while
  // |U| > limit, so the k-th accepted round saw limit < |U| <= reps / 2^(k-1).
  if (accepted >= 1) {
    REQUIRE(accepted < 64);
    CHECK((std::uint64_t{1} << (accepted - 1)) * (limit + 1) <= reps);
  }
  // Proof of Lemma 6.8(4): one record per qualifying group of an accepted round (a group is a
  // table index of that trial's Estimate) plus one per terminal edge.
  const std::uint64_t records = counter(d, "name_routes.dictionary_records");
  CHECK(records <= counter(d, "estimate.table_entries") + terminal);
  CHECK(counter(d, "name_routes.dictionary_size") == r.names.dictionary.size());
  CHECK(r.names.dictionary.size() <= records);
  CHECK(r.names.dictionary.size() <= reps);
  // Q-51 (b), Sec 2 l.303-304: one logarithm table over 1..|U| (one build step per index) for
  // all sample draws; none if no trial runs.
  CHECK(counter(d, "name_routes.ilog_table_entries") == (trials == 0 ? 0 : reps));
#endif
}

pe::PaperParams with_outer(const pe::OuterParams& outer) {
  return pe::PaperParams({}, outer);
}

}  // namespace

TEST_CASE("Mode B NameRoutes: certificates under D-9 + D-17 + D-10 on many seeds",
          "[paper_experimental][name_routes][property]") {
  // Small factors make the loop run at small t; small P makes trials fail so that the cap is
  // reached on some seeds. Lemma 6.8(1) must hold on every run either way.
  const unsigned p_values[] = {0, 1, 3, 6, 12};
  std::uint64_t runs = 0;
  std::uint64_t capped = 0;
  std::uint64_t looped = 0;
  for (std::uint64_t seed = 1; seed <= 150; ++seed) {
    const auto inst = random_instance(seed, 2, 16, 1.5);
    const auto dist = sssp::test::floyd_warshall<W>(inst.t, inst.edges);
    pe::OuterParams outer;
    outer.h_rule = pe::HRule::kProofMinimum;
    outer.sample_factor = 1;
    outer.residual_factor = 2;
    outer.p_log2 = p_values[seed % 5];
    outer.retry_cap = 1 + seed % 6;
    WordRng rng(seed * 104729);
    const auto r = pe::name_routes(inst, rng, with_outer(outer));
    const auto v = pv::check_route_names(inst, r.names, dist);
    INFO("seed " << seed << " t " << inst.t << " M " << inst.edges.size() << " P 2^"
                 << r.setting.p_log2 << (v.empty() ? "" : ": " + v.front()));
    REQUIRE(v.empty());
    check_run(inst, r);
    ++runs;
    capped += r.retry_cap_reached ? 1 : 0;
    looped += counter(r.names.diag, "name_routes.trials") > 0 ? 1 : 0;
  }
  // Both D-10 outcomes were exercised.
  INFO("runs " << runs << " with trials " << looped << " capped " << capped);
  CHECK(capped >= 10);
  CHECK(looped - capped >= 10);
}

TEST_CASE("Mode B NameRoutes: D-9 alone keeps acceptance >= 1/2 (measured beside Mode A)",
          "[paper_experimental][name_routes][property]") {
  pe::OuterParams d9;
  d9.h_rule = pe::HRule::kProofMinimum;
  const auto params_d9 = with_outer(d9);
  const pe::PaperParams paper;
  struct Totals {
    std::uint64_t trials = 0;
    std::uint64_t accepted = 0;
    std::uint64_t worst_round = 0;
    std::uint64_t dictionary = 0;
  };
  Totals mode_a;
  Totals mode_b;
  const auto add = [](Totals& to, const pe::ExperimentalRouteNames<W>& r) {
    to.trials += counter(r.names.diag, "name_routes.trials");
    to.accepted += counter(r.names.diag, "name_routes.accepted_rounds");
    to.worst_round =
        std::max(to.worst_round, counter(r.names.diag, "name_routes.max_trials_per_round"));
    to.dictionary += r.names.dictionary.size();
  };
  for (std::uint64_t seed = 1; seed <= 6; ++seed) {
    const auto inst = dense_instance(seed, 66 + 2 * seed);
    const auto dist = bellman_ford_apsp(inst);
    for (std::uint64_t rep = 0; rep < 2; ++rep) {
      WordRng rb(seed * 31 + rep);
      const auto b = pe::name_routes(inst, rb, params_d9);
      WordRng ra(seed * 31 + rep);
      const auto a = pe::name_routes(inst, ra, paper);
      INFO("seed " << seed << " rep " << rep);
      REQUIRE(pv::check_route_names(inst, b.names, dist).empty());
      REQUIRE(pv::check_route_names(inst, a.names, dist).empty());
      REQUIRE(counter(b.names.diag, "name_routes.trials") >= 1);
      CHECK(b.setting.h_log2 < b.setting.paper_h_log2);
      CHECK(a.setting.h_log2 == a.setting.paper_h_log2);
      check_run(inst, b);
      check_run(inst, a);
      CHECK_FALSE(b.retry_cap_reached);
      add(mode_b, b);
      add(mode_a, a);
    }
  }
  INFO("Mode A: accepted " << mode_a.accepted << " of " << mode_a.trials << ", worst round "
                           << mode_a.worst_round << ", dictionary total " << mode_a.dictionary);
  INFO("D-9:    accepted " << mode_b.accepted << " of " << mode_b.trials << ", worst round "
                           << mode_b.worst_round << ", dictionary total " << mode_b.dictionary);
  // Lemma 6.8(2) (with Q-42): each trial is accepted with probability >= 1/2.
  CHECK(2 * mode_b.accepted >= mode_b.trials);
  CHECK(2 * mode_a.accepted >= mode_a.trials);
  CHECK(mode_b.worst_round <= 20);
}

TEST_CASE("Mode B NameRoutes: double weights, certificates under D-9 + D-17",
          "[paper_experimental][name_routes][property]") {
  // Dyadic weights (k / 4): sums and doublings are exact, so the certificate test is exact.
  for (std::uint64_t seed = 1; seed <= 40; ++seed) {
    const auto ints = random_instance(seed + 900, 3, 14, 1.2);
    pv::PivotRoutingInstance<double> inst;
    inst.t = ints.t;
    for (const auto& e : ints.edges) {
      inst.edges.push_back({e.u, e.v, static_cast<double>(e.weight) / 4.0});
    }
    const auto dist = sssp::test::floyd_warshall<double>(inst.t, inst.edges);
    pe::OuterParams outer;
    outer.h_rule = pe::HRule::kProofMinimum;
    outer.sample_factor = 1;
    outer.residual_factor = 2;
    outer.retry_cap = 8;
    WordRng rng(seed);
    const auto r = pe::name_routes(inst, rng, with_outer(outer));
    const auto v = pv::check_route_names(inst, r.names, dist);
    INFO("seed " << seed << (v.empty() ? "" : ": " + v.front()));
    REQUIRE(v.empty());
    CHECK(r.names.diag.fallbacks.size() == (r.retry_cap_reached ? 1U : 0U));
  }
}
