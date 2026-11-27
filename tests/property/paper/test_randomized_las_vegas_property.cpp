// Task P4-RANDOMIZED: Las Vegas checks of Sec 6 of arXiv 2611.34911v1 in Mode A (the paper's
// parameters): NameRoutes (Algorithm 5, Lemma 6.8) and Estimate (Algorithm 6, Lemmas 6.5-6.7)
// on many seeds per instance family.
//
//   * Lemma 6.8(1) / Lemma 6.7(1), deterministic: on EVERY seed each edge has a certificate
//     (lambda_e <= a_e, dist_GP(x, y) <= P lambda_e) and each connected demand has Z >= d_S.
//     Checked with the independent checkers of randomized_support.hpp against Floyd-Warshall,
//     which is cross-checked against the baseline Bellman-Ford.
//   * Lemma 6.8(2), randomized: a trial is accepted with probability >= 1/2. The acceptance rate
//     and the distribution of the trial count are measured and printed (WARN); the tests reject
//     only if the measurement contradicts the lemma at level kDelta = 1e-9 (derivations in
//     randomized_support.hpp). No distribution is assumed.
//   * Lemmas 6.5, 6.6, 6.7(2), randomized: separation frequency at one scale, adjacency scans
//     and E[Z / d_S], measured in the same way.
//
// Mode A runs the sampling loop only when more than 32t distinct pairs exist, i.e. for t >= 66
// (C(66, 2) = 2145 > 2112). The NameRoutes families therefore use t >= 66; the sweep at small t
// checks the path without trials.
//
// The [slow] cases (ctest label property-slow) run 1000 seeds per family and weight type; the
// uint64 NameRoutes sweeps at t >= 66 do so in optimised builds and run 125 seeds in unoptimised
// ones (randomized_support.hpp: kHeavySweepSeeds). The other cases run the same code on about
// 100 seeds per test.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/pivots/estimate.hpp"
#include "paper/pivots/name_routes.hpp"
#include "property/paper/randomized_support.hpp"

namespace rz = sssp::test::randomized;
namespace pv = sssp::paper::pivots;
using rz::Family;
using rz::SampleFamily;
using sssp::paper::model::WordRng;

namespace {

using rz::FamilySpec;
using rz::kLoopFamilies;

/// Mode A NameRoutes on `seeds` seeds over `instances` instances of one family.
template <class W>
rz::Tally sweep_mode_a(const FamilySpec& spec, std::uint64_t seeds, std::size_t instances,
                       std::uint64_t salt) {
  std::vector<rz::Instance<W>> inst;
  std::vector<std::vector<std::vector<W>>> dist;
  for (std::size_t i = 0; i < instances; ++i) {
    inst.push_back(rz::make_instance<W>(spec.family, spec.t, 1000 * salt + i));
    dist.push_back(rz::reference_distances<W>(spec.t, inst.back().edges));
  }
  rz::RunSetting setting;
  setting.residual_limit = 32 * spec.t;       // Alg 5 l.5
  setting.p_log2 = rz::paper_p_log2(spec.t);  // Alg 5 l.3
  rz::Tally tally;
  for (std::uint64_t seed = 1; seed <= seeds; ++seed) {
    const std::size_t i = static_cast<std::size_t>(seed % instances);
    const std::uint64_t reps = rz::distinct_pairs(inst[i]);
    REQUIRE(reps > 32 * spec.t);
    REQUIRE(reps <= 64 * spec.t);
    WordRng rng(rz::rng_seed(seed, salt));
    const auto names = pv::name_routes(inst[i], rng);
    INFO(rz::family_name(spec.family) << " " << rz::weight_name<W>() << " t " << spec.t
                                      << " instance " << i << " seed " << seed);
    rz::check_run(inst[i], dist[i], names, setting, tally);
    CHECK(names.diag.fallbacks.empty());  // Mode A has no retry cap
    // reps <= 64t: after one accepted round at most 32t pairs remain.
    CHECK(rz::counter(names.diag, "name_routes.accepted_rounds") == 1);
  }
  return tally;
}

/// The statistical assertions on a Mode A tally in which every run has exactly one round.
void check_single_round_tally(const rz::Tally& t) {
  REQUIRE(t.runs_with_trials == t.runs);
  // Lemma 6.8(2) over all trials, valid for any stopping rule (acceptance_test).
  const auto a = rz::acceptance_test(t.accepted, t.trials);
  CHECK_FALSE(a.rejects_half);
  // The first trial of each run alone: the runs are independent and each first trial is
  // accepted with probability >= 1/2, so the fraction of runs without a rejection is at least
  // 1/2 - hoeffding_radius(runs) except with probability kDelta. (With one round per run, "no
  // rejection" is "the first trial was accepted".)
  const double first = static_cast<double>(t.runs_without_rejection) / static_cast<double>(t.runs);
  CHECK(first >= 0.5 - rz::hoeffding_radius(t.runs));
}

template <class W>
void run_mode_a_families(std::uint64_t seeds_per_family, std::size_t instances,
                         std::uint64_t salt_base) {
  rz::Tally total;
  std::uint64_t salt = salt_base;
  for (const auto& spec : kLoopFamilies) {
    const auto t = sweep_mode_a<W>(spec, seeds_per_family, instances, ++salt);
    WARN(rz::describe(std::string("Mode A NameRoutes, ") + rz::family_name(spec.family) +
                          " t=" + std::to_string(spec.t) + ", " + rz::weight_name<W>(),
                      t));
    total.merge(t);
  }
  WARN(
      rz::describe(std::string("Mode A NameRoutes, all families, ") + rz::weight_name<W>(), total));
  check_single_round_tally(total);
}

template <class W>
void run_mode_a_family(const FamilySpec& spec, std::uint64_t salt) {
  const std::uint64_t seeds = rz::loop_sweep_seeds<W>();
  const auto t = sweep_mode_a<W>(spec, seeds, 4, salt);
  WARN(rz::describe(std::string("Mode A NameRoutes, ") + rz::family_name(spec.family) +
                        " t=" + std::to_string(spec.t) + ", " + rz::weight_name<W>() + ", " +
                        std::to_string(seeds) + " seeds",
                    t));
  check_single_round_tally(t);
}

template <class W>
void run_estimate_families(std::uint64_t seeds, std::uint64_t salt_base) {
  const auto paper_params = [](std::size_t) { return pv::EstimateParams{}; };
  const auto paper_eta = [](std::size_t t) { return rz::paper_h_log2(t); };
  std::uint64_t salt = salt_base;
  for (const auto family :
       {SampleFamily::kSparseDecades, SampleFamily::kSparseEqual, SampleFamily::kMultiEdge,
        SampleFamily::kComplete, SampleFamily::kDenseTies}) {
    const auto t = rz::sweep_estimate<W>(family, seeds, ++salt, paper_params, paper_eta);
    WARN(rz::describe(std::string("Mode A Estimate, ") + rz::sample_family_name(family) + ", " +
                          rz::weight_name<W>(),
                      t));
    rz::check_estimate_tally(t);
  }
}

void run_separation(std::uint64_t seeds, std::uint64_t salt_base) {
  const auto paper_params = [](std::size_t) { return pv::EstimateParams{}; };
  std::uint64_t salt = salt_base;
  for (const auto family :
       {SampleFamily::kSparseDecades, SampleFamily::kSparseEqual, SampleFamily::kComplete}) {
    for (const std::uint64_t sample_seed : {3U, 11U}) {
      const auto s = rz::make_sample<std::uint64_t>(family, sample_seed);
      const auto r = rz::separation_test(s, seeds, ++salt, paper_params);
      WARN(rz::describe(std::string("Mode A Lemma 6.5, ") + rz::sample_family_name(family) +
                            " sample " + std::to_string(sample_seed) + " t=" + std::to_string(s.t),
                        r));
      CHECK(r.violations == 0);
    }
  }
}

/// Mode A at small t: no more than 32t pairs exist, so no trial runs and every edge is certified
/// by the lightest weight of its pair (Alg 5 l.4, l.15).
template <class W>
void run_small_t(std::uint64_t seeds, std::uint64_t salt) {
  rz::Tally tally;
  for (std::uint64_t seed = 1; seed <= seeds; ++seed) {
    const Family family = seed % 3 == 0   ? Family::kMultiEdge
                          : seed % 3 == 1 ? Family::kSparseRandom
                                          : Family::kCompleteFewWeights;
    const std::size_t t = 2 + static_cast<std::size_t>(seed % 23);
    const auto inst = rz::make_instance<W>(family, t, seed + salt);
    const auto dist = rz::reference_distances<W>(t, inst.edges);
    rz::RunSetting setting;
    setting.residual_limit = 32 * t;
    setting.p_log2 = rz::paper_p_log2(t);
    WordRng rng(rz::rng_seed(seed, salt));
    const auto names = pv::name_routes(inst, rng);
    INFO(rz::family_name(family) << " t " << t << " seed " << seed);
    rz::check_run(inst, dist, names, setting, tally);
    // Each edge is named by the lightest weight among the edges of its pair.
    for (std::size_t i = 0; i < inst.edges.size(); ++i) {
      W lightest = inst.edges[i].weight;
      for (const auto& e : inst.edges) {
        const bool same = (e.u == inst.edges[i].u && e.v == inst.edges[i].v) ||
                          (e.u == inst.edges[i].v && e.v == inst.edges[i].u);
        if (same && e.weight < lightest) lightest = e.weight;
      }
      REQUIRE(names.dictionary.at(names.name.at(i)) == lightest);
    }
  }
  CHECK(tally.trials == 0);
  WARN(rz::describe(
      std::string("Mode A NameRoutes, t in 2..24 (no trials), ") + rz::weight_name<W>(), tally));
}

}  // namespace

// ----------------------------------------------------------------------- fast (every CI run)

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes: certificates and acceptance, uint64",
          "[randomized][name_routes][property]") {
  run_mode_a_families<std::uint64_t>(5, 2, 100);  // 7 families x 5 seeds
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes: certificates and acceptance, dyadic double",
          "[randomized][name_routes][property]") {
  run_mode_a_families<double>(15, 2, 200);  // 7 families x 15 seeds = 105 runs
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes: small t, no trials",
          "[randomized][name_routes][property]") {
  run_small_t<std::uint64_t>(100, 300);
  run_small_t<double>(100, 400);
}

TEST_CASE("P4-RANDOMIZED Mode A Estimate: Z >= d_S, scans and E[Z/d_S], uint64",
          "[randomized][estimate][property]") {
  run_estimate_families<std::uint64_t>(100, 500);
}

TEST_CASE("P4-RANDOMIZED Mode A Estimate: Z >= d_S, scans and E[Z/d_S], dyadic double",
          "[randomized][estimate][property]") {
  run_estimate_families<double>(100, 600);
}

TEST_CASE("P4-RANDOMIZED Mode A Estimate: separation at one scale (Lemma 6.5)",
          "[randomized][estimate][property]") {
  run_separation(100, 700);
}

// --------------------------------------------------- slow (seed sweeps, one family per case)

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: complete, uint64",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<std::uint64_t>(kLoopFamilies[0], 1100);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: complete, dyadic double",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<double>(kLoopFamilies[0], 1150);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: complete, equal weights, uint64",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<std::uint64_t>(kLoopFamilies[1], 1200);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: complete, equal weights, dyadic double",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<double>(kLoopFamilies[1], 1250);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: complete, few weights, uint64",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<std::uint64_t>(kLoopFamilies[2], 1300);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: complete, few weights, dyadic double",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<double>(kLoopFamilies[2], 1350);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: complete, wide weight range, uint64",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<std::uint64_t>(kLoopFamilies[3], 1400);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: complete, wide weight range, dyadic double",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<double>(kLoopFamilies[3], 1450);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: dense, uint64",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<std::uint64_t>(kLoopFamilies[4], 1500);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: dense, dyadic double",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<double>(kLoopFamilies[4], 1550);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: multi-edge, uint64",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<std::uint64_t>(kLoopFamilies[5], 1600);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: multi-edge, dyadic double",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<double>(kLoopFamilies[5], 1650);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: two cliques (disconnected), uint64",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<std::uint64_t>(kLoopFamilies[6], 1700);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: two cliques (disconnected), dyadic double",
          "[randomized][name_routes][property][slow]") {
  run_mode_a_family<double>(kLoopFamilies[6], 1750);
}

TEST_CASE("P4-RANDOMIZED Mode A NameRoutes, seed sweep: small t, no trials",
          "[randomized][name_routes][property][slow]") {
  run_small_t<std::uint64_t>(rz::kSweepSeeds, 1900);
  run_small_t<double>(rz::kSweepSeeds, 1950);
}

TEST_CASE("P4-RANDOMIZED Mode A Estimate, seed sweep, uint64",
          "[randomized][estimate][property][slow]") {
  run_estimate_families<std::uint64_t>(rz::kSweepSeeds, 2000);
}

TEST_CASE("P4-RANDOMIZED Mode A Estimate, seed sweep, dyadic double",
          "[randomized][estimate][property][slow]") {
  run_estimate_families<double>(rz::kSweepSeeds, 2100);
}

TEST_CASE("P4-RANDOMIZED Mode A Estimate, seed sweep: separation at one scale (Lemma 6.5)",
          "[randomized][estimate][property][slow]") {
  run_separation(rz::kSweepSeeds, 2200);
}
