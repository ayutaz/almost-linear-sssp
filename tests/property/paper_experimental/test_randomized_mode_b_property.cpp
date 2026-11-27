// Task P4-RANDOMIZED: Las Vegas checks of Sec 6 of arXiv 2611.34911v1 under the Mode B
// deviations of NameRoutes and Estimate, on many seeds per instance family:
//   D-9   H = least power of two >= 4tT^2 instead of >= (16t)^5 (Sec 6.2 l.1491, Q-42);
//   D-10  retry cap (not in the paper, Q-45): the number of runs that reach it is counted;
//   D-17  sample factor, residual factor and P other than 16t, 32t and the P of Alg 5 l.3.
//
// What is asserted and what is only measured:
//   * Lemma 6.8(1) / Lemma 6.7(1) hold for every parameter value, so certificates and Z >= d_S
//     are REQUIREd on every seed (independent checkers of randomized_support.hpp, Floyd-Warshall
//     cross-checked against the baseline Bellman-Ford).
//   * D-9 alone keeps the proof of Lemma 6.8(2) (it needs only H >= 4tT^2, Q-42), so the
//     acceptance rate is tested against 1/2 exactly as in Mode A.
//   * Under D-17, Lemma 6.8(2) does not apply. The acceptance rate, the trial-count distribution
//     and the cap hits are measured and printed per value of P; no rate is asserted.
//   * Q-64: the Estimate scan-count bound of Lemma 6.6 and the separation bound of Lemma 6.5 are
//     checked under D-9, including the configuration of the P3-ESTIMATE test.
// Statistical thresholds: randomized_support.hpp (false-positive probability <= 1e-9 each).

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "experimental/paper_experimental/name_routes.hpp"
#include "experimental/paper_experimental/params.hpp"
#include "paper/model/random.hpp"
#include "paper/pivots/estimate.hpp"
#include "paper/pivots/name_routes.hpp"
#include "property/paper/randomized_support.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

namespace rz = sssp::test::randomized;
namespace pe = sssp::paper_experimental;
namespace pv = sssp::paper::pivots;
using rz::Family;
using rz::FamilySpec;
using rz::kLoopFamilies;
using rz::SampleFamily;
using sssp::paper::model::WordRng;

namespace {

pe::PaperParams with_outer(const pe::OuterParams& outer) {
  return pe::PaperParams({}, outer);
}

pe::PaperParams d9_params() {
  pe::OuterParams outer;
  outer.h_rule = pe::HRule::kProofMinimum;
  return with_outer(outer);
}

/// Checks the reported setting against values computed in the test and the run itself.
template <class W>
void check_mode_b_run(const rz::Instance<W>& inst, const std::vector<std::vector<W>>& dist,
                      const pe::ExperimentalRouteNames<W>& r, const pe::OuterParams& outer,
                      rz::Tally& tally) {
  const std::uint64_t t = inst.t;
  rz::RunSetting setting;
  setting.residual_limit = outer.residual_factor * t;
  setting.retry_cap = outer.retry_cap;
  setting.p_log2 = outer.p_log2.has_value() ? *outer.p_log2 : rz::paper_p_log2(t);
  CHECK(r.setting.t == t);
  CHECK(r.setting.p_log2 == setting.p_log2);
  CHECK(r.setting.paper_p_log2 == rz::paper_p_log2(t));
  CHECK(r.setting.paper_h_log2 == rz::paper_h_log2(t));
  CHECK(r.setting.h_log2 == (outer.h_rule == pe::HRule::kProofMinimum ? rz::proof_minimum_h_log2(t)
                                                                      : rz::paper_h_log2(t)));
  CHECK(r.setting.sample_size == outer.sample_factor * t);
  CHECK(r.setting.residual_limit == setting.residual_limit);
  CHECK(r.setting.paper_sample_size == 16 * t);
  CHECK(r.setting.paper_residual_limit == 32 * t);
  CHECK(r.setting.retry_cap == outer.retry_cap);
  const std::uint64_t hits_before = tally.cap_hits;
  const std::uint64_t fallback_before = tally.fallback_edges;
  rz::check_run(inst, dist, r.names, setting, tally);
  // The entry's report of D-10 agrees with the FallbackRecord of the shared component.
  CHECK(r.retry_cap_reached == (tally.cap_hits == hits_before + 1));
  CHECK(r.fallback_edges == tally.fallback_edges - fallback_before);
}

/// D-9 (and optionally D-10 with cap 1 + seed % 3) with the paper's P and factors on one family
/// of kLoopFamilies.
template <class W>
rz::Tally sweep_d9(const FamilySpec& spec, std::uint64_t seeds, std::size_t instances,
                   std::uint64_t salt, bool with_cap) {
  std::vector<rz::Instance<W>> inst;
  std::vector<std::vector<std::vector<W>>> dist;
  for (std::size_t i = 0; i < instances; ++i) {
    inst.push_back(rz::make_instance<W>(spec.family, spec.t, 1000 * salt + i));
    dist.push_back(rz::reference_distances<W>(spec.t, inst.back().edges));
  }
  rz::Tally tally;
  for (std::uint64_t seed = 1; seed <= seeds; ++seed) {
    const std::size_t i = static_cast<std::size_t>(seed % instances);
    pe::OuterParams outer;
    outer.h_rule = pe::HRule::kProofMinimum;
    if (with_cap) outer.retry_cap = 1 + seed % 3;
    WordRng rng(rz::rng_seed(seed, salt));
    const auto r = pe::name_routes(inst[i], rng, with_outer(outer));
    INFO(rz::family_name(spec.family) << " " << rz::weight_name<W>() << " t " << spec.t
                                      << " instance " << i << " seed " << seed);
    CHECK(r.setting.h_log2 < r.setting.paper_h_log2);  // D-9 is in effect
    check_mode_b_run(inst[i], dist[i], r, outer, tally);
    if (!with_cap) CHECK_FALSE(r.retry_cap_reached);
  }
  return tally;
}

/// Lemma 6.8(2) under D-9 (Q-42): the same assertions as in Mode A. A run of kLoopFamilies has
/// one round unless a cap stops it, so "no rejection" is "the first trial was accepted".
void check_d9_tally(const rz::Tally& t) {
  REQUIRE(t.runs_with_trials == t.runs);
  const auto a = rz::acceptance_test(t.accepted, t.trials);
  CHECK_FALSE(a.rejects_half);
  const double first = static_cast<double>(t.runs_without_rejection) / static_cast<double>(t.runs);
  CHECK(first >= 0.5 - rz::hoeffding_radius(t.runs));
}

std::string label(const char* what, const FamilySpec& spec, const char* weight,
                  std::uint64_t seeds) {
  return std::string(what) + ", " + rz::family_name(spec.family) + " t=" + std::to_string(spec.t) +
         ", " + weight + ", " + std::to_string(seeds) + " seeds";
}

template <class W>
void run_d9_family(const FamilySpec& spec, std::uint64_t seeds, std::size_t instances,
                   std::uint64_t salt) {
  const auto t = sweep_d9<W>(spec, seeds, instances, salt, false);
  WARN(rz::describe(label("Mode B D-9 NameRoutes", spec, rz::weight_name<W>(), seeds), t));
  check_d9_tally(t);
  CHECK(t.cap_hits == 0);
}

template <class W>
void run_d9_d10_family(const FamilySpec& spec, std::uint64_t seeds, std::size_t instances,
                       std::uint64_t salt) {
  const auto t = sweep_d9<W>(spec, seeds, instances, salt, true);
  WARN(rz::describe(
      label("Mode B D-9 + D-10 (cap 1..3) NameRoutes", spec, rz::weight_name<W>(), seeds), t));
  // A capped run stops after its rejections, which is a stopping rule; the anytime test of
  // Lemma 6.8(2) still applies.
  CHECK_FALSE(rz::acceptance_test(t.accepted, t.trials).rejects_half);
}

/// D-17 at small t (sample t, residual threshold 2t, P from a list), with D-9 on every second
/// seed and always with a cap (D-10), because with a small P a round may never be accepted.
/// Returns the tallies per log2 P (-1: the paper's P).
template <class W>
std::map<int, rz::Tally> sweep_d17(Family family, std::uint64_t seeds, std::uint64_t salt) {
  const std::optional<unsigned> p_values[] = {0U, 1U, 3U, 6U, 12U, std::nullopt};
  const std::uint64_t caps[] = {1, 2, 3, 4, 6, 64};
  std::map<int, rz::Tally> tallies;
  for (std::uint64_t seed = 1; seed <= seeds; ++seed) {
    const std::size_t t = 4 + static_cast<std::size_t>(seed % 17);
    const auto inst = rz::make_instance<W>(family, t, seed + 7919 * salt);
    const auto dist = rz::reference_distances<W>(t, inst.edges);
    pe::OuterParams outer;
    outer.h_rule = seed % 2 == 0 ? pe::HRule::kProofMinimum : pe::HRule::kPaper;
    outer.sample_factor = 1;
    outer.residual_factor = 2;
    outer.p_log2 = p_values[seed % 6];
    outer.retry_cap = caps[(seed / 6) % 6];
    WordRng rng(rz::rng_seed(seed, salt));
    const auto r = pe::name_routes(inst, rng, with_outer(outer));
    INFO(rz::family_name(family) << " " << rz::weight_name<W>() << " t " << t << " seed " << seed
                                 << " cap " << *outer.retry_cap);
    const int key = outer.p_log2.has_value() ? static_cast<int>(*outer.p_log2) : -1;
    check_mode_b_run(inst, dist, r, outer, tallies[key]);
  }
  return tallies;
}

template <class W>
void run_d17(std::uint64_t seeds, std::uint64_t salt_base) {
  std::uint64_t salt = salt_base;
  rz::Tally total;
  for (const auto family : {Family::kSparseRandom, Family::kComplete, Family::kCompleteEqual,
                            Family::kCompleteFewWeights, Family::kCompleteWideRange, Family::kDense,
                            Family::kMultiEdge, Family::kTwoCliques}) {
    const auto tallies = sweep_d17<W>(family, seeds, ++salt);
    for (const auto& [p, t] : tallies) {
      WARN(rz::describe(std::string("Mode B D-17 (+D-9 on even seeds, D-10 caps 1..64), ") +
                            rz::family_name(family) + " t in 4..20, " + rz::weight_name<W>() +
                            ", P = " + (p < 0 ? "paper" : "2^" + std::to_string(p)),
                        t));
      total.merge(t);
    }
  }
  WARN(
      rz::describe(std::string("Mode B D-17, all families and P, ") + rz::weight_name<W>(), total));
  // Both outcomes of D-10 occurred, so both branches of the bookkeeping checks were exercised.
  CHECK(total.cap_hits > 0);
  CHECK(total.runs_with_trials > total.cap_hits);
  CHECK(total.accepted > 0);
}

pv::EstimateParams d9_estimate_params(std::size_t t) {
  return d9_params().name_routes_params(t).estimate;
}

template <class W>
void run_estimate_d9(std::uint64_t seeds, std::uint64_t salt_base) {
  const auto eta = [](std::size_t t) { return rz::proof_minimum_h_log2(t); };
  std::uint64_t salt = salt_base;
  for (const auto family :
       {SampleFamily::kSparseDecades, SampleFamily::kSparseEqual, SampleFamily::kMultiEdge,
        SampleFamily::kComplete, SampleFamily::kDenseTies}) {
    const auto t = rz::sweep_estimate<W>(family, seeds, ++salt, d9_estimate_params, eta);
    WARN(rz::describe(std::string("Mode B D-9 Estimate, ") + rz::sample_family_name(family) + ", " +
                          rz::weight_name<W>(),
                      t));
    rz::check_estimate_tally(t);
  }
}

void run_separation_d9(std::uint64_t seeds, std::uint64_t salt_base) {
  std::uint64_t salt = salt_base;
  for (const auto family :
       {SampleFamily::kSparseDecades, SampleFamily::kSparseEqual, SampleFamily::kComplete}) {
    for (const std::uint64_t sample_seed : {3U, 11U}) {
      const auto s = rz::make_sample<std::uint64_t>(family, sample_seed);
      const auto r = rz::separation_test(s, seeds, ++salt, d9_estimate_params);
      WARN(rz::describe(std::string("Mode B D-9 Lemma 6.5, ") + rz::sample_family_name(family) +
                            " sample " + std::to_string(sample_seed) + " t=" + std::to_string(s.t),
                        r));
      CHECK(r.violations == 0);
    }
  }
}

}  // namespace

// ----------------------------------------------------------------------- fast (every CI run)

TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes: certificates and acceptance",
          "[randomized][paper_experimental][name_routes][property]") {
  std::uint64_t salt = 3000;
  for (const auto& spec : kLoopFamilies) {
    run_d9_family<std::uint64_t>(spec, 3, 1, ++salt);
    run_d9_family<double>(spec, 15, 2, ++salt);  // 7 families x 15 seeds = 105 runs
  }
}

TEST_CASE("P4-RANDOMIZED Mode B D-9 + D-10 NameRoutes: cap hits with the paper's P",
          "[randomized][paper_experimental][name_routes][property]") {
  run_d9_d10_family<double>(kLoopFamilies[0], 60, 2, 3100);
  run_d9_d10_family<double>(kLoopFamilies[5], 40, 2, 3101);
}

TEST_CASE("P4-RANDOMIZED Mode B D-17 NameRoutes: certificates, measured acceptance, cap hits",
          "[randomized][paper_experimental][name_routes][property]") {
  run_d17<std::uint64_t>(54, 3200);
  run_d17<double>(54, 3300);
}

TEST_CASE("P4-RANDOMIZED Mode B D-9 Estimate: Z >= d_S and the scan bound (Q-64)",
          "[randomized][paper_experimental][estimate][property]") {
  run_estimate_d9<std::uint64_t>(100, 3400);
  run_estimate_d9<double>(100, 3500);
}

TEST_CASE("P4-RANDOMIZED Mode B D-9 Estimate: the P3-ESTIMATE scan bound configuration (Q-64)",
          "[randomized][paper_experimental][estimate][property]") {
  // The samples, demands, seeds and the bound of "estimate: measured E[Z / d_S] and scans per
  // vertex" in tests/property/paper/test_estimate_property.cpp, run with the paper's H (control)
  // and with the H of D-9. The generator is repeated here because it is local to that file.
  using W = std::uint64_t;
  for (const bool d9 : {false, true}) {
    std::uint64_t scans = 0;
    double allowed = 0;
    for (std::uint64_t seed = 1; seed <= 60; ++seed) {
      sssp::test::SplitMix64 g(seed);
      const std::size_t t = 2 + static_cast<std::size_t>(g.below(30));
      const std::size_t m = t + static_cast<std::size_t>(g.below(2 * t));
      std::vector<sssp::Edge<W>> edges;
      for (std::size_t i = 0; i < m; ++i) {
        const auto u = static_cast<sssp::Vertex>(g.below(t));
        auto v = static_cast<sssp::Vertex>(g.below(t));
        if (u == v) v = static_cast<sssp::Vertex>((v + 1) % t);
        const auto decade = static_cast<unsigned>(g.below(4));
        W w = 1 + g.below(9);
        for (unsigned k = 0; k < decade; ++k) w *= 37;
        edges.push_back({u, v, w});
      }
      const auto q = rz::all_pairs(t);
      const auto params = d9 ? d9_estimate_params(t) : pv::EstimateParams{};
      WordRng rng(seed ^ 0xabcdefULL);
      const auto r = pv::estimate<W>(t, edges, q, rng, params);
      scans += rz::counter(r.diag, "estimate.adjacency_scans");
      allowed +=
          rz::harmonic(t) * static_cast<double>(rz::counter(r.diag, "estimate.vertex_scales"));
    }
    WARN((d9 ? "D-9 H" : "paper H")
         << ": adjacency scans " << scans << ", H_t-weighted vertex-scales " << allowed
         << ", ratio " << static_cast<double>(scans) / allowed);
    REQUIRE(scans > 0);
    CHECK(static_cast<double>(scans) <= allowed);
  }
}

TEST_CASE("P4-RANDOMIZED Mode B D-9 Estimate: separation at one scale (Lemma 6.5, Q-64)",
          "[randomized][paper_experimental][estimate][property]") {
  run_separation_d9(100, 3600);
}

// --------------------------------------------------- slow (seed sweeps, one family per case)

TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: complete, uint64",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<std::uint64_t>(kLoopFamilies[0], rz::loop_sweep_seeds<std::uint64_t>(), 4, 4100);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: complete, dyadic double",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<double>(kLoopFamilies[0], rz::loop_sweep_seeds<double>(), 4, 4150);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: complete, equal weights, uint64",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<std::uint64_t>(kLoopFamilies[1], rz::loop_sweep_seeds<std::uint64_t>(), 4, 4200);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: complete, equal weights, dyadic double",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<double>(kLoopFamilies[1], rz::loop_sweep_seeds<double>(), 4, 4250);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: complete, few weights, uint64",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<std::uint64_t>(kLoopFamilies[2], rz::loop_sweep_seeds<std::uint64_t>(), 4, 4300);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: complete, few weights, dyadic double",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<double>(kLoopFamilies[2], rz::loop_sweep_seeds<double>(), 4, 4350);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: complete, wide weight range, uint64",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<std::uint64_t>(kLoopFamilies[3], rz::loop_sweep_seeds<std::uint64_t>(), 4, 4400);
}
TEST_CASE(
    "P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: complete, wide weight range, dyadic "
    "double",
    "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<double>(kLoopFamilies[3], rz::loop_sweep_seeds<double>(), 4, 4450);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: dense, uint64",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<std::uint64_t>(kLoopFamilies[4], rz::loop_sweep_seeds<std::uint64_t>(), 4, 4500);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: dense, dyadic double",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<double>(kLoopFamilies[4], rz::loop_sweep_seeds<double>(), 4, 4550);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: multi-edge, uint64",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<std::uint64_t>(kLoopFamilies[5], rz::loop_sweep_seeds<std::uint64_t>(), 4, 4600);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: multi-edge, dyadic double",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<double>(kLoopFamilies[5], rz::loop_sweep_seeds<double>(), 4, 4650);
}
TEST_CASE("P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: two cliques (disconnected), uint64",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<std::uint64_t>(kLoopFamilies[6], rz::loop_sweep_seeds<std::uint64_t>(), 4, 4700);
}
TEST_CASE(
    "P4-RANDOMIZED Mode B D-9 NameRoutes, seed sweep: two cliques (disconnected), dyadic "
    "double",
    "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_family<double>(kLoopFamilies[6], rz::loop_sweep_seeds<double>(), 4, 4750);
}

TEST_CASE(
    "P4-RANDOMIZED Mode B D-9 + D-10 NameRoutes, seed sweep: cap hits with the paper's P, "
    "complete",
    "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_d10_family<double>(kLoopFamilies[0], rz::kSweepSeeds, 4, 4800 + kLoopFamilies[0].t);
}

TEST_CASE(
    "P4-RANDOMIZED Mode B D-9 + D-10 NameRoutes, seed sweep: cap hits with the paper's P, "
    "complete, wide weight range",
    "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_d10_family<double>(kLoopFamilies[3], rz::kSweepSeeds, 4, 4800 + kLoopFamilies[3].t);
}

TEST_CASE(
    "P4-RANDOMIZED Mode B D-9 + D-10 NameRoutes, seed sweep: cap hits with the paper's P, "
    "multi-edge",
    "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d9_d10_family<double>(kLoopFamilies[5], rz::kSweepSeeds, 4, 4800 + kLoopFamilies[5].t);
}

TEST_CASE("P4-RANDOMIZED Mode B D-17 NameRoutes, seed sweep, uint64",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d17<std::uint64_t>(rz::kSweepSeeds, 5000);
}

TEST_CASE("P4-RANDOMIZED Mode B D-17 NameRoutes, seed sweep, dyadic double",
          "[randomized][paper_experimental][name_routes][property][slow]") {
  run_d17<double>(rz::kSweepSeeds, 5100);
}

TEST_CASE("P4-RANDOMIZED Mode B D-9 Estimate, seed sweep with the scan bound (Q-64), uint64",
          "[randomized][paper_experimental][estimate][property][slow]") {
  run_estimate_d9<std::uint64_t>(rz::kSweepSeeds, 5200);
}

TEST_CASE("P4-RANDOMIZED Mode B D-9 Estimate, seed sweep with the scan bound (Q-64), dyadic double",
          "[randomized][paper_experimental][estimate][property][slow]") {
  run_estimate_d9<double>(rz::kSweepSeeds, 5300);
}

TEST_CASE("P4-RANDOMIZED Mode B D-9 Estimate, seed sweep: separation at one scale (Lemma 6.5)",
          "[randomized][paper_experimental][estimate][property][slow]") {
  run_separation_d9(rz::kSweepSeeds, 5400);
}
