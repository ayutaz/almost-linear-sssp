// Mode B thresholds (task P5-EXP-THRESHOLDS; arXiv 2611.34911v1 Sec 2 l.286, Alg 4 l.2,
// Alg 8 l.4, Sec 5.4 l.1206; deviation D-7, limits of Q-03).
//
// Checked here:
//   * the least size at which a size test "k > c * L_k^e" holds, against values computed
//     outside this code (docs/research/thresholds.md §2.2-2.3; the others by hand from
//     L_k = ceil(log2(k + 2)), shown next to each value) and against a brute-force scan of
//     SizeThreshold::applies;
//   * the least local size at which Alg 4 would take the radial branch, and the least size at
//     which the paper's grid assembly is unsupported (Sec 5.4 l.1206: k >= 32 and h = b + 1 <= k),
//     against a brute-force scan of PaperParams::local_branch;
//   * the Q-03 conservative local setting k >= max(32, b + 1) and the named Q-03 settings;
//   * paper_sssp rejects, before running, a local setting under which some atom size would be
//     unsupported, and records in Diagnostics which branch it took (also with
//     SSSP_DIAGNOSTICS = 0).
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "experimental/paper_experimental/params.hpp"
#include "experimental/paper_experimental/sssp.hpp"
#include "experimental/paper_experimental/thresholds.hpp"
#include "paper/outer/base_case.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "unit/paper/outer/shortest_paths_test_support.hpp"

using sssp::CsrGraph;
using sssp::Vertex;
namespace pe = sssp::paper_experimental;
namespace tp = sssp::test::pipeline;

namespace {

/// The least k in [1, limit] at which pred(k) holds, by scanning (the oracle).
template <class Pred>
std::optional<std::uint64_t> scan_least(std::uint64_t limit, Pred pred) {
  for (std::uint64_t k = 1; k <= limit; ++k) {
    if (pred(k)) return k;
  }
  return std::nullopt;
}

pe::PaperParams with_local(const pe::LocalParams& l) {
  return {l, pe::OuterParams{}};
}

pe::LocalParams local_of(pe::SizeThreshold t, pe::SeparatorRule b) {
  pe::LocalParams l;
  l.threshold = t;
  l.separator = b;
  return l;
}

pe::OuterParams forced_outer(std::uint64_t b) {
  pe::OuterParams o;
  o.threshold = pe::SizeThreshold{0, 0};
  o.separator = pe::SeparatorRule::constant(b);
  return o;
}

}  // namespace

TEST_CASE("Mode B thresholds: least size at which k > c * L_k^e holds",
          "[paper_experimental][thresholds]") {
  // thresholds.md §2.2-2.3 (exact values).
  CHECK(pe::least_applying_size({2, 4}) == std::uint64_t{209953});  // L = 18, 2 * 18^4 = 209952
  CHECK(pe::least_applying_size({1, 4}) == std::uint64_t{83522});   // L = 17, 17^4 = 83521
  CHECK(pe::least_applying_size({1, 8}) == std::uint64_t{14048223625217});  // L = 44, 44^8
  // By hand: c = 0 makes every k >= 1 apply; k > 31 * L^0 starts at 32.
  CHECK(pe::least_applying_size({0, 0}) == std::uint64_t{1});
  CHECK(pe::least_applying_size({31, 0}) == std::uint64_t{32});
  // k > L_k^2: L = 5 on [15, 30] gives 26..30; L = 6 on [31, 62] gives 37..62 only, so the set
  // is not upward closed (k = 31..36 does not apply).
  CHECK(pe::least_applying_size({1, 2}) == std::uint64_t{26});
  for (std::uint64_t k = 31; k <= 36; ++k) CHECK_FALSE(pe::SizeThreshold{1, 2}.applies(k));
  // k > L_k^3: L = 10 on [511, 1022] gives 1001..1022; L = 11 on [1023, 2046] only from 1332.
  CHECK(pe::least_applying_size({1, 3}) == std::uint64_t{1001});
  // The paper's outer test never holds below 2^64 (Sec 2 l.286; L^64 >= 2^64 > N).
  CHECK_FALSE(pe::least_applying_size({1, 64}).has_value());
  CHECK_FALSE(pe::least_applying_size({1, 16}).has_value());  // least N ~ 3.97e32

  for (const pe::SizeThreshold t :
       {pe::SizeThreshold{1, 2}, pe::SizeThreshold{2, 2}, pe::SizeThreshold{3, 1},
        pe::SizeThreshold{1, 3}, pe::SizeThreshold{100, 1}, pe::SizeThreshold{7, 0}}) {
    INFO(t.describe());
    CHECK(pe::least_applying_size(t) ==
          scan_least(5000, [&](std::uint64_t k) { return t.applies(k); }));
  }
}

TEST_CASE("Mode B thresholds: least radial and least unsupported local size",
          "[paper_experimental][thresholds]") {
  // The paper's local setting: radial from 209953, never unsupported (k > 2L^4 implies k >= 32
  // and k > L^4 = b).
  CHECK(pe::least_radial_local_size(pe::PaperParams{}) == std::uint64_t{209953});
  CHECK_FALSE(pe::least_unsupported_local_size(pe::PaperParams{}).has_value());

  // Lowered thresholds that the paper's grid cannot serve at every size (Sec 5.4 l.1206).
  // k > L^2 with b = 2: k = 26 applies but 26 < 32.
  const auto low = with_local(local_of({1, 2}, pe::SeparatorRule::constant(2)));
  CHECK(pe::least_unsupported_local_size(low) == std::uint64_t{26});
  CHECK(pe::least_radial_local_size(low) == std::uint64_t{37});
  // k > 31 with b = L_k^4: b = 6^4 = 1296 >= k at k = 32.
  const auto l4 = with_local(local_of({31, 0}, pe::SeparatorRule::power_of_l(4)));
  CHECK(pe::least_unsupported_local_size(l4) == std::uint64_t{32});
  CHECK(pe::least_radial_local_size(l4) == std::uint64_t{83522});  // k > L_k^4, thresholds.md §2.3
  // D-8 with every k >= 1 radial: k = 1 violates N >= 2.
  auto lambda = local_of({0, 0}, pe::SeparatorRule::constant(3));
  lambda.alignment = pe::LocalAlignment::kLambda;
  CHECK(pe::least_unsupported_local_size(with_local(lambda)) == std::uint64_t{1});
  CHECK(pe::least_radial_local_size(with_local(lambda)) == std::uint64_t{2});

  // Brute force against local_branch on k <= 3000 for settings whose interesting sizes lie
  // there.
  std::vector<pe::PaperParams> cases{low, l4, with_local(lambda)};
  for (const std::uint64_t b : {1, 4, 31, 32, 100, 700}) {
    cases.push_back(with_local(local_of({1, 2}, pe::SeparatorRule::constant(b))));
    cases.push_back(with_local(local_of({3, 1}, pe::SeparatorRule::constant(b))));
    cases.push_back(with_local(local_of({b, 0}, pe::SeparatorRule::constant(b))));
  }
  // A value above the scanned range must show up as "none in [1, 3000]" in the scan.
  constexpr std::uint64_t kScan = 3000;
  const auto within = [&](std::optional<std::uint64_t> v) {
    return v && *v <= kScan ? v : std::nullopt;
  };
  for (const auto& p : cases) {
    INFO(p.local().threshold.describe() << " b = " << p.local().separator.describe());
    CHECK(within(pe::least_unsupported_local_size(p)) == scan_least(kScan, [&](std::uint64_t k) {
            return p.local_branch(k) == pe::LocalBranch::kUnsupported;
          }));
    CHECK(within(pe::least_radial_local_size(p)) == scan_least(kScan, [&](std::uint64_t k) {
            return p.local_branch(k) == pe::LocalBranch::kRadial;
          }));
  }
}

TEST_CASE("Mode B thresholds: the Q-03 conservative local setting",
          "[paper_experimental][thresholds]") {
  // Q-03 decision: the paper's formulas with k >= max(32, b + 1), i.e. k > max(31, b) * L^0.
  for (const std::uint64_t b : {1, 4, 31, 32, 40, 1000}) {
    INFO("b = " << b);
    const pe::LocalParams l = pe::conservative_local_params(b);
    CHECK(l.threshold == pe::SizeThreshold{b > 31 ? b : 31, 0});
    CHECK(l.separator == pe::SeparatorRule::constant(b));
    CHECK(l.alignment == pe::LocalAlignment::kPaper);
    CHECK(l.assembly == pe::LocalAssembly::kGrid);
    const pe::PaperParams p = with_local(l);
    CHECK_FALSE(pe::least_unsupported_local_size(p).has_value());
    CHECK(pe::least_radial_local_size(p) == std::uint64_t{b + 1 > 32 ? b + 1 : 32});
    // Only D-7 entries deviate (threshold and separator).
    for (const auto& d : p.deviations()) CHECK(d.deviation == "D-7");
  }
  CHECK_THROWS_AS(pe::conservative_local_params(0), std::invalid_argument);
}

TEST_CASE("Mode B thresholds: the named Q-03 settings", "[paper_experimental][thresholds]") {
  const auto settings = pe::q03_threshold_settings();
  REQUIRE(settings.size() >= 4);
  std::set<std::string> names;
  bool some_local = false;
  for (const auto& s : settings) {
    INFO(s.name);
    CHECK(names.insert(s.name).second);
    // Every setting can run the outer branch at every size without an unsupported atom.
    CHECK_FALSE(pe::least_unsupported_local_size(s.params).has_value());
    // The outer construction applies on feasible sizes.
    const auto n0 = pe::least_applying_size(s.params.outer().threshold);
    REQUIRE(n0.has_value());
    CHECK(*n0 == s.least_outer_size);
    CHECK(*n0 <= 1001);
    CHECK(pe::least_radial_local_size(s.params) == s.least_local_radial_size);
    // A setting changes thresholds and separators only (D-7).
    REQUIRE_FALSE(s.params.deviations().empty());
    for (const auto& d : s.params.deviations()) CHECK(d.deviation == "D-7");
    if (s.least_local_radial_size && *s.least_local_radial_size <= 1024) some_local = true;
  }
  CHECK(some_local);
}

TEST_CASE("Mode B paper_sssp rejects an unsupported local setting before the outer branch",
          "[paper_experimental][thresholds][fallback]") {
  using W = std::uint64_t;
  // 12 vertices: no atom can reach k = 26, so the run would never meet the unsupported size;
  // the parameter set is rejected anyway, since it cannot serve every input (Mode B failure
  // rule, Q-03).
  const auto g = CsrGraph<W>::from_edges(12, tp::connected_edges<W>(12, 10, 1, 9, 4));
  const pe::LocalParams low = local_of({1, 2}, pe::SeparatorRule::constant(2));
  CHECK_THROWS_AS(pe::paper_sssp(g, Vertex{0}, pe::PaperParams{low, forced_outer(4)}),
                  std::invalid_argument);
  const pe::LocalParams l4 = local_of({31, 0}, pe::SeparatorRule::power_of_l(4));
  CHECK_THROWS_AS(pe::paper_sssp(g, Vertex{0}, pe::PaperParams{l4, forced_outer(4)}),
                  std::invalid_argument);
  // While the outer branch is not taken, Alg 4 never runs, so nothing unsupported runs.
  pe::OuterParams base = forced_outer(4);
  base.threshold = pe::SizeThreshold{1, 64};
  const auto r = pe::paper_sssp(g, Vertex{0}, pe::PaperParams{low, base});
  CHECK(r.dist == sssp::baseline::bellman_ford(g, Vertex{0}).dist);
  CHECK(r.diag.counters.at("mode_b.branch.outer_base_case") == 1);
  CHECK(r.diag.counters.count("mode_b.branch.outer_construction") == 0);
}

TEST_CASE("Mode B paper_sssp records the branch of Alg 8 and of every Alg 4 call",
          "[paper_experimental][thresholds][fallback]") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(40, tp::connected_edges<W>(40, 60, 1, 30, 3));
  // The paper's parameters: the recorded outer base case, no Alg 4 call.
  const auto a = pe::paper_sssp(g, Vertex{1}, pe::PaperParams{});
  CHECK(a.diag.counters.at("mode_b.branch.outer_base_case") == 1);
  CHECK(a.diag.counters.count("mode_b.branch.outer_construction") == 0);
  CHECK(a.diag.counters.count("mode_b.local.calls") == 0);

  // Forced outer branch with the conservative local setting: every Alg 4 call is counted,
  // radial + base case = calls, and the base cases are exactly the recorded FallbackRecords.
  const pe::PaperParams p{pe::conservative_local_params(4), forced_outer(8)};
  const auto r = pe::paper_sssp(g, Vertex{1}, p);
  CHECK(r.dist == sssp::baseline::bellman_ford(g, Vertex{1}).dist);
  CHECK(r.diag.counters.at("mode_b.branch.outer_construction") == 1);
  CHECK(r.diag.counters.count("mode_b.branch.outer_base_case") == 0);
  const auto calls = r.diag.counters.at("mode_b.local.calls");
  const auto radial = r.diag.counters.at("mode_b.local.radial");
  const auto base = r.diag.counters.at("mode_b.local.base_case");
  CHECK(calls > 0);
  CHECK(radial + base == calls);
  CHECK(base == r.diag.fallbacks.size());
#if SSSP_DIAGNOSTICS
  CHECK(r.diag.counters.at("local_hierarchy_calls") == calls);
  CHECK(r.diag.counters.at("local_hierarchy_base_case") == base);
  CHECK(r.diag.counters.count("local_hierarchy_radial") == (radial > 0 ? 1U : 0U));
#endif
}
