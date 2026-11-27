// Mode B parameter set (src/experimental/paper_experimental/params.hpp), task P5-EXP-PARAMS.
//
// Expected values come from the paper (arXiv 2611.34911v1) and from the research note
// docs/research/thresholds.md (§2.2, §5.1-5.3), whose numbers were recomputed independently with
// Python big integers (exact products, no floating point) while writing this test:
//   * L_N = ceil(log2(N + 2)) (Sec 2 l.284); local threshold k > 2 L_k^4 first holds at 209953.
//   * h = b + 1; local K = least power of two >= 8h (Sec 5.2 l.1029).
//   * u = ceil(log2(t + 2)), g = 4u + 3, P = least power of two >= 8*256*g*u^2 (Alg 5 l.2-3).
//   * T = least power of two >= 4t, top scale l with 2^l >= tT^2, H = least power of two
//     >= (16t)^5 (Sec 6.2 l.1473-1491); D-9: H = least power of two >= 4tT^2.
//   * F = least power of two > 2P(2h+1)(8K^2h)^2 with the outer K = 16 (Sec 8.2 l.2080).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "experimental/paper_experimental/params.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "paper/local/radial_states.hpp"
#include "sssp/result.hpp"

namespace pe = sssp::paper_experimental;
using Catch::Matchers::ContainsSubstring;

namespace {

constexpr std::uint64_t kU64Max = std::numeric_limits<std::uint64_t>::max();

std::map<std::string, pe::ParamEntry> by_key(const std::vector<pe::ParamEntry>& log) {
  std::map<std::string, pe::ParamEntry> out;
  for (const auto& e : log) out.emplace(e.key, e);
  return out;
}

/// A Mode B set that changes every tunable parameter.
pe::PaperParams all_tuned() {
  pe::LocalParams local;
  local.threshold = {1, 1};
  local.separator = pe::SeparatorRule::constant(3);
  local.alignment = pe::LocalAlignment::kLambda;
  local.assembly = pe::LocalAssembly::kExactHeight;
  pe::OuterParams outer;
  outer.threshold = {0, 0};
  outer.separator = pe::SeparatorRule::power_of_l(3);
  outer.sample_factor = 1;
  outer.residual_factor = 2;
  outer.p_log2 = 1;
  outer.h_rule = pe::HRule::kProofMinimum;
  outer.retry_cap = 64;
  return pe::PaperParams(local, outer);
}

}  // namespace

TEST_CASE("experimental params: the defaults are the paper's values",
          "[paper_experimental][params]") {
  const pe::LocalParams local;
  CHECK(local.threshold.coeff == 2);  // Sec 2 l.286: k > 2 L_k^4
  CHECK(local.threshold.exponent == 4);
  CHECK(local.separator == pe::SeparatorRule::power_of_l(4));  // Sec 5.1 l.967: b = L^4
  CHECK(local.alignment == pe::LocalAlignment::kPaper);
  CHECK(local.assembly == pe::LocalAssembly::kGrid);

  const pe::OuterParams outer;
  CHECK(outer.threshold.coeff == 1);  // Sec 2 l.286: N > L_N^64
  CHECK(outer.threshold.exponent == 64);
  CHECK(outer.separator == pe::SeparatorRule::power_of_l(64));  // Sec 7.1 l.1644: b = L^64
  CHECK(outer.sample_factor == 16);                             // Alg 5 l.6
  CHECK(outer.residual_factor == 32);                           // Alg 5 l.5
  CHECK_FALSE(outer.p_log2.has_value());                        // Alg 5 l.3
  CHECK(outer.h_rule == pe::HRule::kPaper);                     // Sec 6.2 l.1491
  CHECK_FALSE(outer.retry_cap.has_value());                     // not in the paper (D-10)

  const pe::PaperParams p;
  CHECK(p.is_paper());
  CHECK(p.deviations().empty());
  CHECK(pe::PaperParams(local, outer).is_paper());
}

TEST_CASE("experimental params: constants that cannot be changed even in Mode B are rejected",
          "[paper_experimental][params]") {
  struct Case {
    const char* field;
    std::function<void(pe::LocalParams&, pe::OuterParams&, pe::HierarchyFixed&)> change;
  };
  // PAPER_MAPPING §2 rows marked "Cannot be changed even in Mode B" (and Sec 5.4 c = 10, N^16,
  // which Mode B may only replace as a whole through the D-8 variant). Both directions are tried.
  const std::vector<Case> cases = {
      {"local.shell_ratio_factor", [](auto& l, auto&, auto&) { l.fixed.shell_ratio_factor = 4; }},
      {"local.shell_ratio_factor", [](auto& l, auto&, auto&) { l.fixed.shell_ratio_factor = 16; }},
      {"local.near_ratio_exponent", [](auto& l, auto&, auto&) { l.fixed.near_ratio_exponent = 1; }},
      {"local.macro_exponent", [](auto& l, auto&, auto&) { l.fixed.macro_exponent = 2; }},
      {"local.align_exponent", [](auto& l, auto&, auto&) { l.fixed.align_exponent = 9; }},
      {"local.support_exponent", [](auto& l, auto&, auto&) { l.fixed.support_exponent = 17; }},
      {"outer.shell_ratio", [](auto&, auto& o, auto&) { o.fixed.shell_ratio = 8; }},
      {"outer.shell_ratio", [](auto&, auto& o, auto&) { o.fixed.shell_ratio = 32; }},
      {"outer.near_ratio", [](auto&, auto& o, auto&) { o.fixed.near_ratio = 64; }},
      {"outer.scale_factor_coeff", [](auto&, auto& o, auto&) { o.fixed.scale_factor_coeff = 1; }},
      {"outer.support_cube", [](auto&, auto& o, auto&) { o.fixed.support_cube = 2; }},
      {"outer.band_factor", [](auto&, auto& o, auto&) { o.estimate.band_factor = 2; }},
      {"outer.window_exponent", [](auto&, auto& o, auto&) { o.estimate.window_exponent = 1; }},
      {"outer.grid_octaves_below", [](auto&, auto& o, auto&) { o.calendars.octaves_below = 6; }},
      {"outer.grid_octaves_above", [](auto&, auto& o, auto&) { o.calendars.octaves_above = 1; }},
      {"outer.width_offset", [](auto&, auto& o, auto&) { o.calendars.width_offset = 4; }},
      {"outer.overlap_factor", [](auto&, auto& o, auto&) { o.calendars.overlap_factor = 2; }},
      {"outer.crossing_factor", [](auto&, auto& o, auto&) { o.calendars.crossing_factor = 8; }},
      {"outer.atom_factor", [](auto&, auto& o, auto&) { o.calendars.atom_factor = 3; }},
      {"hierarchy.normalization_link_factor",
       [](auto&, auto&, auto& h) { h.normalization_link_factor = 2; }},
  };
  for (const auto& c : cases) {
    INFO(c.field);
    pe::LocalParams local;
    pe::OuterParams outer;
    pe::HierarchyFixed hier;
    c.change(local, outer, hier);
    CHECK_THROWS_AS(pe::PaperParams(local, outer, hier), std::invalid_argument);
    CHECK_THROWS_WITH(pe::PaperParams(local, outer, hier), ContainsSubstring(c.field));
  }
}

TEST_CASE("experimental params: every fixed constant is rejected above and below the paper",
          "[paper_experimental][params]") {
  // Every field of LocalFixed, OuterFixed, EstimateFixed, CalendarFixed and HierarchyFixed,
  // with the value the paper fixes (PAPER_MAPPING §2). Each field is moved one step down and one
  // step up, and to the extremes 0 and UINT_MAX of its type. Every paper value is >= 1, so the
  // downward direction exists for every field; the only direction the type cannot express is
  // below 0, which is why 0 is the lowest value tried. The paper value itself must be accepted,
  // which also checks that each setter writes the field it names.
  struct Field {
    const char* name;
    unsigned paper;
    std::function<void(pe::LocalParams&, pe::OuterParams&, pe::HierarchyFixed&, unsigned)> set;
  };
  const std::vector<Field> fields = {
      {"local.shell_ratio_factor", 8,  // Sec 5.2 l.1029
       [](auto& l, auto&, auto&, unsigned v) { l.fixed.shell_ratio_factor = v; }},
      {"local.near_ratio_exponent", 2,  // Sec 5.2
       [](auto& l, auto&, auto&, unsigned v) { l.fixed.near_ratio_exponent = v; }},
      {"local.macro_exponent", 3,  // Sec 5.2 l.1021
       [](auto& l, auto&, auto&, unsigned v) { l.fixed.macro_exponent = v; }},
      {"local.align_exponent", 10,  // Sec 5.4 l.1206
       [](auto& l, auto&, auto&, unsigned v) { l.fixed.align_exponent = v; }},
      {"local.support_exponent", 16,  // Sec 5.4 l.1216
       [](auto& l, auto&, auto&, unsigned v) { l.fixed.support_exponent = v; }},
      {"outer.shell_ratio", 16,  // Sec 8.1 l.1977
       [](auto&, auto& o, auto&, unsigned v) { o.fixed.shell_ratio = v; }},
      {"outer.near_ratio", 256,  // Sec 8.1 l.1977
       [](auto&, auto& o, auto&, unsigned v) { o.fixed.near_ratio = v; }},
      {"outer.scale_factor_coeff", 2,  // Sec 8.2 l.2080
       [](auto&, auto& o, auto&, unsigned v) { o.fixed.scale_factor_coeff = v; }},
      {"outer.support_cube", 3,  // Sec 8.2
       [](auto&, auto& o, auto&, unsigned v) { o.fixed.support_cube = v; }},
      {"outer.band_factor", 4,  // Sec 6.2 l.1473
       [](auto&, auto& o, auto&, unsigned v) { o.estimate.band_factor = v; }},
      {"outer.window_exponent", 2,  // Sec 6.2 l.1490
       [](auto&, auto& o, auto&, unsigned v) { o.estimate.window_exponent = v; }},
      {"outer.grid_octaves_below", 7,  // Sec 7.3 l.1697
       [](auto&, auto& o, auto&, unsigned v) { o.calendars.octaves_below = v; }},
      {"outer.grid_octaves_above", 2,  // Sec 7.3 l.1697
       [](auto&, auto& o, auto&, unsigned v) { o.calendars.octaves_above = v; }},
      {"outer.width_offset", 5,  // Sec 7.3 l.1698
       [](auto&, auto& o, auto&, unsigned v) { o.calendars.width_offset = v; }},
      {"outer.overlap_factor", 4,  // Lemma 7.6
       [](auto&, auto& o, auto&, unsigned v) { o.calendars.overlap_factor = v; }},
      {"outer.crossing_factor", 16,  // Lemma 7.8
       [](auto&, auto& o, auto&, unsigned v) { o.calendars.crossing_factor = v; }},
      {"outer.atom_factor", 2,  // Lemma 7.5
       [](auto&, auto& o, auto&, unsigned v) { o.calendars.atom_factor = v; }},
      {"hierarchy.normalization_link_factor", 4,  // Sec 4.1
       [](auto&, auto&, auto& h, unsigned v) { h.normalization_link_factor = v; }},
  };
  constexpr unsigned kUMax = std::numeric_limits<unsigned>::max();
  for (const auto& f : fields) {
    INFO(f.name);
    REQUIRE(f.paper >= 1);  // the downward direction exists
    {
      pe::LocalParams local;
      pe::OuterParams outer;
      pe::HierarchyFixed hier;
      f.set(local, outer, hier, f.paper);
      CHECK_NOTHROW(pe::PaperParams(local, outer, hier));
    }
    for (const unsigned v : {f.paper - 1, f.paper + 1, 0U, kUMax}) {
      INFO("value " << v);
      pe::LocalParams local;
      pe::OuterParams outer;
      pe::HierarchyFixed hier;
      f.set(local, outer, hier, v);
      CHECK_THROWS_AS(pe::PaperParams(local, outer, hier), std::invalid_argument);
      CHECK_THROWS_WITH(pe::PaperParams(local, outer, hier), ContainsSubstring(f.name));
    }
  }
}

TEST_CASE("experimental params: separator exponents above 64 are rejected",
          "[paper_experimental][params]") {
  // The paper's largest separator is the outer b = L^64 (Sec 7.1 l.1644); PaperParams documents
  // "a separator exponent above 64" as outside the domain, for the local and the outer rule.
  for (const unsigned e : {65U, 66U, 128U, std::numeric_limits<unsigned>::max()}) {
    INFO("e = " << e);
    pe::LocalParams l;
    l.separator = pe::SeparatorRule::power_of_l(e);
    CHECK_THROWS_AS(pe::PaperParams(l, pe::OuterParams{}), std::invalid_argument);
    CHECK_THROWS_WITH(
        pe::PaperParams(l, pe::OuterParams{}),
        ContainsSubstring("local.separator") && ContainsSubstring("L^" + std::to_string(e)));
    pe::OuterParams o;
    o.separator = pe::SeparatorRule::power_of_l(e);
    CHECK_THROWS_AS(pe::PaperParams(pe::LocalParams{}, o), std::invalid_argument);
    CHECK_THROWS_WITH(
        pe::PaperParams(pe::LocalParams{}, o),
        ContainsSubstring("outer.separator") && ContainsSubstring("L^" + std::to_string(e)));
  }
  // The boundary: L^64 (the paper's outer b) and L^0 are accepted for both rules.
  for (const unsigned e : {0U, 1U, 4U, 63U, 64U}) {
    INFO("e = " << e);
    pe::LocalParams l;
    l.separator = pe::SeparatorRule::power_of_l(e);
    CHECK_NOTHROW(pe::PaperParams(l, pe::OuterParams{}));
    pe::OuterParams o;
    o.separator = pe::SeparatorRule::power_of_l(e);
    CHECK_NOTHROW(pe::PaperParams(pe::LocalParams{}, o));
  }
}

TEST_CASE("experimental params: invalid tunable values are rejected",
          "[paper_experimental][params]") {
  {
    pe::LocalParams l;
    l.separator = pe::SeparatorRule::constant(0);  // Def 5.1 needs b >= 1
    CHECK_THROWS_AS(pe::PaperParams(l, pe::OuterParams{}), std::invalid_argument);
  }
  {
    pe::OuterParams o;
    o.separator = pe::SeparatorRule::constant(0);
    CHECK_THROWS_AS(pe::PaperParams(pe::LocalParams{}, o), std::invalid_argument);
  }
  {
    pe::OuterParams o;
    o.sample_factor = 0;
    CHECK_THROWS_AS(pe::PaperParams(pe::LocalParams{}, o), std::invalid_argument);
  }
  {
    pe::OuterParams o;
    o.sample_factor = 8;
    o.residual_factor = 7;  // the loop must stop before |U| < |S|
    CHECK_THROWS_AS(pe::PaperParams(pe::LocalParams{}, o), std::invalid_argument);
  }
  {
    pe::OuterParams o;
    o.retry_cap = 0;
    CHECK_THROWS_AS(pe::PaperParams(pe::LocalParams{}, o), std::invalid_argument);
  }
  {
    pe::OuterParams o;
    o.p_log2 = 257;  // NameRoutes accepts p_log2 <= 256
    CHECK_THROWS_AS(pe::PaperParams(pe::LocalParams{}, o), std::invalid_argument);
  }
}

TEST_CASE("experimental params: size thresholds follow N > c * L_N^e",
          "[paper_experimental][params]") {
  const pe::PaperParams paper;
  // Sec 2 l.286 local: least k is 209953 (L = 18, 2 * 18^4 = 209952).
  CHECK_FALSE(paper.local_applies(209952));
  CHECK(paper.local_applies(209953));
  // Outer: L_N^64 >= 2^64 > N for every N < 2^64.
  for (const std::uint64_t n :
       {std::uint64_t{1}, std::uint64_t{1} << 40, (std::uint64_t{1} << 63) + 5, kU64Max}) {
    CHECK_FALSE(paper.outer_applies(n));
  }

  pe::OuterParams always;
  always.threshold = {0, 0};  // N > 0: always the outer construction (thresholds.md §7.2)
  pe::OuterParams four;
  four.threshold = {4, 0};  // N > 4
  CHECK(pe::PaperParams(pe::LocalParams{}, always).outer_applies(1));
  CHECK_FALSE(pe::PaperParams(pe::LocalParams{}, four).outer_applies(4));
  CHECK(pe::PaperParams(pe::LocalParams{}, four).outer_applies(5));
}

TEST_CASE("experimental params: separator rules", "[paper_experimental][params]") {
  // L(209953) = 18, L(10^6) = 20, L(10^8) = 27.
  CHECK(pe::SeparatorRule::power_of_l(4).evaluate(209953).value == 104976);
  CHECK_FALSE(pe::SeparatorRule::power_of_l(4).evaluate(209953).saturated);
  CHECK(pe::SeparatorRule::power_of_l(3).evaluate(1000000).value == 8000);
  CHECK(pe::SeparatorRule::power_of_l(4).evaluate(100000000).value == 531441);
  CHECK(pe::SeparatorRule::power_of_l(0).evaluate(5).value == 1);
  CHECK(pe::SeparatorRule::constant(3).evaluate(1000000).value == 3);
  // 20^64 >= 2^64: saturates and says so.
  const auto sat = pe::SeparatorRule::power_of_l(64).evaluate(1000000);
  CHECK(sat.saturated);
  CHECK(sat.value == kU64Max);
}

TEST_CASE("experimental params: local derived constants are recomputed from b",
          "[paper_experimental][params]") {
  // Sec 5.3 of thresholds.md: k = 209953, b = 104976, h = 104977, K = 2^20.
  const auto paper = pe::PaperParams{}.local_derived(209953);
  CHECK(paper.separator == 104976);
  CHECK(paper.hop_bound == 104977);
  CHECK_FALSE(paper.saturated);
  CHECK(paper.shell_ratio_log2 == 20);

  // K = least power of two >= 8(b + 1).
  const std::vector<std::pair<std::uint64_t, unsigned>> expected = {
      {1, 4}, {2, 5}, {3, 5}, {4, 6}, {7, 6}, {8, 7}, {160000, 21}};
  for (const auto& [b, k_log2] : expected) {
    INFO("b = " << b);
    pe::LocalParams l;
    l.separator = pe::SeparatorRule::constant(b);
    const auto d = pe::PaperParams(l, pe::OuterParams{}).local_derived(1000);
    CHECK(d.separator == b);
    CHECK(d.hop_bound == b + 1);
    CHECK(d.shell_ratio_log2 == k_log2);
    // Same K as the shared Mode A component computes from h.
    CHECK(d.shell_ratio_log2 == sssp::paper::local_shell_ratio_log2(b + 1));
  }

  // b = 2^64 - 1: h = 2^64 does not fit, K = 2^67 is still exact.
  pe::LocalParams huge;
  huge.separator = pe::SeparatorRule::constant(kU64Max);
  const auto d = pe::PaperParams(huge, pe::OuterParams{}).local_derived(1000);
  CHECK(d.saturated);
  CHECK(d.hop_bound == kU64Max);
  CHECK(d.shell_ratio_log2 == 67);
}

TEST_CASE("experimental params: outer derived constants are recomputed from b and t",
          "[paper_experimental][params]") {
  struct Row {
    std::uint64_t n;
    pe::SeparatorRule b;
    std::uint64_t t;
    unsigned u, g, p_log2, t_log2, top_scale, h_log2, f_log2;
  };
  // thresholds.md §5.2 (N, b, t, u, g, P, T, H, l_a, F).
  const std::vector<Row> rows = {
      {1000000, pe::SeparatorRule::power_of_l(1), 99999, 17, 71, 26, 19, 55, 104, 64},
      {1000000, pe::SeparatorRule::constant(160000), 11, 4, 19, 20, 6, 16, 38, 96},
      {100000000, pe::SeparatorRule::power_of_l(4), 375, 9, 39, 23, 11, 31, 63, 105},
  };
  for (const auto& r : rows) {
    INFO("N = " << r.n << ", t = " << r.t);
    pe::OuterParams o;
    o.separator = r.b;
    const auto d = pe::PaperParams(pe::LocalParams{}, o).outer_derived(r.n, r.t);
    CHECK(d.hop_bound == d.separator + 1);
    CHECK(d.u == r.u);
    CHECK(d.g == r.g);
    CHECK(d.p_log2 == r.p_log2);
    CHECK(d.t_log2 == r.t_log2);
    CHECK(d.top_scale == r.top_scale);
    CHECK(d.h_log2 == r.h_log2);
    CHECK(d.f_log2 == r.f_log2);
  }

  // The paper's b = L^64 at N = 10^6 is 20^64 (saturated as a word); F is still exact:
  // t = 1 gives u = 2, g = 11, P = 2^17 (thresholds.md §5.1) and F = 2^871 (Python big ints).
  const auto sat = pe::PaperParams{}.outer_derived(1000000, 1);
  CHECK(sat.saturated);
  CHECK(sat.p_log2 == 17);
  CHECK(sat.f_log2 == 871);

  // D-9: H = least power of two >= 4tT^2; t = 99999 gives 2^57.
  pe::OuterParams d9;
  d9.h_rule = pe::HRule::kProofMinimum;
  CHECK(pe::PaperParams(pe::LocalParams{}, d9).outer_derived(1000000, 99999).h_log2 == 57);

  // D-17: F is derived from the P actually used (P = 2, b = 160000, t = 11 gives F = 2^77).
  pe::OuterParams small_p;
  small_p.separator = pe::SeparatorRule::constant(160000);
  small_p.p_log2 = 1;
  const auto sp = pe::PaperParams(pe::LocalParams{}, small_p).outer_derived(1000000, 11);
  CHECK(sp.p_log2 == 1);
  CHECK(sp.f_log2 == 77);

  CHECK_THROWS_AS(pe::PaperParams{}.outer_derived(1000000, 0), std::invalid_argument);
}

TEST_CASE("experimental params: local branch decision per instance",
          "[paper_experimental][params]") {
  const pe::PaperParams paper;
  CHECK(paper.local_branch(209952) == pe::LocalBranch::kBaseCase);
  CHECK(paper.local_branch(209953) == pe::LocalBranch::kRadial);

  // Conservative variant (Q-03): the radial branch needs k >= max(32, b + 1); otherwise the
  // instance is unsupported (a failure, never a silent substitution).
  pe::LocalParams l;
  l.threshold = {1, 1};  // k > L_k
  l.separator = pe::SeparatorRule::constant(3);
  const pe::PaperParams small(l, pe::OuterParams{});
  CHECK(small.local_branch(3) == pe::LocalBranch::kBaseCase);  // 3 <= L_3 = 3
  CHECK(small.local_branch(31) == pe::LocalBranch::kUnsupported);
  CHECK(small.local_branch(32) == pe::LocalBranch::kRadial);
  l.separator = pe::SeparatorRule::constant(40);
  const pe::PaperParams wide(l, pe::OuterParams{});
  CHECK(wide.local_branch(40) == pe::LocalBranch::kUnsupported);
  CHECK(wide.local_branch(41) == pe::LocalBranch::kRadial);

  // D-8 (Lambda variant): the claim of research note 02 is N >= 2 and b >= 1.
  l.threshold = {0, 0};
  l.separator = pe::SeparatorRule::constant(1);
  l.alignment = pe::LocalAlignment::kLambda;
  const pe::PaperParams lambda(l, pe::OuterParams{});
  CHECK(lambda.local_branch(1) == pe::LocalBranch::kUnsupported);
  CHECK(lambda.local_branch(2) == pe::LocalBranch::kRadial);
}

TEST_CASE("experimental params: D-11 accepts the radial branch for k >= 2 under either alignment",
          "[paper_experimental][params]") {
  // Q-60 (lead, 2026-11-27): the exact-height assembly (D-11) does not use the alignment, so the
  // Sec 5.4 l.1206 precondition k >= max(32, b + 1) of the paper's alignment does not apply; D-11
  // needs N >= 2 only (local_assembly_variant_applies, local_assembly.hpp), like D-8.
  for (const std::uint64_t b : {std::uint64_t{1}, std::uint64_t{3}, std::uint64_t{40}}) {
    for (const auto alignment : {pe::LocalAlignment::kPaper, pe::LocalAlignment::kLambda}) {
      INFO("b = " << b << " alignment " << static_cast<int>(alignment));
      pe::LocalParams l;
      l.threshold = {0, 0};  // every k >= 1 passes the size test
      l.separator = pe::SeparatorRule::constant(b);
      l.alignment = alignment;
      l.assembly = pe::LocalAssembly::kExactHeight;
      const pe::PaperParams p(l, pe::OuterParams{});
      CHECK(p.local_branch(1) == pe::LocalBranch::kUnsupported);  // N >= 2 (D-11)
      for (std::uint64_t k = 2; k <= 64; ++k) {
        INFO("k = " << k);
        CHECK(p.local_branch(k) == pe::LocalBranch::kRadial);
      }
    }
  }
  // The size test still decides first: k <= L_k takes the base case.
  pe::LocalParams l;
  l.threshold = {1, 1};
  l.separator = pe::SeparatorRule::constant(3);
  l.assembly = pe::LocalAssembly::kExactHeight;
  const pe::PaperParams exact(l, pe::OuterParams{});
  CHECK(exact.local_branch(3) == pe::LocalBranch::kBaseCase);  // 3 <= L_3 = 3
  CHECK(exact.local_branch(4) == pe::LocalBranch::kRadial);
  // The paper's alignment with the grid keeps the Sec 5.4 l.1206 precondition.
  l.assembly = pe::LocalAssembly::kGrid;
  const pe::PaperParams grid(l, pe::OuterParams{});
  CHECK(grid.local_branch(4) == pe::LocalBranch::kUnsupported);
  CHECK(grid.local_branch(31) == pe::LocalBranch::kUnsupported);
  CHECK(grid.local_branch(32) == pe::LocalBranch::kRadial);
}

TEST_CASE("experimental params: the alignment is logged as ignored under D-11",
          "[paper_experimental][params]") {
  // Q-61 (lead): both alignments stay selectable with D-11, but D-11 does not use the alignment;
  // the parameter log says so in the value. The paper's alignment under D-11 is not a deviation
  // (the only deviation is the assembly); the Lambda alignment keeps its D-8 tag.
  pe::LocalParams l;
  l.assembly = pe::LocalAssembly::kExactHeight;
  {
    const pe::PaperParams p(l, pe::OuterParams{});
    const auto log = by_key(p.parameter_log());
    CHECK(log.at("local.alignment").value == "paper (ignored under D-11)");
    CHECK(log.at("local.alignment").paper == "paper");
    CHECK(log.at("local.alignment").deviation.empty());
    const auto dev = p.deviations();
    REQUIRE(dev.size() == 1);
    CHECK(dev[0].key == "local.assembly");
    CHECK(dev[0].deviation == "D-11");
#if SSSP_DIAGNOSTICS
    sssp::Diagnostics diag;
    p.record(diag);
    CHECK(diag.counters.at("params.local.alignment=paper (ignored under D-11)") == 1);
    CHECK(diag.counters.count("params.local.alignment=paper") == 0);
    CHECK(diag.counters.count("params.deviation.D-8.local.alignment") == 0);
    CHECK(diag.counters.at("params.deviation.D-11.local.assembly") == 1);
    CHECK(diag.counters.at("params.deviations") == 1);
#endif
  }
  l.alignment = pe::LocalAlignment::kLambda;
  {
    const pe::PaperParams p(l, pe::OuterParams{});
    const auto log = by_key(p.parameter_log());
    CHECK(log.at("local.alignment").value == "lambda (ignored under D-11)");
    CHECK(log.at("local.alignment").deviation == "D-8");
#if SSSP_DIAGNOSTICS
    sssp::Diagnostics diag;
    p.record(diag);
    CHECK(diag.counters.at("params.local.alignment=lambda (ignored under D-11)") == 1);
    CHECK(diag.counters.at("params.deviation.D-8.local.alignment") == 1);
#endif
  }
  // Under the grid the alignment is used: no suffix.
  l.assembly = pe::LocalAssembly::kGrid;
  CHECK(by_key(pe::PaperParams(l, pe::OuterParams{}).parameter_log()).at("local.alignment").value ==
        "lambda");
  CHECK(by_key(pe::PaperParams{}.parameter_log()).at("local.alignment").value == "paper");
}

TEST_CASE("experimental params: bridge to the shared local component",
          "[paper_experimental][params]") {
  // The paper's values map to the paper's LocalHierarchyParams exactly.
  for (const std::uint64_t k : {std::uint64_t{1}, std::uint64_t{209953}, std::uint64_t{1} << 40}) {
    CHECK(pe::PaperParams{}.local_hierarchy_params(k) == sssp::paper::kPaperLocalParams);
  }
  pe::LocalParams l;
  l.threshold = {1, 1};
  l.separator = pe::SeparatorRule::constant(3);
  CHECK(pe::PaperParams(l, pe::OuterParams{}).local_hierarchy_params(100) ==
        sssp::paper::LocalHierarchyParams{1, 1, 3});
  l.separator = pe::SeparatorRule::power_of_l(2);  // L(100) = 7
  CHECK(pe::PaperParams(l, pe::OuterParams{}).local_hierarchy_params(100) ==
        sssp::paper::LocalHierarchyParams{1, 1, 49});
  // D-8 and D-11 are not provided by the shared Mode A component (P5-EXP-LOCAL).
  l.alignment = pe::LocalAlignment::kLambda;
  CHECK_THROWS_AS(pe::PaperParams(l, pe::OuterParams{}).local_hierarchy_params(100),
                  std::logic_error);
  l.alignment = pe::LocalAlignment::kPaper;
  l.assembly = pe::LocalAssembly::kExactHeight;
  CHECK_THROWS_AS(pe::PaperParams(l, pe::OuterParams{}).local_hierarchy_params(100),
                  std::logic_error);
}

TEST_CASE("experimental params: bridge to NameRoutes", "[paper_experimental][params]") {
  const auto paper = pe::PaperParams{}.name_routes_params(99999);
  CHECK(paper.sample_factor == 16);
  CHECK(paper.residual_factor == 32);
  CHECK_FALSE(paper.p_log2.has_value());
  CHECK_FALSE(paper.retry_cap.has_value());
  CHECK_FALSE(paper.estimate.h_log2.has_value());

  const auto tuned = all_tuned().name_routes_params(99999);
  CHECK(tuned.sample_factor == 1);
  CHECK(tuned.residual_factor == 2);
  CHECK(tuned.p_log2 == 1U);
  CHECK(tuned.retry_cap == 64U);
  CHECK(tuned.estimate.h_log2 == 57U);  // D-9 at t = 99999
}

TEST_CASE("experimental params: every deviation from the paper is listed with its ID",
          "[paper_experimental][params]") {
  const pe::PaperParams p = all_tuned();
  CHECK_FALSE(p.is_paper());
  const auto log = by_key(p.parameter_log());
  const std::vector<std::pair<std::string, std::string>> expected = {
      {"local.threshold", "D-7"},      {"local.separator", "D-7"},
      {"local.alignment", "D-8"},      {"local.assembly", "D-11"},
      {"outer.threshold", "D-7"},      {"outer.separator", "D-7"},
      {"outer.sample_factor", "D-17"}, {"outer.residual_factor", "D-17"},
      {"outer.p_log2", "D-17"},        {"outer.h_rule", "D-9"},
      {"outer.retry_cap", "D-10"},
  };
  REQUIRE(log.size() == expected.size());
  for (const auto& [key, id] : expected) {
    INFO(key);
    REQUIRE(log.count(key) == 1);
    CHECK(log.at(key).deviation == id);
    CHECK(log.at(key).value != log.at(key).paper);
  }
  CHECK(log.at("local.threshold").value == "1*L^1");
  CHECK(log.at("local.threshold").paper == "2*L^4");
  CHECK(log.at("local.separator").value == "3");
  CHECK(log.at("outer.separator").value == "L^3");
  CHECK(log.at("outer.separator").paper == "L^64");
  CHECK(log.at("outer.retry_cap").paper == "off");
  CHECK(log.at("outer.retry_cap").value == "64");
  CHECK(p.deviations().size() == expected.size());

  // The paper's set: every entry equals the paper and has no deviation ID.
  for (const auto& e : pe::PaperParams{}.parameter_log()) {
    INFO(e.key);
    CHECK(e.value == e.paper);
    CHECK(e.deviation.empty());
  }

  // Only the changed entries deviate.
  pe::LocalParams l;
  l.threshold = {1, 4};
  const auto one = pe::PaperParams(l, pe::OuterParams{}).deviations();
  REQUIRE(one.size() == 1);
  CHECK(one[0].key == "local.threshold");
  CHECK(one[0].deviation == "D-7");
}

TEST_CASE("experimental params: a run records its parameter set in Diagnostics",
          "[paper_experimental][params]") {
#if SSSP_DIAGNOSTICS
  sssp::Diagnostics paper_diag;
  pe::PaperParams{}.record(paper_diag);
  CHECK(paper_diag.counters.at("params.deviations") == 0);
  CHECK(paper_diag.counters.at("params.local.threshold=2*L^4") == 1);
  CHECK(paper_diag.counters.at("params.outer.separator=L^64") == 1);
  CHECK(paper_diag.counters.at("params.outer.retry_cap=off") == 1);
  for (const auto& [key, value] : paper_diag.counters) {
    INFO(key);
    CHECK(key.rfind("params.deviation.", 0) == std::string::npos);
  }
  CHECK(paper_diag.fallbacks.empty());  // a parameter choice is not a fallback

  sssp::Diagnostics diag;
  all_tuned().record(diag);
  CHECK(diag.counters.at("params.deviations") == 11);
  CHECK(diag.counters.at("params.local.separator=3") == 1);
  CHECK(diag.counters.at("params.deviation.D-7.local.threshold") == 1);
  CHECK(diag.counters.at("params.deviation.D-8.local.alignment") == 1);
  CHECK(diag.counters.at("params.deviation.D-10.outer.retry_cap") == 1);
  CHECK(diag.fallbacks.empty());
#else
  SUCCEED("parameter counters are checked in builds with SSSP_DIAGNOSTICS=1");
#endif
}
