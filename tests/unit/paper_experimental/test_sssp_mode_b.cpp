// Mode B (paper-experimental) entry to Algorithm 8 (arXiv 2611.34911v1, Sec 8.3):
// paper_experimental::paper_sssp runs zero-edge contraction (l.2), the size test of l.4 under
// PaperParams, then either the recorded base case (l.5) or the shared outer construction
// (l.8-15, src/paper/outer/shortest_paths.hpp) and the restoration of Lemma 8.5 (l.16-17).
//
// Checked here: the branch taken and its records (PL-6), the parameter log, the mapping of
// PaperParams to the shared component's parameters, the hand levels L0-L4 end to end (PL-3),
// and that an unsupported parameter set is a failure (an exception), never another search.
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "experimental/paper_experimental/params.hpp"
#include "experimental/paper_experimental/sssp.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "paper/outer/base_case.hpp"
#include "paper/outer/shortest_paths.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/paper.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"
#include "unit/paper/outer/shortest_paths_test_support.hpp"

using sssp::CsrGraph;
using sssp::Vertex;
using sssp::WeightTraits;
namespace pe = sssp::paper_experimental;
namespace tp = sssp::test::pipeline;

namespace {

/// The minimal forcing set of note 04 Sec 3.4: outer threshold 0 * L^0 and a constant b;
/// everything else the paper's.
pe::PaperParams forced(std::uint64_t b) {
  pe::OuterParams o;
  o.threshold = pe::SizeThreshold{0, 0};
  o.separator = pe::SeparatorRule::constant(b);
  return {pe::LocalParams{}, o};
}

}  // namespace

TEMPLATE_TEST_CASE("Mode B paper_sssp: forced outer branch on the levels L0-L4 (PL-3)",
                   "[paper_experimental][sssp][pipeline]", std::uint64_t, double) {
  using W = TestType;
  for (const auto& lv : tp::levels<W>()) {
    INFO(lv.name);
    const auto g = CsrGraph<W>::from_edges(lv.n, lv.edges);
    const auto params = forced(lv.b);
    for (Vertex s = 0; s < lv.n; ++s) {
      const auto r = pe::paper_sssp(g, s, params);
      if (s == 0) CHECK(r.dist == lv.dist_from_0);
      const auto bf = sssp::baseline::bellman_ford(g, s);
      CHECK(r.dist == bf.dist);
      const auto report = sssp::verify_sssp(g, s, r, &bf.dist);
      INFO(report.summary());
      REQUIRE(report.ok);
      // PL-6: no outer base case; only Alg 4's recorded base cases on atoms.
      const std::string local = sssp::paper::local_base_case_reason(sssp::paper::kPaperLocalParams);
      for (const auto& fb : r.diag.fallbacks) CHECK(fb.reason == local);
      // The parameter set is part of the result, also with SSSP_DIAGNOSTICS = 0.
      CHECK(r.diag.counters.at("params.deviations") == 2);
      CHECK(r.diag.counters.at("params.outer.threshold=0*L^0") == 1);
      CHECK(r.diag.counters.at("params.outer.separator=" + std::to_string(lv.b)) == 1);
      CHECK(r.diag.counters.at("params.deviation.D-7.outer.threshold") == 1);
#if SSSP_DIAGNOSTICS
      CHECK(r.diag.counters.at("outer.constructions") == 1);
      CHECK(r.diag.counters.at("outer.vertices") == lv.n);
      CHECK(r.diag.counters.count("outer_base_cases") == 0);
#endif
    }
  }
}

TEMPLATE_TEST_CASE("Mode B paper_sssp with the paper's parameters takes the recorded base case",
                   "[paper_experimental][sssp][fallback]", std::uint64_t, double) {
  using W = TestType;
  const auto g = CsrGraph<W>::from_edges(40, tp::connected_edges<W>(40, 60, 1, 30, 3));
  const pe::PaperParams paper;
  REQUIRE(paper.is_paper());
  const auto r = pe::paper_sssp(g, Vertex{7}, paper);
  const auto a = sssp::paper::paper_sssp(g, Vertex{7});
  CHECK(r.dist == a.dist);
  CHECK(r.parent == a.parent);
  REQUIRE(r.diag.fallbacks.size() == 1);
  CHECK(r.diag.fallbacks[0].reason == sssp::paper::outer::kOuterBaseCaseScanReason);
  CHECK(r.diag.fallbacks[0].vertices == 40);
  CHECK(r.diag.counters.at("params.deviations") == 0);
  CHECK(r.diag.counters.count("outer.constructions") == 0);
}

TEST_CASE("Mode B paper_sssp: a threshold between sizes selects the branch per input",
          "[paper_experimental][sssp][fallback]") {
  using W = std::uint64_t;
  // N > 5 * L_N^1: false for N = 20 (L = 5, 25), true for N = 60 (L = 6, 30).
  pe::OuterParams o;
  o.threshold = pe::SizeThreshold{5, 1};
  o.separator = pe::SeparatorRule::constant(4);
  const pe::PaperParams params{pe::LocalParams{}, o};
  const auto small = CsrGraph<W>::from_edges(20, tp::connected_edges<W>(20, 20, 1, 9, 1));
  const auto large = CsrGraph<W>::from_edges(60, tp::connected_edges<W>(60, 80, 1, 9, 2));
  const auto rs = pe::paper_sssp(small, Vertex{0}, params);
  const auto rl = pe::paper_sssp(large, Vertex{0}, params);
  CHECK(rs.dist == sssp::baseline::bellman_ford(small, Vertex{0}).dist);
  CHECK(rl.dist == sssp::baseline::bellman_ford(large, Vertex{0}).dist);
  const auto shares_small = sssp::paper::outer::fallback_shares(rs.diag);
  REQUIRE(shares_small.size() == 1);
  CHECK(shares_small[0].reason == sssp::paper::outer::kOuterBaseCaseScanReason);
  CHECK(shares_small[0].vertices == 20);
  for (const auto& share : sssp::paper::outer::fallback_shares(rl.diag)) {
    CHECK(share.reason != sssp::paper::outer::kOuterBaseCaseScanReason);
  }
}

TEST_CASE("Mode B: PaperParams map to the parameters of the shared outer construction",
          "[paper_experimental][sssp]") {
  using sssp::paper::outer::OuterSeparator;
  // The paper's values map to the paper's parameters (separator left to b = L_N^64).
  const auto paper = pe::outer_construction_params(pe::PaperParams{}, 1000);
  CHECK_FALSE(paper.separator.has_value());
  CHECK(paper.local == sssp::paper::kPaperLocalParams);
  CHECK(paper.name_routes.sample_factor == 16);
  CHECK(paper.name_routes.residual_factor == 32);
  CHECK_FALSE(paper.name_routes.p_log2.has_value());
  CHECK_FALSE(paper.name_routes.retry_cap.has_value());

  pe::OuterParams o;
  o.threshold = pe::SizeThreshold{0, 0};
  o.separator = pe::SeparatorRule::constant(3);
  o.sample_factor = 1;
  o.residual_factor = 2;
  o.p_log2 = 0;
  o.retry_cap = 5;
  pe::LocalParams l;
  l.threshold = pe::SizeThreshold{1, 2};
  l.separator = pe::SeparatorRule::constant(2);
  const auto p = pe::outer_construction_params(pe::PaperParams{l, o}, 1000);
  CHECK(p.separator == OuterSeparator{3, 1});
  CHECK(p.name_routes.sample_factor == 1);
  CHECK(p.name_routes.residual_factor == 2);
  CHECK(p.name_routes.p_log2 == 0U);
  CHECK(p.name_routes.retry_cap == std::uint64_t{5});
  CHECK(p.local.threshold_coeff == 1);
  CHECK(p.local.threshold_exponent == 2);
  CHECK(p.local.separator == 2);

  // b = L_N^e stays a power, so that F is formed from the exact h (Q-71).
  o = pe::OuterParams{};
  o.separator = pe::SeparatorRule::power_of_l(2);
  const auto q = pe::outer_construction_params(pe::PaperParams{pe::LocalParams{}, o}, 1000);
  CHECK(q.separator == OuterSeparator{10, 2});  // L_1000 = ceil(log2 1002) = 10
}

TEST_CASE("Mode B paper_sssp: an unsupported parameter set is a failure, not another search",
          "[paper_experimental][sssp][fallback]") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(12, tp::connected_edges<W>(12, 10, 1, 9, 4));
  pe::OuterParams o;
  o.threshold = pe::SizeThreshold{0, 0};
  o.separator = pe::SeparatorRule::constant(2);

  // D-8 / D-11 local variants are not available for the atoms of the outer construction.
  pe::LocalParams lambda;
  lambda.alignment = pe::LocalAlignment::kLambda;
  CHECK_THROWS_AS(pe::paper_sssp(g, Vertex{0}, pe::PaperParams{lambda, o}), std::invalid_argument);
  pe::LocalParams exact_height;
  exact_height.assembly = pe::LocalAssembly::kExactHeight;
  CHECK_THROWS_AS(pe::paper_sssp(g, Vertex{0}, pe::PaperParams{exact_height, o}),
                  std::invalid_argument);
  // A local separator L_k^e other than the paper's depends on each atom's size.
  pe::LocalParams power;
  power.separator = pe::SeparatorRule::power_of_l(2);
  CHECK_THROWS_AS(pe::paper_sssp(g, Vertex{0}, pe::PaperParams{power, o}), std::invalid_argument);
  // D-9 (H from t) needs the pivot count before NameRoutes is parameterised.
  pe::OuterParams d9 = o;
  d9.h_rule = pe::HRule::kProofMinimum;
  CHECK_THROWS_AS(pe::paper_sssp(g, Vertex{0}, pe::PaperParams{pe::LocalParams{}, d9}),
                  std::invalid_argument);
  // The same sets are fine while the outer branch is not taken (nothing unsupported runs).
  pe::OuterParams base = d9;
  base.threshold = pe::SizeThreshold{1, 64};
  const auto r = pe::paper_sssp(g, Vertex{0}, pe::PaperParams{lambda, base});
  CHECK(r.dist == sssp::baseline::bellman_ford(g, Vertex{0}).dist);
  REQUIRE(r.diag.fallbacks.size() == 1);
  CHECK(r.diag.fallbacks[0].reason == sssp::paper::outer::kOuterBaseCaseScanReason);

  CHECK_THROWS_AS(pe::paper_sssp(g, Vertex{12}, pe::PaperParams{pe::LocalParams{}, o}),
                  std::out_of_range);
}

TEMPLATE_TEST_CASE("Mode B paper_sssp: the seed is explicit and has a fixed default",
                   "[paper_experimental][sssp]", std::uint64_t, double) {
  using W = TestType;
  const auto g = CsrGraph<W>::from_edges(30, tp::connected_edges<W>(30, 40, 1, 64, 6));
  const auto params = forced(2);
  const auto a = pe::paper_sssp(g, Vertex{3}, params);
  const auto b = pe::paper_sssp(g, Vertex{3}, params, sssp::paper::kDefaultPaperSeed);
  CHECK(a.dist == b.dist);
  CHECK(a.parent == b.parent);
  CHECK(a.diag.fallbacks.size() == b.diag.fallbacks.size());
  const auto bf = sssp::baseline::bellman_ford(g, Vertex{3});
  for (const std::uint64_t seed : {std::uint64_t{0}, std::uint64_t{1}, std::uint64_t{99}}) {
    const auto r = pe::paper_sssp(g, Vertex{3}, params, seed);
    CHECK(r.dist == bf.dist);
    REQUIRE(sssp::verify_sssp(g, Vertex{3}, r, &bf.dist).ok);
  }
}
