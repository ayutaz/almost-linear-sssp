// Unit tests of the Mode B local hierarchy (task P5-EXP-LOCAL, arXiv 2611.34911v1):
//   * the entry point branches on PaperParams::local_branch (Alg 4 l.2): the base case is
//     recorded as a FallbackRecord with paper::local_base_case_reason, an unsupported instance
//     is a failure (exception), and the parameter log is written first;
//   * the paper's alignment with the grid assembly reproduces paper::local_hierarchy exactly
//     (hierarchy, maps, fallbacks and counters);
//   * the Mode B assembly engine with the paper's constants reproduces
//     paper::assemble_local_hierarchy exactly (this ties the D-8 engine to Mode A);
//   * D-8: the Lambda formulas of research note 02 §4.3, computed by hand;
//   * D-11: exact heights, links not longer than the height (c0 = 1).

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "experimental/paper_experimental/local_assembly.hpp"
#include "experimental/paper_experimental/local_hierarchy.hpp"
#include "experimental/paper_experimental/params.hpp"
#include "paper/local/assemble.hpp"
#include "paper/local/base_case.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/paper_experimental/local_mode_b_support.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::WeightTraits;
namespace pe = sssp::paper_experimental;
namespace lt = sssp::test::local;
namespace mb = sssp::test::mode_b;

namespace {

constexpr const char* kEntryPrefix = "paper_experimental_local_";

template <class W>
CsrGraph<W> graph_of(lt::Family f, lt::Weights wd, std::size_t n, std::uint64_t seed) {
  return CsrGraph<W>::from_edges(n, lt::make_graph<W>(f, wd, n, seed));
}

}  // namespace

TEST_CASE("Mode B local: the base case is a recorded fallback", "[paper_experimental][local]") {
  using W = std::uint64_t;
  const auto g = graph_of<W>(lt::Family::kSparse, lt::Weights::kUniform, 10, 7);
  // D-8 parameters with the paper's size test 2 L^4: k = 10 takes the base case.
  const auto params = mb::make_params(pe::LocalAlignment::kLambda, pe::LocalAssembly::kGrid, 2,
                                      pe::SizeThreshold{2, 4});
  const auto run = pe::local_hierarchy(g, params);
  CHECK(run.branch == pe::LocalBranch::kBaseCase);
  REQUIRE(run.diag.fallbacks.size() == 1);
  CHECK(run.diag.fallbacks[0].reason ==
        sssp::paper::local_base_case_reason(sssp::paper::LocalHierarchyParams{2, 4, 2}));
  CHECK(run.diag.fallbacks[0].vertices == 10);
  CHECK(run.diag.fallbacks[0].edges == g.num_edges());
  // The parameter log is part of the result (written even without diagnostics).
  CHECK(run.diag.counters.at("params.local.alignment=lambda") == 1);
  CHECK(run.diag.counters.at("params.deviation.D-8.local.alignment") == 1);
  mb::require_same_result(run.local, sssp::paper::local_hierarchy_base_case(g));
}

TEST_CASE("Mode B local: an unsupported instance is a failure, not a substitution",
          "[paper_experimental][local]") {
  using W = std::uint64_t;
  // D-8 needs N >= 2 (research note 02 §4.3); k = 1 passes the size test 0 * L^0.
  const auto one = CsrGraph<W>::from_edges(1, {});
  const auto lambda = mb::make_params(pe::LocalAlignment::kLambda, pe::LocalAssembly::kGrid, 1);
  REQUIRE(lambda.local_branch(1) == pe::LocalBranch::kUnsupported);
  CHECK_THROWS_AS(pe::local_hierarchy(one, lambda), std::invalid_argument);
  // The paper's formulas need N >= 32 and b + 1 <= N (Sec 5.4 l.1206, Q-03 conservative).
  const auto g31 = graph_of<W>(lt::Family::kPath, lt::Weights::kUniform, 31, 3);
  for (const auto assembly : {pe::LocalAssembly::kGrid, pe::LocalAssembly::kExactHeight}) {
    const auto p = mb::make_params(pe::LocalAlignment::kPaper, assembly, 3, {1, 1});
    REQUIRE(p.local_branch(31) == pe::LocalBranch::kUnsupported);
    CHECK_THROWS_AS(pe::local_hierarchy(g31, p), std::invalid_argument);
  }
}

TEMPLATE_TEST_CASE("Mode B local: the paper's alignment and grid reproduce paper::local_hierarchy",
                   "[paper_experimental][local]", std::uint64_t, double) {
  using W = TestType;
  for (const std::size_t n : {std::size_t{32}, std::size_t{40}}) {
    const auto g = graph_of<W>(lt::Family::kSparse, lt::Weights::kExponential, n, 11 + n);
    INFO("n=" << n);
    // Mode A: the paper's parameters (base case for k <= 209952, recorded by Mode A).
    {
      Diagnostics want_diag;
      const auto want = sssp::paper::local_hierarchy(g, &want_diag);
      const auto run = pe::local_hierarchy(g, pe::PaperParams{});
      CHECK(run.branch == pe::LocalBranch::kBaseCase);
      mb::require_same_result(run.local, want);
      REQUIRE(run.diag.fallbacks.size() == want_diag.fallbacks.size());
      CHECK(run.diag.fallbacks[0].reason == want_diag.fallbacks[0].reason);
      CHECK(mb::counters_without(run.diag, kEntryPrefix) == mb::counters_without(want_diag, ""));
    }
    // Mode B size test and separator with the paper's Sec 5.4 formulas: the radial branch.
    for (const std::uint64_t b : {std::uint64_t{1}, std::uint64_t{3}, std::uint64_t{6}}) {
      INFO("b=" << b);
      const auto p = mb::make_params(pe::LocalAlignment::kPaper, pe::LocalAssembly::kGrid, b,
                                     pe::SizeThreshold{1, 1});
      Diagnostics want_diag;
      const auto want = sssp::paper::local_hierarchy(g, &want_diag, p.local_hierarchy_params(n));
      const auto run = pe::local_hierarchy(g, p);
      CHECK(run.branch == pe::LocalBranch::kRadial);
      mb::require_same_result(run.local, want);
      CHECK(run.diag.fallbacks.empty());
      CHECK(want_diag.fallbacks.empty());
      CHECK(mb::counters_without(run.diag, kEntryPrefix) == mb::counters_without(want_diag, ""));
    }
  }
}

TEMPLATE_TEST_CASE("Mode B assembly engine with the paper's constants equals Mode A",
                   "[paper_experimental][local]", std::uint64_t, double) {
  using W = TestType;
  const auto g = graph_of<W>(lt::Family::kGrid, lt::Weights::kExponential, 36, 5);
  const auto od = sssp::paper::owner_decomposition(g, 3);
  const auto rs = sssp::paper::local_radial_states(g, od);
  const auto joins = sssp::paper::compute_joins(od, rs);
  Diagnostics want_diag;
  Diagnostics got_diag;
  const auto want = sssp::paper::assemble_local_hierarchy(od, rs, joins, &want_diag);
  const auto got = pe::assemble_local_hierarchy(od, rs, joins, pe::LocalAlignment::kPaper,
                                                pe::LocalAssembly::kGrid, &got_diag);
  mb::require_same_assembly(got, want);
  CHECK(got_diag.counters == want_diag.counters);
  CHECK(got_diag.fallbacks.empty());
}

TEST_CASE("Mode B assembly engine rejects inputs outside its preconditions",
          "[paper_experimental][local]") {
  using W = std::uint64_t;
  const auto g = graph_of<W>(lt::Family::kPath, lt::Weights::kUniform, 31, 2);
  const auto od = sssp::paper::owner_decomposition(g, 2);
  const auto rs = sssp::paper::local_radial_states(g, od);
  const auto joins = sssp::paper::compute_joins(od, rs);
  // Paper formulas: N = 31 < 32 (Sec 5.4 l.1206).
  CHECK_THROWS_AS(pe::assemble_local_hierarchy(od, rs, joins, pe::LocalAlignment::kPaper,
                                               pe::LocalAssembly::kGrid),
                  std::invalid_argument);
  // D-8 and D-11 accept it (N >= 2).
  CHECK_NOTHROW(pe::assemble_local_hierarchy(od, rs, joins, pe::LocalAlignment::kLambda,
                                             pe::LocalAssembly::kGrid));
  CHECK_NOTHROW(pe::assemble_local_hierarchy(od, rs, joins, pe::LocalAlignment::kLambda,
                                             pe::LocalAssembly::kExactHeight));
}

TEST_CASE("D-8: Lambda = 2N * 3hK^2 and the support end 2(3hK^2)^2 N^3 sigma",
          "[paper_experimental][local]") {
  using W = std::uint64_t;
  // N = 2, b = 1: h = 2, K = 16 (least power of two >= 8h), 3hK^2 = 1536.
  CHECK(pe::lambda_scaled<W>(1, 2, 1, 4) == 6144);           // 2 * 2 * 1536
  CHECK(pe::lambda_support_end<W>(1, 2, 1, 4) == 37748736);  // 2 * 1536^2 * 8
  CHECK(pe::lambda_scaled<W>(3, 2, 1, 4) == 18432);
  // N = 40, b = 3: h = 4, K = 32, 3hK^2 = 12288.
  CHECK(pe::lambda_scaled<W>(1, 40, 3, 5) == 983040);                  // 80 * 12288
  CHECK(pe::lambda_support_end<W>(1, 40, 3, 5) == 19327352832000ULL);  // 2 * 12288^2 * 64000
  CHECK(pe::lambda_scaled<double>(0.5, 40, 3, 5) == 491520.0);
  // Saturation at the type's infinity (Q-38), also when h = b + 1 overflows.
  CHECK(pe::lambda_support_end<W>(std::uint64_t{1} << 40, 1U << 20, 1U << 20, 24) ==
        WeightTraits<W>::infinity());
  CHECK(pe::lambda_scaled<W>(1, 2, WeightTraits<W>::infinity(), 4) == WeightTraits<W>::infinity());
  CHECK(pe::lambda_support_end<double>(1e300, 1U << 20, 1U << 20, 24) ==
        WeightTraits<double>::infinity());
}

TEMPLATE_TEST_CASE("D-11: exact heights with links no longer than the height",
                   "[paper_experimental][local]", std::uint64_t, double) {
  using W = TestType;
  const auto g = graph_of<W>(lt::Family::kSparse, lt::Weights::kTies, 12, 9);
  const auto od = sssp::paper::owner_decomposition(g, 2);
  const auto rs = sssp::paper::local_radial_states(g, od);
  const auto joins = sssp::paper::compute_joins(od, rs);
  const auto a = pe::assemble_local_hierarchy(od, rs, joins, pe::LocalAlignment::kLambda,
                                              pe::LocalAssembly::kExactHeight);
  const auto& h = a.result.hierarchy;
  CHECK(h.c0 == 1);
  const auto report = sssp::paper::validate_ordered_hierarchy(h, g);
  INFO(report.summary());
  REQUIRE(report.ok());
  std::size_t internal = 0;
  for (std::size_t v = 0; v < h.nodes.size(); ++v) {
    if (h.is_leaf(static_cast<sssp::paper::NodeId>(v))) continue;
    ++internal;
    const W height = h.nodes[v].height;
    // The height is the weight of a join record (exact, no grid).
    bool is_join_weight = false;
    for (const auto& j : joins.joins) is_join_weight = is_join_weight || j.weight == height;
    CHECK(is_join_weight);
    CHECK(a.node_upper[v] == height);
    CHECK(a.node_upper_inclusive[v] == 1);
    CHECK(a.node_origin[v] == sssp::paper::NodeOrigin::kExactLocal);
    for (const auto& link : h.nodes[v].links) CHECK(!(height < link.length));
  }
  CHECK(internal > 0);
  CHECK(a.grid.empty());
}

TEST_CASE("Mode B local: the parameter bridge keeps D-8 and D-11 away from Mode A",
          "[paper_experimental][local][params]") {
  for (const auto alignment : {pe::LocalAlignment::kPaper, pe::LocalAlignment::kLambda}) {
    for (const auto assembly : {pe::LocalAssembly::kGrid, pe::LocalAssembly::kExactHeight}) {
      const auto p = mb::make_params(alignment, assembly, 3, {1, 1});
      // The size test and the separator are shared by every variant.
      CHECK(p.local_size_params(100) == sssp::paper::LocalHierarchyParams{1, 1, 3});
      if (alignment == pe::LocalAlignment::kPaper && assembly == pe::LocalAssembly::kGrid) {
        CHECK(p.local_hierarchy_params(100) == p.local_size_params(100));
      } else {
        // paper::local_hierarchy would silently run the paper's grid: the bridge refuses and
        // names the Mode B entry point.
        try {
          static_cast<void>(p.local_hierarchy_params(100));
          FAIL("local_hierarchy_params accepted a Mode B assembly variant");
        } catch (const std::logic_error& e) {
          CHECK(std::string(e.what()).find("paper_experimental::local_hierarchy") !=
                std::string::npos);
        }
      }
    }
  }
  CHECK(pe::PaperParams{}.local_size_params(209953) == sssp::paper::kPaperLocalParams);
}
