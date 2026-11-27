// Differential tests of the Mode B local hierarchy (task P5-EXP-LOCAL) through Evaluate (Alg 3)
// against a Dijkstra written in the test support (no code shared with src/):
//   * D-8 (Lambda alignment and support, research note 02 §4.3) with the grid assembly, and
//     D-11 (exact-height assembly, Q-39), from a native entry seeded with 0: the least emission
//     at every physical vertex equals d(s, v) (Lemma 5.10 + Lemma 4.6), every emission is at
//     least the distance of its image, and the seeded hierarchy validates (Def 4.1 / 4.2);
//   * the grid and the exact-height versions give identical per-vertex distances (P3-LOCAL-
//     ASSEMBLE criterion carried by this task), for D-8 on every N >= 2 and for the paper's
//     Sec 5.4 formulas where they apply (N >= 32, b + 1 <= N);
//   * separators b in {1, 2, 3, 4, ceil(log2 N)} (acceptance criterion) plus N and N + 1;
//   * the D-8 premise "N >= 2, any b >= 1" is exercised from N = 2 upwards: no instance may
//     throw (a violated scale relation of Sec 5.4 throws std::logic_error).
// Weights: std::uint64_t and quarter-integer double (all sums exact, comparisons are equality).

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "experimental/paper_experimental/local_hierarchy.hpp"
#include "experimental/paper_experimental/params.hpp"
#include "sssp/graph.hpp"
#include "support/paper_experimental/local_mode_b_support.hpp"
#include "unit/paper/local/clustered_graphs.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Vertex;
namespace pe = sssp::paper_experimental;
namespace lt = sssp::test::local;
namespace mb = sssp::test::mode_b;

namespace {

struct Variant {
  pe::LocalAlignment alignment;
  pe::LocalAssembly assembly;
  pe::SizeThreshold threshold;
};

/// Runs every variant on g; checks each against Dijkstra and all variants against each other.
template <class W>
void check_variants(const CsrGraph<W>& g, const std::vector<sssp::Edge<W>>& edges, std::uint64_t b,
                    const std::vector<Variant>& variants, std::size_t& runs) {
  const std::size_t n = g.num_vertices();
  std::vector<std::vector<std::vector<W>>> dist;  // [variant][source] -> per-vertex distances
  for (const auto& v : variants) {
    const auto params = mb::make_params(v.alignment, v.assembly, b, v.threshold);
    INFO("alignment " << static_cast<int>(v.alignment) << " assembly "
                      << static_cast<int>(v.assembly));
    REQUIRE(params.local_branch(n) == pe::LocalBranch::kRadial);
    const auto run = pe::local_hierarchy(g, params);
    REQUIRE(run.branch == pe::LocalBranch::kRadial);
    REQUIRE(run.diag.fallbacks.empty());
    REQUIRE(run.local.hierarchy.c0 == (v.assembly == pe::LocalAssembly::kGrid ? 2U : 1U));
    ++runs;
    auto& d = dist.emplace_back();
    for (const Vertex src : {Vertex{0}, static_cast<Vertex>(n / 2), static_cast<Vertex>(n - 1)}) {
      d.push_back(mb::evaluate_and_check(run.local, g, edges, src));
    }
  }
  for (std::size_t i = 1; i < dist.size(); ++i) REQUIRE(dist[i] == dist[0]);
}

}  // namespace

TEMPLATE_TEST_CASE("Mode B local (D-8, D-11) + Evaluate equals Dijkstra; grid equals exact",
                   "[paper_experimental][local]", std::uint64_t, double) {
  using W = TestType;
  std::size_t runs = 0;
  std::size_t paper_runs = 0;
  for (const auto f : {lt::Family::kPath, lt::Family::kSpider, lt::Family::kRandomTree,
                       lt::Family::kSparse, lt::Family::kGrid, lt::Family::kComplete}) {
    for (const auto wd : {lt::Weights::kTies, lt::Weights::kUniform, lt::Weights::kExponential}) {
      for (const std::size_t n : {2U, 3U, 5U, 9U, 17U, 33U}) {
        const auto edges = lt::make_graph<W>(f, wd, n, 977 * n + 13);
        const auto g = CsrGraph<W>::from_edges(n, edges);
        for (const std::uint64_t b : mb::separators(n)) {
          INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " b=" << b);
          std::vector<Variant> variants{
              {pe::LocalAlignment::kLambda, pe::LocalAssembly::kGrid, {0, 0}},
              {pe::LocalAlignment::kLambda, pe::LocalAssembly::kExactHeight, {0, 0}},
          };
          if (n >= 32 && b < n) {  // Sec 5.4 l.1206 as written (Q-03 conservative variant)
            variants.push_back({pe::LocalAlignment::kPaper, pe::LocalAssembly::kGrid, {1, 1}});
            variants.push_back(
                {pe::LocalAlignment::kPaper, pe::LocalAssembly::kExactHeight, {1, 1}});
            ++paper_runs;
          }
          check_variants(g, edges, b, variants, runs);
        }
      }
    }
  }
  CHECK(runs > 0);
  CHECK(paper_runs > 0);
}

TEST_CASE("Mode B local (D-8, D-11) on clustered graphs with unaligned blocks",
          "[paper_experimental][local]") {
  // Light clusters joined by edges of about 2^54: sigma_U far above Lambda s_U, so D-8 has
  // unaligned blocks (separate grids below a global attachment) also for N < 32.
  using W = std::uint64_t;
  std::size_t runs = 0;
  [[maybe_unused]] std::uint64_t unaligned = 0;
  [[maybe_unused]] std::uint64_t unaligned_small = 0;  // N < 32: outside the paper's formulas
  [[maybe_unused]] std::uint64_t aligned = 0;
  for (const std::size_t n : {6U, 12U, 20U, 32U, 40U}) {
    for (const std::uint64_t b : {1U, 2U, 3U, 5U}) {
      for (std::uint64_t seed = 1; seed <= 3; ++seed) {
        const auto edges = lt::clustered(n, seed * 31 + b);
        const auto g = CsrGraph<W>::from_edges(n, edges);
        INFO("clustered n=" << n << " b=" << b << " seed=" << seed);
        const auto params =
            mb::make_params(pe::LocalAlignment::kLambda, pe::LocalAssembly::kGrid, b);
        const auto run = pe::local_hierarchy(g, params);
#if SSSP_DIAGNOSTICS
        CHECK(run.diag.counters.at("local_assemble_calls") == 1);
        unaligned += run.diag.counters.at("local_assemble_unaligned_blocks");
        if (n < 32) unaligned_small += run.diag.counters.at("local_assemble_unaligned_blocks");
        aligned += run.diag.counters.at("local_assemble_aligned_blocks");
#endif
        static_cast<void>(run);
        check_variants(g, edges, b,
                       {{pe::LocalAlignment::kLambda, pe::LocalAssembly::kGrid, {0, 0}},
                        {pe::LocalAlignment::kLambda, pe::LocalAssembly::kExactHeight, {0, 0}}},
                       runs);
      }
    }
  }
  CHECK(runs == 120);
#if SSSP_DIAGNOSTICS
  CHECK(unaligned > 0);
  CHECK(unaligned_small > 0);
  CHECK(aligned > 0);
#endif
}
