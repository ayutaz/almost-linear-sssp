// Differential tests of Algorithm 4 (LocalHierarchy, Lemma 5.10) through Evaluate (Alg 3,
// Sec 4.3) against a Dijkstra written in the test support (no code shared with src/paper):
//   * Lemma 5.10: every geodesic has an equal-length state realization from its native entry,
//     so seeding the native entry of a source s with 0 makes the least emission at every
//     physical vertex equal to d(s, v) (Lemma 4.6 on the seeded hierarchy of Def 4.2);
//   * every emission is at least d(s, pi(x)) and every vertex has a protected state;
//   * the seeded hierarchy validates (Def 4.1 / 4.2);
//   * size bounds of Lemma 5.10 / 5.7 by counters: states <= 4N (at most four copies per
//     vertex, Sec 5.2), internal nodes < states, arcs <= 8m + 2N (one arc per copy, edge and
//     direction, plus at most two switches per vertex).
// Both branches of Alg 4 are exercised: the base case with Mode A parameters (k <= 209952) and
// the radial construction with Mode B parameters (threshold N > L_N, separator b passed in,
// N >= 32 and b + 1 <= N as Sec 5.4 l.1206 requires).
//
// Weights: std::uint64_t and double. The double families use quarter-integers whose sums stay
// below 2^53, so every sum is exact and the comparison is exact equality.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/paper/eval_oracle.hpp"
#include "unit/paper/local/clustered_graphs.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::LocalHierarchyParams;
using sssp::paper::StateId;
namespace lt = sssp::test::local;
namespace tp = sssp::test::paper;

namespace {

struct Seen {
  std::size_t radial = 0;
  std::size_t base = 0;
};

template <class W>
void check_against_dijkstra(const CsrGraph<W>& g, const std::vector<sssp::Edge<W>>& edges,
                            const LocalHierarchyParams& p, Seen& seen) {
  const std::size_t n = g.num_vertices();
  const std::size_t m = g.num_edges();
  Diagnostics diag;
  const auto lh = sssp::paper::local_hierarchy(g, &diag, p);
  const auto& h = lh.hierarchy;
  if (diag.fallbacks.empty()) {
    ++seen.radial;
  } else {
    ++seen.base;
  }
  REQUIRE(lh.native_entry.size() == n);
  REQUIRE(h.num_states() <= 4 * n);
  REQUIRE(h.nodes.size() < 2 * h.num_states());
  REQUIRE(h.graph.arcs.size() <= 8 * m + 2 * n);

  for (const Vertex src : {Vertex{0}, static_cast<Vertex>(n / 2), static_cast<Vertex>(n - 1)}) {
    INFO("source " << src);
    sssp::paper::SeededOrderedHierarchy<W> seeded{h, {{src, W{0}}}, {{lh.native_entry[src], W{0}}}};
    const auto valid = sssp::paper::validate_seeded_ordered_hierarchy(seeded, g);
    INFO(valid.summary());
    REQUIRE(valid.ok());
    const auto ref = tp::compute_eval_reference(seeded, g);
    const auto r = sssp::paper::evaluate(seeded);
    const auto want = lt::dijkstra_naive<W>(n, edges, {{src, W{0}}});
    const auto got = tp::physical_minima<W>(n, h.projection, r.emissions);
    const auto served = tp::served_targets<W>(n, h.projection, ref);
    for (Vertex v = 0; v < n; ++v) {
      INFO("vertex " << v);
      REQUIRE(want[v] != WeightTraits<W>::infinity());
      REQUIRE(ref.target[v] == want[v]);
      REQUIRE(served[v]);
      REQUIRE(got[v] == want[v]);
    }
    for (const auto& [x, val] : r.emissions) REQUIRE(!(val < want[h.projection[x]]));
  }
}

}  // namespace

TEMPLATE_TEST_CASE("Alg 4 + Evaluate from a native entry equals Dijkstra at every vertex", "",
                   std::uint64_t, double) {
  using W = TestType;
  Seen seen;
  for (const auto f : {lt::Family::kPath, lt::Family::kSpider, lt::Family::kRandomTree,
                       lt::Family::kSparse, lt::Family::kGrid, lt::Family::kComplete}) {
    for (const auto wd : {lt::Weights::kTies, lt::Weights::kUniform, lt::Weights::kExponential}) {
      for (const std::size_t n : {32U, 41U}) {
        for (std::uint64_t seed = 1; seed <= 2; ++seed) {
          const auto edges = lt::make_graph<W>(f, wd, n, seed * 131 + n);
          const auto g = CsrGraph<W>::from_edges(n, edges);
          INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " seed=" << seed);
          // Mode A: k <= 209952, the base case (with its FallbackRecord).
          check_against_dijkstra(g, edges, sssp::paper::kPaperLocalParams, seen);
          // Mode B: the radial construction with b in {1, 2, 3, 6, N - 1}.
          for (const std::size_t b :
               {std::size_t{1}, std::size_t{2}, std::size_t{3}, std::size_t{6}, n - 1}) {
            INFO("b=" << b);
            check_against_dijkstra(g, edges, LocalHierarchyParams{1, 1, b}, seen);
          }
        }
      }
    }
  }
  CHECK(seen.radial > 0);
  CHECK(seen.base > 0);
}

TEST_CASE("Alg 4 + Evaluate equals Dijkstra on clustered graphs with unaligned blocks") {
  // Light clusters joined by edges of about 2^54 (Sec 5.4 l.1206-1213: sigma_U > N^10 s_U, so
  // some blocks use a separate grid below a global attachment).
  using W = std::uint64_t;
  Seen seen;
  [[maybe_unused]] std::uint64_t unaligned = 0;
  for (const std::size_t n : {32U, 40U, 64U}) {
    for (const std::size_t b : {2U, 3U, 5U, 8U}) {
      for (std::uint64_t seed = 1; seed <= 3; ++seed) {
        const auto edges = lt::clustered(n, seed * 17 + b);
        const auto g = CsrGraph<W>::from_edges(n, edges);
        INFO("clustered n=" << n << " b=" << b << " seed=" << seed);
        Diagnostics diag;
        static_cast<void>(sssp::paper::local_hierarchy(g, &diag, LocalHierarchyParams{1, 1, b}));
#if SSSP_DIAGNOSTICS
        CHECK(diag.counters["local_assemble_calls"] == 1);
        unaligned += diag.counters["local_assemble_unaligned_blocks"];
#endif
        check_against_dijkstra(g, edges, LocalHierarchyParams{1, 1, b}, seen);
      }
    }
  }
  CHECK(seen.radial == 36);
#if SSSP_DIAGNOSTICS
  CHECK(unaligned > 0);
#endif
}
