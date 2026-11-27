// Differential tests of Algorithm 4 (LocalHierarchy, Sec 5.4, Lemma 5.10) through Evaluate
// (Alg 3, Sec 4.3, Lemma 4.6) with the paper's Mode A constants (kPaperLocalParams), at the
// sizes where the radial branch applies: k > 2 L_k^4, first reached at k = 209953 (L = 18,
// b = 18^4); k = 262143 is the first size with L = 19. k = 209952 is the last base-case size.
//
// For every graph and source s:
//   * the seeded hierarchy (native entry of s offered 0, source (s, 0); Def 4.2) passes the
//     project's validator;
//   * Lemma 5.10 + Lemma 4.6: the least emission at each physical vertex equals d(s, v) from a
//     binary-heap Dijkstra written in the test support (shares no code with src/paper), and no
//     emission is below d(s, pi(x));
//   * fallback policy: the radial branch records no FallbackRecord, the base case records
//     exactly {"local: k <= 2*L_k^4 (paper base case)", k, m}; Evaluate records none.
//
// Weights: std::uint64_t (exact) and quarter-integer doubles with sums below 2^53 (exact), so
// every comparison is exact equality.
//
// Counterexamples found here (minimised, kept in tests/regression/, not fixed by this task):
//   * 20261127-evaluate-mass-overflow: lt::clustered(209953, 107) (uint64 bridges of 2^54 ..
//     2^56, MST total about 2^70) makes Evaluate throw std::overflow_error although every
//     distance fits; the clustered family below therefore uses bridges of 2^40 .. 2^42.
//   * 20261127-evaluate-double-slot-rounding: heavy_clusters<double>(209953, 112, 200) (bridges
//     of 2^200 .. 2^202) makes Evaluate throw std::logic_error; no family with non-representable
//     double sums is kept here.
//
// All cases are tagged [slow]; they are not reduced, because smaller inputs fall below the
// Mode A threshold and take the base case.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "property/paper/local_hierarchy_large_support.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/prng.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::Vertex;
namespace lt = sssp::test::local;
namespace ll = sssp::test::local_large;

namespace {

template <class W>
struct Run {
  sssp::paper::LocalHierarchyResult<W> lh;
  Diagnostics diag;
};

template <class W>
Run<W> build(const CsrGraph<W>& g) {
  Run<W> r;
  r.lh = sssp::paper::local_hierarchy(g, &r.diag);
  return r;
}

void check_radial_diag(const Diagnostics& diag) {
  CHECK(diag.fallbacks.empty());
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters.at("local_hierarchy_radial") == 1);
  CHECK(diag.counters.count("local_hierarchy_base_case") == 0);
  CHECK(diag.counters.at("local_assemble_calls") == 1);
#endif
}

/// Seeds the native entry of `src`, validates Def 4.1 / 4.2 and runs Evaluate. Returns the
/// least emission per vertex; checks emissions against `lower` (nothing below the target).
template <class W>
std::vector<W> evaluate_from(const CsrGraph<W>& g, const sssp::paper::LocalHierarchyResult<W>& lh,
                             Vertex src, const std::vector<W>& lower) {
  const auto& h = lh.hierarchy;
  const std::size_t n = g.num_vertices();
  sssp::paper::SeededOrderedHierarchy<W> seeded{h, {{src, W{0}}}, {{lh.native_entry[src], W{0}}}};
  const auto valid = sssp::paper::validate_seeded_ordered_hierarchy(seeded, g);
  INFO(valid.summary());
  CHECK(valid.ok());
  const auto r = sssp::paper::evaluate(seeded);
  CHECK(r.diag.fallbacks.empty());
  std::vector<W> best(n, ll::infinity<W>());
  std::size_t below = 0;
  for (const auto& [x, val] : r.emissions) {
    const Vertex v = h.projection[x];
    if (val < lower[v]) ++below;
    if (val < best[v]) best[v] = val;
  }
  CHECK(below == 0);
  return best;
}

/// Exact families: the least emission equals the Dijkstra distance at every vertex.
template <class W>
void check_exact(std::size_t n, const std::vector<Edge<W>>& edges, bool expect_radial) {
  const auto g = CsrGraph<W>::from_edges(n, edges);
  REQUIRE(sssp::paper::local_construction_applies(n, sssp::paper::kPaperLocalParams) ==
          expect_radial);
  const auto run = build(g);
  if (expect_radial) {
    check_radial_diag(run.diag);
  } else {
    REQUIRE(run.diag.fallbacks.size() == 1);
    CHECK(run.diag.fallbacks[0].reason == "local: k <= 2*L_k^4 (paper base case)");
    CHECK(run.diag.fallbacks[0].vertices == n);
    CHECK(run.diag.fallbacks[0].edges == edges.size());
#if SSSP_DIAGNOSTICS
    CHECK(run.diag.counters.at("local_hierarchy_base_case") == 1);
    CHECK(run.diag.counters.count("local_hierarchy_radial") == 0);
#endif
  }
  REQUIRE(run.lh.native_entry.size() == n);
  const auto phys = ll::undirected(n, edges);
  for (const Vertex src : {Vertex{0}, static_cast<Vertex>(n / 2), static_cast<Vertex>(n - 1)}) {
    INFO("source " << src);
    const auto want = ll::dijkstra(phys, {{src, W{0}}});
    const auto got = evaluate_from(g, run.lh, src, want);
    std::size_t mismatches = 0;
    Vertex first_bad = sssp::kNoVertex;
    for (Vertex v = 0; v < n; ++v) {
      REQUIRE(want[v] < ll::infinity<W>());
      if (got[v] != want[v]) {
        if (mismatches == 0) first_bad = v;
        ++mismatches;
      }
    }
    INFO("first mismatching vertex " << first_bad);
    CHECK(mismatches == 0);
  }
}

/// Light clusters of 3..8 vertices (weights 1..3: a random spanning tree plus one chord) joined
/// by a random spanning tree of bridges and about (clusters / 2) extra bridges of length
/// 2^e, 2^(e+1) or 2^(e+2). With uint64 and e = 40 every distance and the MST total stay far
/// below 2^64.
template <class W>
std::vector<Edge<W>> heavy_clusters(std::size_t n, std::uint64_t seed, int e) {
  sssp::test::SplitMix64 rng(seed);
  std::vector<Edge<W>> out;
  const auto light = [&](std::size_t u, std::size_t v) {
    out.push_back(
        {static_cast<Vertex>(u), static_cast<Vertex>(v), static_cast<W>(rng.between(1, 3))});
  };
  const auto heavy = [&](Vertex u, Vertex v) {
    const int x = e + static_cast<int>(rng.between(0, 2));
    if constexpr (std::is_floating_point_v<W>) {
      out.push_back({u, v, std::ldexp(W{1}, x)});
    } else {
      out.push_back({u, v, W{1} << x});
    }
  };
  std::vector<std::size_t> start{0};
  while (start.back() < n) {
    start.push_back(std::min<std::size_t>(n, start.back() + 3 + rng.below(6)));
  }
  const std::size_t k = start.size() - 1;
  for (std::size_t c = 0; c < k; ++c) {
    for (std::size_t v = start[c] + 1; v < start[c + 1]; ++v) {
      light(start[c] + rng.below(v - start[c]), v);
    }
    if (start[c + 1] - start[c] >= 3) light(start[c], start[c + 1] - 1);
  }
  const auto pick = [&](std::size_t c) {
    return static_cast<Vertex>(start[c] + rng.below(start[c + 1] - start[c]));
  };
  for (std::size_t c = 1; c < k; ++c) heavy(pick(rng.below(c)), pick(c));
  for (std::size_t i = 0; i < k / 2; ++i) heavy(pick(rng.below(k)), pick(rng.below(k)));
  return out;
}

}  // namespace

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209953: path, exponential weights",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_exact<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kPath, lt::Weights::kExponential, n, 101), true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209953: spider, tied weights", "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_exact<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kSpider, lt::Weights::kTies, n, 102), true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209953: random tree, uniform weights",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_exact<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kRandomTree, lt::Weights::kUniform, n, 103),
      true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209953: sparse multigraph, tied weights",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_exact<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kSparse, lt::Weights::kTies, n, 104), true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209953: sparse multigraph, exponential weights",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_exact<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kSparse, lt::Weights::kExponential, n, 105),
      true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209953: grid, uniform weights", "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_exact<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kGrid, lt::Weights::kUniform, n, 106), true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209953: clusters joined by 2^40 .. 2^42 bridges",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_exact<std::uint64_t>(n, heavy_clusters<std::uint64_t>(n, 107, 40), true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209953: grid, quarter-integer doubles with ties",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_exact<double>(n, lt::make_graph<double>(lt::Family::kGrid, lt::Weights::kTies, n, 108),
                      true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209953: random tree, quarter-integer doubles",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_exact<double>(
      n, lt::make_graph<double>(lt::Family::kRandomTree, lt::Weights::kUniform, n, 109), true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 262143 (first L = 19): grid, exponential weights",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstL19;
  check_exact<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kGrid, lt::Weights::kExponential, n, 110), true);
}

TEST_CASE("Alg 4 Mode A + Evaluate at k = 209952: base case with its FallbackRecord",
          "[slow][local-large]") {
  const std::size_t n = ll::kLastBaseCase;
  check_exact<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kSparse, lt::Weights::kUniform, n, 111), false);
}
