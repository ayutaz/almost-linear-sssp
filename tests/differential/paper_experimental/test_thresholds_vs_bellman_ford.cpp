// Mode B with lowered thresholds (task P5-EXP-THRESHOLDS; D-7 with the limits of Q-03):
// the named Q-03 settings run the outer construction (Alg 8 l.8-15, Sec 8.3) and, inside its
// atoms (Alg 7 l.6, Lemma 5.15), the local radial construction (Alg 4 l.6-12, Sec 5.1-5.4) on
// feasible sizes, and the distances are exact.
//
// Evidence that the constructions ran and were not replaced by a direct search comes from
// Diagnostics: the Mode B branch counters (written also with SSSP_DIAGNOSTICS = 0), the
// FallbackRecords (only Alg 4's recorded base case on small atoms; never the outer base case on
// the outer branch) and, with SSSP_DIAGNOSTICS, the counters of the shared components.
// Oracles: Bellman-Ford on the original graph (computed once per (g, s)) and verify_sssp
// (root, every arc feasible, tight tree arcs, parents reach s, unreachable = infinity).
// Inputs: clusters of light edges (weights 1..4) joined by heavy edges, so that atoms reach the
// local radial size; and library generator families. Integer weights only, also as integral
// (dyadic) doubles: Q-06 / Q-55 plan A; uint64 weights stay far below the D-20 range.
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "experimental/paper_experimental/params.hpp"
#include "experimental/paper_experimental/sssp.hpp"
#include "experimental/paper_experimental/thresholds.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "paper/outer/base_case.hpp"
#include "paper/outer/shortest_paths.hpp"
#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"
#include "support/prng.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::SsspResult;
using sssp::Vertex;
namespace pe = sssp::paper_experimental;
namespace po = sssp::paper::outer;

namespace {

/// `clusters` clusters of `size` vertices: a random recursive tree plus size/2 chords, weights
/// 1..4, the clusters joined in a chain by one edge of weight `heavy` each. Connected, no loops.
template <class W>
CsrGraph<W> clustered(std::size_t clusters, std::size_t size, std::uint64_t heavy,
                      std::uint64_t seed) {
  sssp::test::SplitMix64 rng(seed);
  std::vector<Edge<W>> edges;
  for (std::size_t c = 0; c < clusters; ++c) {
    const std::size_t base = c * size;
    for (std::size_t i = 1; i < size; ++i) {
      edges.push_back({static_cast<Vertex>(base + rng.below(i)), static_cast<Vertex>(base + i),
                       static_cast<W>(rng.between(1, 4))});
    }
    for (std::size_t i = 0; i < size / 2; ++i) {
      const auto u = static_cast<Vertex>(base + rng.below(size));
      const auto v = static_cast<Vertex>(base + rng.below(size));
      if (u != v) edges.push_back({u, v, static_cast<W>(rng.between(1, 4))});
    }
    if (c > 0) {
      edges.push_back({static_cast<Vertex>(base - 1 - rng.below(size)),
                       static_cast<Vertex>(base + rng.below(size)), static_cast<W>(heavy)});
    }
  }
  return CsrGraph<W>::from_edges(clusters * size, edges);
}

/// A library family with integer weights in [1, hi], converted to W (integral for double).
template <class W>
CsrGraph<W> family(std::string_view type, std::size_t n, std::uint64_t seed, double hi) {
  sssp::gen::GraphSpec spec;
  spec.type = std::string(type);
  spec.n = n;
  spec.avg_degree = 3.0;
  spec.degree = 3;
  spec.options.seed = seed;
  spec.options.weights.dist = sssp::gen::WeightDist::kUniformInt;
  spec.options.weights.lo = 1.0;
  spec.options.weights.hi = hi;
  spec.options.connected = true;
  const auto base = sssp::to_csr(sssp::gen::generate<std::uint64_t>(spec));
  std::vector<Edge<W>> edges;
  for (const auto& e : base.edges()) edges.push_back({e.u, e.v, static_cast<W>(e.weight)});
  return CsrGraph<W>::from_edges(base.num_vertices(), edges);
}

std::uint64_t counter(const sssp::Diagnostics& d, const std::string& key) {
  const auto it = d.counters.find(key);
  return it == d.counters.end() ? 0 : it->second;
}

/// One run of `setting` on (g, s): exact against Bellman-Ford, a valid tree, the outer branch
/// taken, and every Alg 4 call accounted for. Returns the number of radial Alg 4 calls.
template <class W>
std::uint64_t run_outer(const pe::ThresholdSetting& setting, const CsrGraph<W>& g, Vertex s,
                        const SsspResult<W>& bf, std::uint64_t seed) {
  INFO(setting.name << ", N = " << g.num_vertices() << ", s = " << s << ", seed = " << seed);
  REQUIRE(setting.params.outer_applies(g.num_vertices()));
  const auto r = pe::paper_sssp(g, s, setting.params, seed);
  REQUIRE(r.dist == bf.dist);
  const auto report = sssp::verify_sssp(g, s, r, &bf.dist);
  INFO(report.summary());
  REQUIRE(report.ok);
  // Alg 8 l.8-15 ran; the outer base case did not.
  REQUIRE(counter(r.diag, "mode_b.branch.outer_construction") == 1);
  REQUIRE(counter(r.diag, "mode_b.branch.outer_base_case") == 0);
  // Every FallbackRecord is Alg 4's recorded base case on an atom, one per base-case call.
  const std::string local =
      sssp::paper::local_base_case_reason(setting.params.local_size_params(g.num_vertices()));
  for (const auto& fb : r.diag.fallbacks) REQUIRE(fb.reason == local);
  const std::uint64_t calls = counter(r.diag, "mode_b.local.calls");
  const std::uint64_t radial = counter(r.diag, "mode_b.local.radial");
  const std::uint64_t base = counter(r.diag, "mode_b.local.base_case");
  REQUIRE(radial + base == calls);
  REQUIRE(base == r.diag.fallbacks.size());
#if SSSP_DIAGNOSTICS
  // The shared components' own counts agree (observed, not derived).
  REQUIRE(counter(r.diag, "outer.constructions") == 1);
  REQUIRE(counter(r.diag, "local_hierarchy_calls") == calls);
  REQUIRE(counter(r.diag, "local_hierarchy_radial") == radial);
  REQUIRE(counter(r.diag, "local_hierarchy_base_case") == base);
#endif
  return radial;
}

const pe::ThresholdSetting& setting_named(std::string_view name) {
  static const auto all = pe::q03_threshold_settings();
  for (const auto& s : all) {
    if (s.name == name) return s;
  }
  FAIL("no Q-03 setting named " << name);
  return all.front();
}

}  // namespace

TEMPLATE_TEST_CASE("Q-03 settings: outer and local radial constructions on clustered graphs",
                   "[paper_experimental][thresholds][differential]", std::uint64_t, double) {
  using W = TestType;
  // 10 clusters of 40 vertices (N = 400), heavy chain edges 2^20.
  for (const auto& setting : pe::q03_threshold_settings()) {
    if (!setting.params.outer_applies(400)) continue;  // the scaled setting: separate case
    std::uint64_t radial_total = 0;
    for (const std::uint64_t gseed : {1, 2, 3}) {
      const auto g = clustered<W>(10, 40, std::uint64_t{1} << 20, gseed);
      for (const Vertex s : {Vertex{0}, Vertex{399}}) {
        const auto bf = sssp::baseline::bellman_ford(g, s);
        for (const std::uint64_t seed : {std::uint64_t{0}, std::uint64_t{7}}) {
          radial_total += run_outer(setting, g, s, bf, seed);
        }
      }
    }
    INFO(setting.name << ": radial Alg 4 calls " << radial_total);
    // The local radial construction is expected iff its least radial size fits in a cluster
    // (the settings with a lowered local threshold use an outer b >= 64, so atoms can be that
    // large); with the paper's local threshold (209953) it never runs here.
    const auto least = setting.least_local_radial_size;
    if (least && *least <= 40) {
      CHECK(radial_total > 0);
    } else {
      CHECK(radial_total == 0);
    }
  }
}

TEMPLATE_TEST_CASE("Q-03 settings: generator families stay exact on the outer branch",
                   "[paper_experimental][thresholds][differential]", std::uint64_t, double) {
  using W = TestType;
  std::map<std::string, std::uint64_t> radial;
  for (const std::string_view type : {"er", "grid", "tree", "geometric", "regular", "ba"}) {
    for (const std::uint64_t gseed : {11, 12}) {
      const auto g = family<W>(type, 500, gseed, 1e6);
      const Vertex s{0};
      const auto bf = sssp::baseline::bellman_ford(g, s);
      for (const auto& setting : pe::q03_threshold_settings()) {
        if (!setting.params.outer_applies(g.num_vertices())) continue;
        INFO(type << " gseed " << gseed);
        radial[setting.name] += run_outer(setting, g, s, bf, gseed);
      }
    }
  }
  for (const auto& [name, count] : radial) {
    INFO(name << ": radial Alg 4 calls " << count);
    const auto least = setting_named(name).least_local_radial_size;
    if (!least || *least > 500) CHECK(count == 0);
  }
}

TEMPLATE_TEST_CASE(
    "Q-03 scaled setting: base case below the outer threshold, both "
    "constructions above it",
    "[paper_experimental][thresholds][differential]", std::uint64_t, double) {
  using W = TestType;
  const auto& setting = setting_named("q03-scaled-L3-L2");
  // Below the least outer size: the recorded outer base case, nothing else.
  {
    const auto g = clustered<W>(10, 40, std::uint64_t{1} << 20, 5);
    REQUIRE_FALSE(setting.params.outer_applies(g.num_vertices()));
    const auto bf = sssp::baseline::bellman_ford(g, Vertex{0});
    const auto r = pe::paper_sssp(g, Vertex{0}, setting.params);
    REQUIRE(r.dist == bf.dist);
    REQUIRE(sssp::verify_sssp(g, Vertex{0}, r, &bf.dist).ok);
    CHECK(counter(r.diag, "mode_b.branch.outer_base_case") == 1);
    CHECK(counter(r.diag, "mode_b.local.calls") == 0);
    REQUIRE(r.diag.fallbacks.size() == 1);
    CHECK(r.diag.fallbacks[0].reason == po::kOuterBaseCaseScanReason);
  }
  // Above it (N = 1600 > 1331, inside the upward-closed range): outer and local radial.
  std::uint64_t radial_total = 0;
  for (const std::uint64_t gseed : {6, 7}) {
    const auto g = clustered<W>(40, 40, std::uint64_t{1} << 20, gseed);
    const auto bf = sssp::baseline::bellman_ford(g, Vertex{0});
    radial_total += run_outer(setting, g, Vertex{0}, bf, gseed);
  }
  CHECK(radial_total > 0);
}
