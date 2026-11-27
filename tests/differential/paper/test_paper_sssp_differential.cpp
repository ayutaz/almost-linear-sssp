// Mode A entry point of Algorithm 8 (arXiv 2611.34911v1, Sec 8.3, Alg 8 l.2-7, Theorem 1,
// Lemma 8.5): paper_sssp contracts zero-edge components (l.2), evaluates the outer threshold
// N > L_N^64 on the contracted size N (l.3-4, Sec 2 txt l.284-286), runs the recorded base case
// (l.5; scan Dijkstra for "paper-faithful", binary heap for "paper-faithful-d1" = D-1, ADR 0008)
// and restores the original distances and tree (l.6, Lemma 8.5).
//
// Oracles: Bellman–Ford and baseline binary-heap Dijkstra on the ORIGINAL graph, verify_sssp for
// the tree, and an independent union-find count of the contracted size for the FallbackRecord.
// The fallback tests carry no SSSP_DIAGNOSTICS guard: they also run in the release preset
// (diagnostics off), where a missing record would be a silent fallback.
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/outer/contract_restore.hpp"
#include "paper/sssp.hpp"
#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/paper.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"

using sssp::CsrGraph;
using sssp::kNoVertex;
using sssp::SsspResult;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::PaperOptions;
using sssp::paper::PaperVariant;

namespace {

constexpr PaperVariant kBoth[] = {PaperVariant::kFaithful, PaperVariant::kFaithfulD1};

std::string expected_reason(PaperVariant v) {
  return v == PaperVariant::kFaithful
             ? "outer: N <= L_N^64 (paper base case, scan Dijkstra)"
             : "outer: N <= L_N^64 (paper base case, binary-heap Dijkstra, deviation D-1)";
}

/// Size of the contracted source component (Alg 8 l.2-3), computed independently of the
/// implementation: union-find over zero edges, restricted to the vertices Bellman–Ford reaches;
/// edges are those joining two different reached classes (loops discarded, parallels kept).
struct ContractedSize {
  std::uint64_t vertices = 0;
  std::uint64_t edges = 0;
};

template <class W>
ContractedSize contracted_size(const CsrGraph<W>& g, const std::vector<W>& bf_dist) {
  const std::size_t n = g.num_vertices();
  std::vector<std::size_t> p(n);
  std::iota(p.begin(), p.end(), std::size_t{0});
  auto find = [&](std::size_t x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  };
  for (const auto& e : g.edges()) {
    if (e.weight == W{0}) p[find(e.u)] = find(e.v);
  }
  const W inf = WeightTraits<W>::infinity();
  ContractedSize sz;
  for (std::size_t v = 0; v < n; ++v) {
    if (find(v) == v && bf_dist[v] != inf) ++sz.vertices;
  }
  for (const auto& e : g.edges()) {
    if (bf_dist[e.u] != inf && find(e.u) != find(e.v)) ++sz.edges;
  }
  return sz;
}

template <class W>
void require_one_record(const SsspResult<W>& r, PaperVariant v, const ContractedSize& sz) {
  REQUIRE(r.diag.fallbacks.size() == 1);
  CHECK(r.diag.fallbacks[0].reason == expected_reason(v));
  CHECK(r.diag.fallbacks[0].vertices == sz.vertices);
  CHECK(r.diag.fallbacks[0].edges == sz.edges);
  CHECK(r.diag.counters.at("outer_base_cases") == 1);
  CHECK(r.diag.counters.at("contract_retained_vertices") == sz.vertices);
  CHECK(r.diag.counters.at("contract_retained_edges") == sz.edges);
}

}  // namespace

TEST_CASE("paper algorithm names", "[paper][sssp]") {
  CHECK(sssp::paper::algorithm_name(PaperVariant::kFaithful) == "paper-faithful");
  CHECK(sssp::paper::algorithm_name(PaperVariant::kFaithfulD1) == "paper-faithful-d1");
  REQUIRE(sssp::paper::kPaperVariants.size() == 2);
  for (const PaperVariant v : sssp::paper::kPaperVariants) {
    const auto parsed = sssp::paper::parse_algorithm_name(sssp::paper::algorithm_name(v));
    REQUIRE(parsed.has_value());
    CHECK(*parsed == v);
  }
  CHECK_FALSE(sssp::paper::parse_algorithm_name("paper").has_value());
  CHECK_FALSE(sssp::paper::parse_algorithm_name("dijkstra").has_value());
  CHECK_FALSE(sssp::paper::parse_algorithm_name("paper-faithful-d2").has_value());
  CHECK(PaperOptions{}.variant == PaperVariant::kFaithful);  // the faithful variant is default
}

TEST_CASE("Mode A never takes the outer construction for N < 2^64", "[paper][sssp][threshold]") {
  // Alg 8 l.4, Sec 2 l.286: N > L_N^64 is false whenever L_N >= 2 and N < 2^64.
  for (std::uint64_t n = 1; n <= 100000; ++n) {
    REQUIRE_FALSE(sssp::paper::takes_outer_construction(n));
  }
  for (unsigned b = 17; b < 64; ++b) {
    const std::uint64_t p = std::uint64_t{1} << b;
    REQUIRE_FALSE(sssp::paper::takes_outer_construction(p - 2));
    REQUIRE_FALSE(sssp::paper::takes_outer_construction(p - 1));
    REQUIRE_FALSE(sssp::paper::takes_outer_construction(p));
    REQUIRE_FALSE(sssp::paper::takes_outer_construction(p + 1));
  }
  REQUIRE_FALSE(sssp::paper::takes_outer_construction(std::numeric_limits<std::uint64_t>::max()));
}

TEST_CASE("the missing outer construction throws instead of being replaced by a search",
          "[paper][sssp][outer]") {
  const auto g = CsrGraph<std::uint64_t>::from_edges(3, {{0, 1, 2}, {1, 2, 0}});
  const auto c = sssp::paper::outer::contract_zero_edges(g, 0);
  sssp::Diagnostics diag;
  try {
    (void)sssp::paper::outer_construction(c, diag);
    FAIL("outer_construction returned");
  } catch (const std::logic_error& e) {
    CHECK(std::string(e.what()) == "outer construction not implemented");
  }
  CHECK(diag.fallbacks.empty());  // nothing was solved, nothing is reported as solved
}

TEMPLATE_TEST_CASE("paper_sssp: hand example with zero classes, loops and an unreachable part",
                   "[paper][sssp][fallback]", std::uint64_t, double) {
  using W = TestType;
  // e0 0-1 (0), e1 1-2 (3), e2 2-3 (0), e3 3-0 (5), e4 1-1 (2), e5 0-1 (4), e6 4-5 (1), 6 alone.
  // From 0: classes {0,1}, {2,3}; contracted edges e1, e3 (e0, e2, e4, e5 become loops, e6 is
  // unreachable) -> N = 2, m' = 2. d = {0, 0, 3, 3, inf, inf, inf}.
  const auto g = CsrGraph<W>::from_edges(7, {{0, 1, W{0}},
                                             {1, 2, W{3}},
                                             {2, 3, W{0}},
                                             {3, 0, W{5}},
                                             {1, 1, W{2}},
                                             {0, 1, W{4}},
                                             {4, 5, W{1}}});
  const W inf = WeightTraits<W>::infinity();
  INFO("SSSP_DIAGNOSTICS=" << SSSP_DIAGNOSTICS);
  for (const PaperVariant v : kBoth) {
    const auto r = sssp::paper::paper_sssp(g, 0, PaperOptions{v});
    CHECK(r.dist == std::vector<W>{W{0}, W{0}, W{3}, W{3}, inf, inf, inf});
    CHECK(r.parent[0] == kNoVertex);
    CHECK(r.parent[4] == kNoVertex);
    CHECK(r.parent[6] == kNoVertex);
    require_one_record(r, v, ContractedSize{2, 2});
    REQUIRE(sssp::verify_sssp(g, 0, r).ok);
  }
  // Default options are the faithful variant.
  const auto r = sssp::paper::paper_sssp(g, 0);
  require_one_record(r, PaperVariant::kFaithful, ContractedSize{2, 2});
  // From the isolated vertex: N = 1, no contracted edge.
  for (const PaperVariant v : kBoth) {
    const auto iso = sssp::paper::paper_sssp(g, 6, PaperOptions{v});
    CHECK(iso.dist == std::vector<W>{inf, inf, inf, inf, inf, inf, W{0}});
    require_one_record(iso, v, ContractedSize{1, 0});
  }
}

TEMPLATE_TEST_CASE("paper_sssp: all-zero graph contracts to one vertex", "[paper][sssp][fallback]",
                   std::uint64_t, double) {
  using W = TestType;
  const auto g =
      CsrGraph<W>::from_edges(4, {{0, 1, W{0}}, {1, 2, W{0}}, {2, 0, W{0}}, {2, 3, W{0}}});
  for (const PaperVariant v : kBoth) {
    const auto r = sssp::paper::paper_sssp(g, 2, PaperOptions{v});
    CHECK(r.dist == std::vector<W>(4, W{0}));
    require_one_record(r, v, ContractedSize{1, 0});
    REQUIRE(sssp::verify_sssp(g, 2, r).ok);
  }
}

TEST_CASE("paper_sssp rejects an invalid source", "[paper][sssp]") {
  const auto g = CsrGraph<std::uint64_t>::from_edges(2, {{0, 1, 1}});
  for (const PaperVariant v : kBoth) {
    REQUIRE_THROWS_AS(sssp::paper::paper_sssp(g, 2, PaperOptions{v}), std::out_of_range);
  }
  const auto empty = CsrGraph<double>::from_edges(0, {});
  REQUIRE_THROWS_AS(sssp::paper::paper_sssp(empty, 0), std::out_of_range);
}

TEMPLATE_TEST_CASE("paper_sssp equals Bellman-Ford and Dijkstra on every generator family",
                   "[differential][paper][sssp]", std::uint64_t, double) {
  using W = TestType;
  struct WeightCase {
    sssp::gen::WeightDist dist;
    double zero_probability;
  };
  const WeightCase weight_cases[] = {
      {WeightTraits<W>::is_exact ? sssp::gen::WeightDist::kUniformInt
                                 : sssp::gen::WeightDist::kUniformReal,
       0.0},
      {sssp::gen::WeightDist::kZeroMixed, 0.3},
      {sssp::gen::WeightDist::kZeroMixed, 0.9},
  };
  constexpr std::uint64_t kSeeds = 5;
  std::size_t runs = 0;
  std::uint64_t merged_or_unreached = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    for (const auto& wc : weight_cases) {
      for (std::uint64_t seed = 1; seed <= kSeeds; ++seed) {
        sssp::gen::GraphSpec spec;
        spec.type = std::string(info.name);
        spec.n = 100;
        spec.avg_degree = 2.5;
        spec.degree = 3;
        spec.multi = sssp::gen::MultiEdgePolicy::kKeep;  // parallel edges and loops for regular
        spec.options.seed = seed;
        spec.options.weights.dist = wc.dist;
        spec.options.weights.lo = 1.0;
        spec.options.weights.hi = 40.0;
        spec.options.weights.zero_probability = wc.zero_probability;
        const auto g = sssp::to_csr(sssp::gen::generate<W>(spec));
        const std::size_t n = g.num_vertices();
        const auto tol = sssp::Tolerance::for_graph(g);
        for (const Vertex s : {Vertex{0}, static_cast<Vertex>(n / 2)}) {
          const auto bf = sssp::baseline::bellman_ford(g, s);
          const auto dij = sssp::baseline::dijkstra_binary_heap(g, s);
          const ContractedSize sz = contracted_size(g, bf.dist);
          for (const PaperVariant v : kBoth) {
            INFO(info.name << " " << sssp::gen::to_string(wc.dist) << " p0=" << wc.zero_probability
                           << " seed=" << seed << " s=" << s << " "
                           << sssp::paper::algorithm_name(v));
            const auto r = sssp::paper::paper_sssp(g, s, PaperOptions{v});
            const auto vs_bf = sssp::verify_sssp(g, s, r, &bf.dist, tol);
            INFO(vs_bf.summary());
            REQUIRE(vs_bf.ok);
            REQUIRE(sssp::verify_sssp(g, s, r, &dij.dist, tol).ok);
            if constexpr (WeightTraits<W>::is_exact) {
              REQUIRE(r.dist == bf.dist);
              REQUIRE(r.dist == dij.dist);
            }
            require_one_record(r, v, sz);
            // Fallback share: the base case handled the whole contracted graph.
            CHECK(r.diag.fallbacks[0].vertices <= n);
            CHECK(r.diag.fallbacks[0].edges <= g.num_edges());
            merged_or_unreached += n - sz.vertices;
            ++runs;
          }
        }
      }
    }
  }
  REQUIRE(runs == sssp::gen::generator_catalog().size() * 3 * kSeeds * 2 * 2);
  REQUIRE(merged_or_unreached > 0);  // contraction actually changed N on some inputs
}

TEST_CASE("paper-faithful-d1 at N = 90000 with zero edges equals baseline Dijkstra",
          "[differential][paper][sssp]") {
  sssp::gen::GraphSpec spec;
  spec.type = "grid";
  spec.rows = 300;
  spec.cols = 300;
  spec.n = 300 * 300;
  spec.options.seed = 11;
  spec.options.weights.dist = sssp::gen::WeightDist::kZeroMixed;
  spec.options.weights.lo = 1.0;
  spec.options.weights.hi = 1000.0;
  spec.options.weights.zero_probability = 0.2;
  const auto g = sssp::to_csr(sssp::gen::generate<std::uint64_t>(spec));
  const Vertex s = 12345;
  const auto ref = sssp::baseline::dijkstra_binary_heap(g, s);
  const auto r = sssp::paper::paper_sssp(g, s, PaperOptions{PaperVariant::kFaithfulD1});
  REQUIRE(r.dist == ref.dist);
  REQUIRE(sssp::verify_sssp(g, s, r, &ref.dist).ok);
  REQUIRE(r.diag.fallbacks.size() == 1);
  CHECK(r.diag.fallbacks[0].reason == expected_reason(PaperVariant::kFaithfulD1));
  CHECK(r.diag.fallbacks[0].vertices < g.num_vertices());  // zero edges were contracted
}
