// Differential test of Algorithm 8 l.2 + Lemma 8.5 (arXiv 2611.34911v1, Sec 8.3, txt
// l.2145-2185): contract zero-edge components, solve the contracted positive graph with a
// baseline Dijkstra, restore the original distances and tree, and compare with Bellman–Ford on
// the original graph. Weights are zero_mixed, so zero edges, zero cycles and (for random regular
// graphs with kKeep) parallel edges and self-loops all occur.
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

#include "paper/outer/contract_restore.hpp"
#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/verify.hpp"

using sssp::CsrGraph;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::outer::contract_zero_edges;
using sssp::paper::outer::restore_predecessor_tree;

namespace {

struct Family {
  std::string type;
  std::size_t n;
  sssp::gen::MultiEdgePolicy multi = sssp::gen::MultiEdgePolicy::kSwitchToSimple;
};

std::vector<Family> families() {
  return {
      {"er", 80},        {"er", 300},
      {"grid", 64},      {"regular", 60, sssp::gen::MultiEdgePolicy::kKeep},
      {"ba", 70},        {"tree", 50},
      {"chung_lu", 120},
  };
}

/// Number of zero-edge classes among the vertices that Bellman–Ford reaches (naive labels).
template <class W>
std::size_t zero_classes_reached(const CsrGraph<W>& g, const std::vector<W>& bf_dist) {
  const std::size_t n = g.num_vertices();
  std::vector<std::size_t> p(n);
  std::iota(p.begin(), p.end(), std::size_t{0});
  auto find = [&](std::size_t x) {
    while (p[x] != x) x = p[x];
    return x;
  };
  for (const auto& e : g.edges()) {
    if (e.weight == W{0}) p[find(e.u)] = find(e.v);
  }
  std::size_t classes = 0;
  for (Vertex v = 0; v < n; ++v) {
    if (find(v) == v && bf_dist[v] != sssp::WeightTraits<W>::infinity()) ++classes;
  }
  return classes;
}

}  // namespace

TEMPLATE_TEST_CASE("contract, Dijkstra on the contracted graph, restore equals Bellman-Ford",
                   "[differential][paper][outer][contract]", std::uint64_t, double) {
  using W = TestType;
  std::uint64_t runs = 0;
  std::uint64_t contracted_away = 0;
  for (const auto& fam : families()) {
    for (const double p0 : {0.0, 0.3, 0.7, 0.95, 1.0}) {
      for (std::uint64_t seed = 1; seed <= 6; ++seed) {
        sssp::gen::GraphSpec spec;
        spec.type = fam.type;
        spec.n = fam.n;
        spec.degree = 3;
        spec.multi = fam.multi;
        spec.avg_degree = 2.5;
        spec.options.seed = seed;
        spec.options.weights.dist = sssp::gen::WeightDist::kZeroMixed;
        spec.options.weights.lo = 1.0;
        spec.options.weights.hi = 20.0;
        spec.options.weights.zero_probability = p0;
        const auto g = sssp::to_csr(sssp::gen::generate<W>(spec));
        const std::size_t n = g.num_vertices();
        for (const Vertex s : {Vertex{0}, static_cast<Vertex>(n / 2)}) {
          INFO(fam.type << " n=" << n << " p0=" << p0 << " seed=" << seed << " source=" << s);
          sssp::Diagnostics diag;
          const auto c = contract_zero_edges(g, s, &diag);
          // Alg 8 l.2: only positive non-loop edges remain.
          for (const auto& e : c.graph.edges()) {
            REQUIRE(e.u != e.v);
            REQUIRE(e.weight > W{0});
          }
          const auto r = sssp::baseline::dijkstra_binary_heap(c.graph, c.source);
          const auto res = restore_predecessor_tree(g, s, c, r.dist);
          const auto bf = sssp::baseline::bellman_ford(g, s);
          // N = number of zero classes in the source component.
          REQUIRE(c.graph.num_vertices() == zero_classes_reached(g, bf.dist));
          for (Vertex v = 0; v < n; ++v) {
            REQUIRE((c.vertex_map[v] == kNoVertex) ==
                    (bf.dist[v] == sssp::WeightTraits<W>::infinity()));
          }
          const auto tol = sssp::Tolerance::for_graph(g);
          const auto rep = sssp::verify_sssp(g, s, res, &bf.dist, tol);
          INFO(rep.summary());
          REQUIRE(rep.ok);
          if constexpr (sssp::WeightTraits<W>::is_exact) REQUIRE(res.dist == bf.dist);
          // Lemma 8.5: the restored tree is verified against the original graph and weights.
          const auto self = sssp::verify_sssp(g, s, res, nullptr, tol);
          REQUIRE(self.ok);
          contracted_away += diag.counters.at("contract_removed_vertices");
          ++runs;
        }
      }
    }
  }
  REQUIRE(runs == families().size() * 5 * 6 * 2);
  REQUIRE(contracted_away > 0);
}
