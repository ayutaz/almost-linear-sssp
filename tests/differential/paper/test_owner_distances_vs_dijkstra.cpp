// Differential test of Algorithm 7 (OwnerDistances, Sec 7.4, Lemma 7.9) of arXiv 2611.34911v1:
// J must equal the auxiliary distances d_{w'}(v, a) computed by two independent oracles on G',
// the binary-heap Dijkstra and Bellman-Ford of src/baseline (exact equality: integer weights, or
// quarter-integer doubles whose sums stay below 2^53 so that every sum is exact, Q-06 / Q-55
// plan A).
//
// Inputs: G' of Sec 7.1-7.2 (owner_decomposition -> names -> classify_outer_edges ->
// build_aux_graph) on seeded graph families, separator sizes b in {1, 2, 3, 4, ceil(log2 N), N,
// N + 1} (Mode B style, Q-03; b >= N is the paper's single pivot, Q-48) and two name sources
// (exact names with P = 1 and NameRoutes with the paper's parameters). With the paper's Alg 4
// parameters every atom takes the base case (Q-04); the clustered inputs at the end use the
// Mode B parameters {1, 2, b} of LocalHierarchyParams so that atoms of 37..50 vertices take the
// radial branch, which the test requires to happen.
//
// Beyond equality, J is checked as a potential on G' (J_a = 0, J_v <= J_u + w' on every edge,
// and a tight edge into every nonroot vertex), the shortest-path-tree verification of the
// correctness rules without parent pointers (Alg 7 outputs distances only, l.14).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "paper/local/local_hierarchy.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/owner_distances.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::LocalHierarchyParams;
using sssp::paper::owner::AuxGraph;
using sssp::paper::pivots::PivotRoutingInstance;
using sssp::paper::pivots::RouteNames;
namespace lt = sssp::test::local;
namespace ow = sssp::paper::owner;

namespace {

std::vector<std::size_t> separator_sizes(std::size_t n) {
  std::size_t lg = 0;
  while ((std::size_t{1} << lg) < n) ++lg;
  std::set<std::size_t> s{1, 2, 3, 4, std::max<std::size_t>(lg, 1), n, n + 1};
  return {s.begin(), s.end()};
}

/// Exact names with P = 1 (test_aux_metric_property.cpp): the least instance weight lambda <= a_e
/// with lambda >= d_{G_P}(x, y).
template <class W>
RouteNames<W> exact_names(const PivotRoutingInstance<W>& inst) {
  const auto dist = sssp::test::floyd_warshall<W>(inst.t, inst.edges);
  std::vector<W> weights;
  for (const auto& e : inst.edges) weights.push_back(e.weight);
  std::ranges::sort(weights);
  std::vector<W> chosen;
  for (const auto& e : inst.edges) {
    W best = e.weight;
    for (const W a : weights) {
      if (a > e.weight) break;
      if (a >= dist[e.u][e.v]) {
        best = a;
        break;
      }
    }
    chosen.push_back(best);
  }
  RouteNames<W> names;
  names.p_log2 = 0;
  names.dictionary = chosen;
  std::ranges::sort(names.dictionary);
  names.dictionary.erase(std::unique(names.dictionary.begin(), names.dictionary.end()),
                         names.dictionary.end());
  for (const W c : chosen) {
    const auto it = std::ranges::lower_bound(names.dictionary, c);
    names.name.push_back(static_cast<std::uint32_t>(it - names.dictionary.begin()));
  }
  return names;
}

template <class W>
AuxGraph<W> make_aux_graph(const std::vector<Edge<W>>& edges, std::size_t n, std::size_t b,
                           bool exact_p1, std::uint64_t seed) {
  const auto g = CsrGraph<W>::from_edges(n, edges);
  const auto od = sssp::paper::owner_decomposition(g, b);
  const auto pe = ow::projected_interowner_edges(g, od);
  ow::OuterEdgeClassification<W> cls;
  if (exact_p1) {
    cls = ow::classify_outer_edges(g, od, pe, exact_names(pe.instance));
  } else {
    sssp::paper::model::WordRng rng(seed);
    cls = ow::classify_interowner_edges(g, od, rng);
  }
  return ow::build_aux_graph(g, od, cls);
}

/// Light clusters under heavy spokes (test_calendars_property.cpp): multi-vertex atoms with
/// defining chords, several graph pieces.
template <class W>
std::vector<Edge<W>> heavy_spokes(std::size_t n, std::uint64_t seed) {
  sssp::test::SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  if (n <= 1) return e;
  const std::size_t clusters = 1 + rng.below(3);
  std::vector<std::vector<Vertex>> members(clusters);
  for (Vertex v = 1; v < n; ++v) members[rng.below(clusters)].push_back(v);
  for (const auto& mem : members) {
    if (mem.empty()) continue;
    e.push_back({0, mem[0], lt::to_weight<W>(1000ULL << rng.between(0, 10), rng)});
    for (std::size_t i = 1; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(i)], mem[i], lt::to_weight<W>(rng.between(1, 3), rng)});
    }
    for (std::size_t i = 0; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(mem.size())], mem[rng.below(mem.size())],
                   lt::to_weight<W>(rng.between(1, 3), rng)});
    }
  }
  return e;
}

/// Clusters of 37..50 vertices (light trees of weight 1..3, chords of weight 1..3 and 40..80)
/// hanging from vertex 0 by spokes of 1000 * 2^U[0, 3], or below an earlier cluster by a medium
/// edge of weight 100..300 (same root bottleneck: several atoms in one family), plus root
/// shortcuts and cross edges of weight 100..5000. Every vertex has root bottleneck >= 1000, so
/// every cluster edge is defining (w' <= 3 < 2q, q >= r / 64) and every atom is a union of
/// whole clusters (>= 37 vertices): with LocalHierarchyParams{1, 2, b} each atom takes the
/// radial branch of Alg 4 (k > L_k^2, k >= 32).
template <class W>
std::vector<Edge<W>> radial_clusters(std::size_t clusters, std::uint64_t seed, std::size_t& n) {
  sssp::test::SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  std::vector<std::vector<Vertex>> members(clusters);
  Vertex next = 1;
  for (std::size_t c = 0; c < clusters; ++c) {
    const std::size_t size = rng.between(37, 50);
    for (std::size_t i = 0; i < size; ++i) members[c].push_back(next++);
  }
  n = next;
  for (std::size_t c = 0; c < clusters; ++c) {
    const auto& mem = members[c];
    if (c == 0 || rng.below(2) == 0) {
      e.push_back({0, mem[0], lt::to_weight<W>(1000ULL << rng.between(0, 3), rng)});
    } else {
      const auto& up = members[rng.below(c)];
      e.push_back({up[rng.below(up.size())], mem[0], lt::to_weight<W>(rng.between(100, 300), rng)});
    }
    for (std::size_t i = 1; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(i)], mem[i], lt::to_weight<W>(rng.between(1, 3), rng)});
    }
    for (std::size_t i = 0; i < mem.size() / 3; ++i) {
      e.push_back({mem[rng.below(mem.size())], mem[rng.below(mem.size())],
                   lt::to_weight<W>(rng.between(1, 3), rng)});
      e.push_back({mem[rng.below(mem.size())], mem[rng.below(mem.size())],
                   lt::to_weight<W>(rng.between(40, 80), rng)});
    }
  }
  for (std::size_t i = 0; i < clusters; ++i) {
    const auto& mem = members[rng.below(clusters)];
    e.push_back({0, mem[rng.below(mem.size())], lt::to_weight<W>(rng.between(1000, 9000), rng)});
    const auto& a = members[rng.below(clusters)];
    const auto& b = members[rng.below(clusters)];
    e.push_back({a[rng.below(a.size())], b[rng.below(b.size())],
                 lt::to_weight<W>(rng.between(100, 5000), rng)});
  }
  return e;
}

struct Coverage {
  std::uint64_t cases = 0;
  std::uint64_t atoms = 0;
  std::uint64_t multi_atoms = 0;
  std::uint64_t multi_atom_families = 0;
  std::uint64_t radial_atoms = 0;
  std::uint64_t j_below_r = 0;  ///< vertices with J_v < R_v (J is not just the tree distance)
};

/// J equals both oracles exactly and is a tight potential on G'.
template <class W>
void check_against_oracles(const AuxGraph<W>& aux, const LocalHierarchyParams& lp, Coverage& cov) {
  Diagnostics diag;
  ow::OwnerDistancesOptions opt;
  opt.local = lp;
  const auto res = ow::owner_distances(aux, &diag, opt);
  const std::size_t na = aux.num_vertices();
  const auto dij = sssp::baseline::dijkstra_binary_heap(aux.graph, ow::kAuxRoot);
  const auto bf = sssp::baseline::bellman_ford(aux.graph, ow::kAuxRoot);
  REQUIRE(res.distance.size() == na);
  REQUIRE(res.distance == dij.dist);
  REQUIRE(res.distance == bf.dist);

  // Potential check (independent of the oracles).
  REQUIRE(res.distance[ow::kAuxRoot] == W{0});
  std::vector<std::uint8_t> tight(na, 0);
  for (const auto& ed : aux.graph.edges()) {
    const W ju = res.distance[ed.u];
    const W jv = res.distance[ed.v];
    REQUIRE(jv <= ju + ed.weight);
    REQUIRE(ju <= jv + ed.weight);
    if (jv == ju + ed.weight) tight[ed.v] = 1;
    if (ju == jv + ed.weight) tight[ed.u] = 1;
  }
  for (Vertex v = 1; v < na; ++v) REQUIRE(tight[v] == 1);

  // Fallbacks: exactly the atoms that took Alg 4's base case, each recorded.
  std::size_t base = 0;
  for (std::size_t a = 0; a < res.atoms.num_atoms(); ++a) {
    const bool radial = sssp::paper::local_construction_applies(res.atoms.members(a).size(), lp);
    CHECK(res.radial_atom[a] == (radial ? 1 : 0));
    if (!radial) ++base;
    if (res.atoms.members(a).size() > 1) ++cov.multi_atoms;
  }
  CHECK(diag.fallbacks.size() == base);
  for (const auto& f : diag.fallbacks) CHECK(f.reason == sssp::paper::local_base_case_reason(lp));

  ++cov.cases;
  cov.atoms += res.atoms.num_atoms();
  cov.radial_atoms += res.atoms.num_atoms() - base;
  for (std::size_t f = 0; f < res.families.num_families(); ++f) {
    if (res.families.offsets[f + 1] - res.families.offsets[f] > 1) ++cov.multi_atom_families;
  }
  for (Vertex v = 1; v < na; ++v) {
    if (res.distance[v] < aux.owner_distance[v]) ++cov.j_below_r;
  }
}

}  // namespace

TEMPLATE_TEST_CASE("Alg 7: J equals Dijkstra and Bellman-Ford on G' (Mode A atom parameters)",
                   "[differential][owner_distances]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  const lt::Family families[] = {lt::Family::kPath,   lt::Family::kSpider, lt::Family::kRandomTree,
                                 lt::Family::kSparse, lt::Family::kGrid,   lt::Family::kComplete};
  const lt::Weights weights[] = {lt::Weights::kTies, lt::Weights::kUniform,
                                 lt::Weights::kExponential};
  const std::size_t sizes[] = {1, 2, 4, 7, 12, 20, 33};
  for (const auto f : families) {
    for (const auto wd : weights) {
      for (const std::size_t n : sizes) {
        for (std::uint64_t seed = 1; seed <= 2; ++seed) {
          const std::uint64_t gseed = seed * 1000003 + n * 31 + static_cast<std::uint64_t>(f) * 7 +
                                      static_cast<std::uint64_t>(wd);
          const auto edges = lt::make_graph<W>(f, wd, n, gseed);
          for (const std::size_t b : separator_sizes(n)) {
            for (const bool exact_p1 : {true, false}) {
              INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " seed=" << gseed
                               << " b=" << b << (exact_p1 ? " exact-P1" : " name-routes-paper"));
              const auto aux = make_aux_graph(edges, n, b, exact_p1, gseed + b);
              check_against_oracles(aux, sssp::paper::kPaperLocalParams, cov);
            }
          }
        }
      }
    }
  }
  for (const std::size_t n : sizes) {
    for (std::uint64_t seed = 1; seed <= 4; ++seed) {
      const std::uint64_t gseed = seed * 7919 + n;
      const auto edges = heavy_spokes<W>(n, gseed);
      for (const std::size_t b : separator_sizes(n)) {
        for (const bool exact_p1 : {true, false}) {
          INFO("heavy-spokes n=" << n << " seed=" << gseed << " b=" << b
                                 << (exact_p1 ? " exact-P1" : " name-routes-paper"));
          const auto aux = make_aux_graph(edges, n, b, exact_p1, gseed + b);
          check_against_oracles(aux, sssp::paper::kPaperLocalParams, cov);
        }
      }
    }
  }
  INFO("cases " << cov.cases << ", atoms " << cov.atoms << ", multi-vertex atoms "
                << cov.multi_atoms << ", multi-atom families " << cov.multi_atom_families
                << ", J < R " << cov.j_below_r);
  CHECK(cov.multi_atoms > 0);
  CHECK(cov.multi_atom_families > 0);
  CHECK(cov.j_below_r > 0);
  CHECK(cov.radial_atoms == 0);  // Q-04: the paper's Alg 4 threshold is never met here
}

TEMPLATE_TEST_CASE("Alg 7: J equals the oracles with radial atoms (Mode B atom parameters, Q-04)",
                   "[differential][owner_distances]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  for (std::size_t clusters = 1; clusters <= 4; ++clusters) {
    for (std::uint64_t seed = 1; seed <= 3; ++seed) {
      std::size_t n = 0;
      const std::uint64_t gseed = seed * 104729 + clusters;
      const auto edges = radial_clusters<W>(clusters, gseed, n);
      for (const std::size_t local_b : {std::size_t{2}, std::size_t{5}}) {
        for (const bool exact_p1 : {true, false}) {
          INFO("radial-clusters k=" << clusters << " seed=" << gseed << " local b=" << local_b
                                    << (exact_p1 ? " exact-P1" : " name-routes-paper"));
          // b = n + 1: one pivot (vertex 0), so the atoms are whole clusters.
          const auto aux = make_aux_graph(edges, n, n + 1, exact_p1, gseed);
          check_against_oracles(aux, LocalHierarchyParams{1, 2, local_b}, cov);
        }
      }
    }
  }
  INFO("cases " << cov.cases << ", atoms " << cov.atoms << ", radial atoms " << cov.radial_atoms
                << ", multi-atom families " << cov.multi_atom_families << ", J < R "
                << cov.j_below_r);
  CHECK(cov.radial_atoms == cov.atoms);
  CHECK(cov.multi_atom_families > 0);
  CHECK(cov.j_below_r > 0);
}
