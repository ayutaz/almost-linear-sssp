// Differential tests for the outer state graph of arXiv 2611.34911v1, Sec 8.1 (Lemma 8.1,
// txt l.2002-2026), against baseline binary-heap Dijkstra, Bellman-Ford, a naive exact Dijkstra
// written in the test support, and brute-force enumeration of all simple shortest paths.
//
// Checked statements:
//   * "The ordinary graph realizes every geodesic" (Lemma 8.1; Lemma 5.7 with K = 16): from the
//     native inward copy of every source s, the least state distance over the ordinary copies of
//     v equals d(s, v), both in the ordinary graph and with the special states added (l.1994-1995:
//     the added arcs project to physical walks and preserve all lower bounds, so no state is
//     below the physical distance of its vertex).
//   * Deletion (l.1967-1968): no shortest path uses a deleted or omitted edge (reference trees
//     of Dijkstra and Bellman-Ford, and every simple shortest path for N <= 7).
//   * Swing (Lemma 8.1, l.2004-2005): every maximal radial run of those shortest paths, split at
//     any minimum of J, has J <= 16 x the running minimum on the prefix and 16 J >= the running
//     maximum on the suffix.
// Inputs: the generators of outer_states_test_support.hpp (Mode A and Mode B separator sizes,
// exact and NameRoutes names) and clusters of 100..110 vertices with LocalHierarchyParams
// {2, 2, b}, so that OwnerDistances runs Alg 4's radial branch on its atoms (Q-57).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "paper/local/local_hierarchy.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/outer/outer_states.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/prng.hpp"
#include "unit/paper/local/local_test_support.hpp"
#include "unit/paper/outer/outer_states_test_support.hpp"

using sssp::EdgeId;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::EdgeClass;
using sssp::paper::StateId;
using sssp::paper::outer::build_outer_states;
using sssp::paper::outer::OuterArcKind;
using sssp::paper::outer::OuterStates;
namespace to = sssp::test::outer;
namespace lt = sssp::test::local;
using to::Exact;
using to::exact;
using to::kExactInf;

namespace {

struct Coverage {
  std::uint64_t cases = 0;
  std::uint64_t sources = 0;
  std::uint64_t runs = 0;
  std::uint64_t runs_with_swing = 0;  ///< runs whose J varies by a factor >= 2
  std::uint64_t geodesics = 0;        ///< enumerated simple shortest paths
  std::uint64_t special_minima = 0;   ///< vertices whose least state distance is at a special state
  std::uint64_t radial_atoms = 0;
};

/// Lemma 8.1: swing <= 16 on the maximal radial runs of the path given by vertices and edges.
template <class W>
void check_path(const OuterStates<W>& os, const std::vector<Vertex>& verts,
                const std::vector<EdgeId>& path, Coverage& cov) {
  const auto& cls = os.states.classification.edge_class;
  const auto jj = [&](Vertex v) { return exact(os.states.coordinate[v]); };
  std::size_t i = 0;
  while (i < path.size()) {
    CHECK(cls[path[i]] != EdgeClass::kOmitted);  // deleted / omitted edges are never geodesic
    if (cls[path[i]] != EdgeClass::kRadial) {
      ++i;
      continue;
    }
    std::size_t k = i;
    while (k < path.size() && cls[path[k]] == EdgeClass::kRadial) ++k;
    // Run: vertices verts[i..k].
    Exact lo = kExactInf;
    Exact hi = 0;
    for (std::size_t t = i; t <= k; ++t) {
      lo = std::min(lo, jj(verts[t]));
      hi = std::max(hi, jj(verts[t]));
    }
    ++cov.runs;
    if (hi >= 2 * lo) ++cov.runs_with_swing;
    for (std::size_t mpos = i; mpos <= k; ++mpos) {
      if (jj(verts[mpos]) != lo) continue;
      INFO("run of " << (k - i) << " edges split at position " << (mpos - i));
      Exact run_min = kExactInf;
      for (std::size_t t = i; t <= mpos; ++t) {  // inward prefix
        run_min = std::min(run_min, jj(verts[t]));
        CHECK(jj(verts[t]) <= 16 * run_min);
      }
      Exact run_max = 0;
      for (std::size_t t = mpos; t <= k; ++t) {  // outward suffix
        run_max = std::max(run_max, jj(verts[t]));
        CHECK(16 * jj(verts[t]) >= run_max);
      }
    }
    i = k;
  }
}

/// Every simple path from s (DFS) whose length equals d(s, end): checked by check_path.
template <class W>
void enumerate_geodesics(const to::OuterCase<W>& c, const OuterStates<W>& os, Vertex s,
                         const std::vector<Exact>& d, Coverage& cov) {
  std::vector<Vertex> verts{s};
  std::vector<EdgeId> path;
  std::vector<std::uint8_t> on(c.n, 0);
  on[s] = 1;
  const auto rec = [&](auto&& self, Vertex u, Exact len) -> void {
    if (len == d[u] && !path.empty()) {
      ++cov.geodesics;
      check_path(os, verts, path, cov);
    }
    for (const auto& arc : c.g.arcs(u)) {
      if (on[arc.to] != 0) continue;
      const Exact next = len + exact(arc.weight);
      if (next > d[arc.to]) continue;  // a prefix of a shortest path is shortest
      on[arc.to] = 1;
      verts.push_back(arc.to);
      path.push_back(arc.edge);
      self(self, arc.to, next);
      path.pop_back();
      verts.pop_back();
      on[arc.to] = 0;
    }
  };
  rec(rec, s, 0);
}

/// Reference tree path s -> v from parent pointers.
template <class W>
void tree_path(const sssp::SsspResult<W>& r, Vertex v, std::vector<Vertex>& verts,
               std::vector<EdgeId>& edges) {
  verts.assign(1, v);
  edges.clear();
  while (r.parent[v] != kNoVertex) {
    edges.push_back(r.parent_edge[v]);
    v = r.parent[v];
    verts.push_back(v);
  }
  std::ranges::reverse(verts);
  std::ranges::reverse(edges);
}

template <class W>
void check_case(const to::OuterCase<W>& c, std::size_t max_sources, Coverage& cov) {
  const auto os = build_outer_states(c.g, c.od, c.cls, c.aux, c.jr.scaled_distance);
  const auto& st = os.states;
  const std::size_t s_count = st.num_states();
  std::vector<to::XArc> all, ordinary;
  for (std::size_t a = 0; a < st.graph.arcs.size(); ++a) {
    const auto& arc = st.graph.arcs[a];
    const to::XArc x{arc.tail, arc.head, exact(arc.length)};
    all.push_back(x);
    if (os.arc_kind[a] == OuterArcKind::kOrdinary) ordinary.push_back(x);
  }
  for (const auto r : c.jr.radial_atom) cov.radial_atoms += r;
  const std::size_t step = std::max<std::size_t>(1, c.n / std::max<std::size_t>(1, max_sources));
  std::vector<Vertex> verts;
  std::vector<EdgeId> edges;
  for (Vertex s = 0; s < c.n; s += static_cast<Vertex>(step)) {
    INFO("source " << s);
    const auto d = to::physical_distances(c.n, c.edges, s);
    const auto dj = sssp::baseline::dijkstra_binary_heap(c.g, s);
    const auto bf = sssp::baseline::bellman_ford(c.g, s);
    for (Vertex v = 0; v < c.n; ++v) {
      REQUIRE(d[v] != kExactInf);  // the inputs are connected
      REQUIRE(exact(dj.dist[v]) == d[v]);
      REQUIRE(exact(bf.dist[v]) == d[v]);
    }
    const StateId src = st.native_inward[s];
    const auto da = to::state_distances(s_count, all, src);
    const auto dor = to::state_distances(s_count, ordinary, src);
    std::vector<Exact> best_all(c.n, kExactInf), best_ord(c.n, kExactInf),
        best_ordinary_copy(c.n, kExactInf);
    for (StateId x = 0; x < s_count; ++x) {
      const Vertex v = st.projection[x];
      CHECK(da[x] >= d[v]);  // lower bounds (l.1994)
      best_all[v] = std::min(best_all[v], da[x]);
      if (x < os.num_ordinary_states) {
        best_ord[v] = std::min(best_ord[v], dor[x]);
        best_ordinary_copy[v] = std::min(best_ordinary_copy[v], da[x]);
      }
    }
    for (Vertex v = 0; v < c.n; ++v) {
      INFO("vertex " << v);
      CHECK(best_ord[v] == d[v]);            // the ordinary graph realizes every geodesic
      CHECK(best_ordinary_copy[v] == d[v]);  // ... also with the special states present
      CHECK(best_all[v] == d[v]);
    }
    for (const auto* r : {&dj, &bf}) {
      for (Vertex v = 0; v < c.n; ++v) {
        tree_path(*r, v, verts, edges);
        REQUIRE(verts.front() == s);
        check_path(os, verts, edges, cov);
      }
    }
    if (c.n <= 7) enumerate_geodesics(c, os, s, d, cov);
    ++cov.sources;
  }
  ++cov.cases;
}

/// Clusters of 100..110 vertices (light trees and chords of weight 1..3, heavier chords
/// 40..80) below vertex 0 by spokes 1000 * 2^U[0, 3] or below an earlier cluster by 100..300,
/// plus root shortcuts and cross edges (radial_clusters of the Alg 7 tests, with larger
/// clusters): with one pivot every atom is a union of whole clusters, k >= 100 > 2 L_k^2 = 98,
/// so LocalHierarchyParams{2, 2, b} takes Alg 4's radial branch (Q-57).
template <class W>
std::vector<sssp::Edge<W>> big_clusters(std::size_t clusters, std::uint64_t seed, std::size_t& n) {
  sssp::test::SplitMix64 rng(seed);
  std::vector<sssp::Edge<W>> e;
  std::vector<std::vector<Vertex>> members(clusters);
  Vertex next = 1;
  for (std::size_t c = 0; c < clusters; ++c) {
    const std::size_t size = rng.between(100, 110);
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

}  // namespace

TEMPLATE_TEST_CASE("Lemma 8.1: outer states realize every geodesic; swing <= 16 on shortest paths",
                   "[differential][outer_states]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  to::for_each_case<W>({1, 2, 4, 6, 7, 12, 20}, 2,
                       [&](to::OuterCase<W>& c) { check_case(c, 20, cov); });
  WARN("cases " << cov.cases << ", sources " << cov.sources << ", radial runs " << cov.runs
                << " (swinging " << cov.runs_with_swing << "), simple geodesics " << cov.geodesics);
  CHECK(cov.cases > 0);
  CHECK(cov.runs_with_swing > 0);
  CHECK(cov.geodesics > 0);
}

TEMPLATE_TEST_CASE("Lemma 8.1: outer states on J from radial atoms (Mode B atoms {2, 2, b})",
                   "[differential][outer_states]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  for (std::size_t clusters = 1; clusters <= 2; ++clusters) {
    for (std::uint64_t seed = 1; seed <= 2; ++seed) {
      std::size_t n = 0;
      const std::uint64_t gseed = seed * 104729 + clusters;
      const auto edges = big_clusters<W>(clusters, gseed, n);
      INFO("big-clusters k=" << clusters << " seed=" << gseed << " n=" << n);
      // b = n + 1: one pivot (vertex 0), so the atoms are whole clusters.
      auto c = to::make_outer_case<W>(n, edges, n + 1, to::Names::kExactP1, gseed,
                                      sssp::paper::LocalHierarchyParams{2, 2, 4});
      check_case(c, 3, cov);
    }
  }
  WARN("cases " << cov.cases << ", radial atoms " << cov.radial_atoms << ", radial runs "
                << cov.runs);
  CHECK(cov.cases == 4);
  CHECK(cov.radial_atoms > 0);
}
