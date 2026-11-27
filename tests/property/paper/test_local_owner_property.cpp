// Property tests for the owner decomposition of Sec 5.1 (Def 5.1, Lemmas 5.2 and 5.3) on seeded
// graph families, checked against brute-force tree paths written in the test.
//
// Separator sizes: the paper's local b = L_N^4 exceeds N on every graph small enough to test
// exhaustively (one pivot, trivial construction). The tests therefore pass b in
// {1, 2, 3, 4, ceil(log2 N), N, N + 1} directly, as Mode B does (Q-03); these values are a
// deviation from the paper's b and are recorded in docs/PAPER_MAPPING.md.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/primitives/mst.hpp"
#include "sssp/graph.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeId;
using sssp::Vertex;
using sssp::paper::kNoPivot;
using sssp::paper::OwnerDecomposition;
namespace lt = sssp::test::local;

namespace {

std::vector<std::size_t> separator_sizes(std::size_t n) {
  std::size_t lg = 0;
  while ((std::size_t{1} << lg) < n) ++lg;
  std::set<std::size_t> s{1, 2, 3, 4, std::max<std::size_t>(lg, 1), n, n + 1};
  return {s.begin(), s.end()};
}

// Components of the tree minus a vertex set (`removed[v]`), as vertex lists.
template <class W>
std::vector<std::vector<Vertex>> tree_components_without(std::size_t n,
                                                         const std::vector<Edge<W>>& edges,
                                                         const std::vector<EdgeId>& tree,
                                                         const std::vector<bool>& removed) {
  std::vector<std::vector<Vertex>> adj(n);
  for (const EdgeId id : tree) {
    adj[edges[id].u].push_back(edges[id].v);
    adj[edges[id].v].push_back(edges[id].u);
  }
  std::vector<int> comp(n, -1);
  std::vector<std::vector<Vertex>> out;
  for (Vertex s = 0; s < n; ++s) {
    if (removed[s] || comp[s] >= 0) continue;
    out.emplace_back();
    std::vector<Vertex> st{s};
    comp[s] = static_cast<int>(out.size() - 1);
    while (!st.empty()) {
      const Vertex x = st.back();
      st.pop_back();
      out.back().push_back(x);
      for (const Vertex y : adj[x]) {
        if (!removed[y] && comp[y] < 0) {
          comp[y] = comp[s];
          st.push_back(y);
        }
      }
    }
  }
  return out;
}

template <class W>
void check_decomposition(const std::vector<Edge<W>>& edges, std::size_t n, std::size_t b,
                         const OwnerDecomposition<W>& od) {
  // T is a spanning tree of minimum weight.
  REQUIRE(od.tree_edges.size() == n - 1);
  const auto kruskal = sssp::paper::mst_kruskal_ref<W>(n, edges);
  CHECK(lt::total_weight(edges, od.tree_edges) == lt::total_weight(edges, kruskal));
  const auto t = lt::root_tree<W>(n, edges, od.tree_edges);
  for (Vertex v = 1; v < n; ++v) REQUIRE(t.parent[v] != sssp::kNoVertex);

  // Reduced graph: tree edges, and nontree edges with w <= D_T (and b_T <= w by the cut rule).
  std::vector<bool> in_tree(edges.size(), false);
  for (const EdgeId id : od.tree_edges) in_tree[id] = true;
  for (EdgeId i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    const auto path = lt::tree_path(t, e.u, e.v);
    const bool expect = in_tree[i] || (e.u != e.v && !(lt::path_sum(edges, path) < e.weight));
    CHECK((od.reduced[i] != 0) == expect);
    if (e.u != e.v) CHECK_FALSE(e.weight < lt::path_max(edges, path));
  }

  // A0 meets every connected b-vertex subtree: T - A0 has components of < b vertices.
  std::vector<bool> in_a0(n, false);
  std::size_t a0 = 0;
  for (Vertex v = 0; v < n; ++v) {
    in_a0[v] = od.marked[v] != 0;
    a0 += od.marked[v];
  }
  CHECK(a0 >= 1);
  CHECK(a0 <= std::max<std::size_t>(1, n / b));
  for (const auto& c : tree_components_without(n, edges, od.tree_edges, in_a0)) {
    CHECK(c.size() < b);
  }

  // A = A0 + branch vertices of the minimal subtree spanning A0; pivots sorted, indexed.
  std::vector<bool> in_a(n, false);
  REQUIRE(std::is_sorted(od.pivots.begin(), od.pivots.end()));
  for (std::uint32_t i = 0; i < od.pivots.size(); ++i) {
    in_a[od.pivots[i]] = true;
    CHECK(od.pivot_index[od.pivots[i]] == i);
  }
  for (Vertex v = 0; v < n; ++v) {
    if (in_a0[v]) CHECK(in_a[v]);
    if (!in_a[v]) CHECK(od.pivot_index[v] == kNoPivot);
    if (in_a[v] && !in_a0[v]) {
      // v is a branch vertex: at least three components of T - v contain marks.
      std::vector<bool> rm(n, false);
      rm[v] = true;
      int dirs = 0;
      for (const auto& c : tree_components_without(n, edges, od.tree_edges, rm)) {
        dirs += std::any_of(c.begin(), c.end(), [&](Vertex x) { return in_a0[x]; }) ? 1 : 0;
      }
      CHECK(dirs >= 3);
    }
  }
  CHECK(od.pivots.size() <= std::max<std::size_t>(1, 2 * a0 - 1));

  // Components of T - A: < b vertices, at most two distinct pivot neighbours.
  for (const auto& c : tree_components_without(n, edges, od.tree_edges, in_a)) {
    CHECK(c.size() < b);
    std::set<Vertex> nb;
    for (const EdgeId id : od.tree_edges) {
      const auto& e = edges[id];
      const bool u_in = std::find(c.begin(), c.end(), e.u) != c.end();
      const bool v_in = std::find(c.begin(), c.end(), e.v) != c.end();
      if (u_in && in_a[e.v]) nb.insert(e.v);
      if (v_in && in_a[e.u]) nb.insert(e.u);
    }
    CHECK(nb.size() <= 2);
  }

  // Corridors: |A| - 1 tree edges, each a maximum (smallest id on ties) of its corridor.
  REQUIRE(od.corridors.size() == od.pivots.size() - 1);
  std::size_t cut_count = 0;
  for (EdgeId i = 0; i < edges.size(); ++i) cut_count += od.cut[i];
  CHECK(cut_count == od.corridors.size());
  for (const auto& c : od.corridors) {
    REQUIRE(in_tree[c.edge]);
    CHECK(c.weight == edges[c.edge].weight);
    const Vertex p = od.pivots[c.lower];
    const Vertex q = od.pivots[c.upper];
    const auto path = lt::tree_path(t, p, q);
    // The path between the two pivots has no other pivot inside.
    CHECK(std::find(path.begin(), path.end(), c.edge) != path.end());
    for (const EdgeId id : path) {
      CHECK_FALSE(c.weight < edges[id].weight);
      if (edges[id].weight == c.weight) CHECK(c.edge <= id);
    }
    // T_A joins the owners of the cut edge's endpoints.
    const auto a = od.owner[edges[c.edge].u];
    const auto bq = od.owner[edges[c.edge].v];
    CHECK(((a == c.lower && bq == c.upper) || (a == c.upper && bq == c.lower)));
  }

  // Cells: one pivot each; r, R and hops along T[v, p(v)]; R <= h r.
  for (Vertex v = 0; v < n; ++v) {
    const Vertex p = od.pivots[od.owner[v]];
    const auto path = lt::tree_path(t, v, p);
    for (const EdgeId id : path) CHECK(od.cut[id] == 0);
    CHECK(od.radius[v] == lt::path_max(edges, path));
    CHECK(od.owner_distance[v] == lt::path_sum(edges, path));
    CHECK(od.owner_hops[v] == path.size());
    CHECK(od.owner_hops[v] <= od.hop_bound);
    CHECK(od.hop_bound == b + 1);
    W hr{0};
    for (std::size_t i = 0; i < od.hop_bound; ++i) hr += od.radius[v];
    CHECK_FALSE(hr < od.owner_distance[v]);
    if (od.is_pivot(v)) CHECK(od.radius[v] == W{0});
    if (!od.is_pivot(v)) CHECK(W{0} < od.radius[v]);
  }

  // sigma(p, q) = b_T(p, q) = maximum on T_A[p, q]; Lemma 5.3: b_T(v, q) >= r_v.
  for (std::uint32_t p = 0; p < od.pivots.size(); ++p) {
    for (std::uint32_t q = 0; q < od.pivots.size(); ++q) {
      const auto path = lt::tree_path(t, od.pivots[p], od.pivots[q]);
      CHECK(od.sigma(p, q) == lt::path_max(edges, path));
      if (p != q) CHECK(od.corridors[od.corridor_max(p, q)].weight == od.sigma(p, q));
    }
  }
  for (Vertex v = 0; v < n; ++v) {
    for (const Vertex q : od.pivots) {
      CHECK_FALSE(lt::path_max(edges, lt::tree_path(t, v, q)) < od.radius[v]);
    }
  }

  // Lemma 5.2 on the reduced graph. Sums are exact (weights < 2^42, see local_test_support).
  const W nn = static_cast<W>(n);
  for (EdgeId i = 0; i < edges.size(); ++i) {
    if (od.reduced[i] == 0) continue;
    const auto& e = edges[i];
    const W rx = od.radius[e.u];
    const W ry = od.radius[e.v];
    const auto px = od.owner[e.u];
    const auto py = od.owner[e.v];
    const W sg = od.sigma(px, py);
    if (px != py) {
      CHECK_FALSE(e.weight < std::max({rx, ry, sg}));
    } else if (rx != ry) {
      CHECK_FALSE(e.weight < std::max(rx, ry));
    }
    const W h = static_cast<W>(od.hop_bound);
    CHECK_FALSE(h * (rx + ry) + nn * sg < e.weight);
  }
}

}  // namespace

TEMPLATE_TEST_CASE("owner decomposition satisfies Def 5.1 and Lemmas 5.2-5.3 on graph families",
                   "[property]", std::uint64_t, double) {
  using W = TestType;
  const lt::Family families[] = {lt::Family::kPath,   lt::Family::kSpider, lt::Family::kRandomTree,
                                 lt::Family::kSparse, lt::Family::kGrid,   lt::Family::kComplete};
  const lt::Weights dists[] = {lt::Weights::kTies, lt::Weights::kUniform,
                               lt::Weights::kExponential};
  std::uint64_t runs = 0;
  for (const auto f : families) {
    for (const auto wd : dists) {
      for (std::uint64_t seed = 1; seed <= 6; ++seed) {
        const std::size_t n = f == lt::Family::kComplete ? 2 + seed : 3 + 5 * seed;
        const auto edges = lt::make_graph<W>(f, wd, n, seed * 7919 + static_cast<unsigned>(f));
        const auto g = CsrGraph<W>::from_edges(n, edges);
        for (const std::size_t b : separator_sizes(n)) {
          INFO(lt::name(f) << " / " << lt::name(wd) << " n=" << n << " seed=" << seed
                           << " b=" << b);
          const auto od = sssp::paper::owner_decomposition(g, b);
          check_decomposition(edges, n, b, od);
          ++runs;
        }
      }
    }
  }
  CHECK(runs > 500);
}
