// Sec 3.1, Lemma 3.1: depth-first search, breadth-first search, connected components and
// topological sorting on k vertices and l arcs. References are written independently in this
// file (recursive DFS, transitive closure by Floyd-Warshall style iteration).
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "paper/primitives/search.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::DiArc;
using sssp::paper::Digraph;

namespace {

using Adj = std::vector<std::vector<std::pair<Vertex, EdgeId>>>;

Adj adjacency_of(const Digraph& g) {
  Adj adj(g.num_vertices());
  for (Vertex v = 0; v < g.num_vertices(); ++v) {
    const auto heads = g.heads(v);
    const auto ids = g.arc_ids(v);
    for (std::size_t i = 0; i < heads.size(); ++i) adj[v].emplace_back(heads[i], ids[i]);
  }
  return adj;
}

template <class W>
Adj adjacency_of(const CsrGraph<W>& g) {
  Adj adj(g.num_vertices());
  for (Vertex v = 0; v < g.num_vertices(); ++v) {
    for (const auto& a : g.arcs(v)) adj[v].emplace_back(a.to, a.edge);
  }
  return adj;
}

// Recursive reference: preorder, children in adjacency order.
void ref_dfs(const Adj& adj, Vertex v, std::vector<char>& seen, std::vector<Vertex>& order,
             std::vector<Vertex>& parent, std::vector<EdgeId>& parent_edge) {
  seen[v] = 1;
  order.push_back(v);
  for (const auto& [x, id] : adj[v]) {
    if (!seen[x]) {
      parent[x] = v;
      parent_edge[x] = id;
      ref_dfs(adj, x, seen, order, parent, parent_edge);
    }
  }
}

// reach[u][v] = 1 iff v is reachable from u (Floyd-Warshall transitive closure).
std::vector<std::vector<char>> closure(const Adj& adj) {
  const std::size_t k = adj.size();
  std::vector<std::vector<char>> r(k, std::vector<char>(k, 0));
  for (std::size_t u = 0; u < k; ++u) {
    r[u][u] = 1;
    for (const auto& [x, id] : adj[u]) r[u][x] = 1;
  }
  for (std::size_t m = 0; m < k; ++m)
    for (std::size_t u = 0; u < k; ++u)
      if (r[u][m])
        for (std::size_t v = 0; v < k; ++v)
          if (r[m][v]) r[u][v] = 1;
  return r;
}

// Hop distances by repeated relaxation (independent of any queue order).
std::vector<std::size_t> hop_distances(const Adj& adj, Vertex s) {
  const std::size_t inf = adj.size() + 1;
  std::vector<std::size_t> h(adj.size(), inf);
  h[s] = 0;
  for (bool changed = true; changed;) {
    changed = false;
    for (Vertex u = 0; u < adj.size(); ++u) {
      if (h[u] == inf) continue;
      for (const auto& [x, id] : adj[u]) {
        if (h[u] + 1 < h[x]) {
          h[x] = h[u] + 1;
          changed = true;
        }
      }
    }
  }
  return h;
}

std::vector<DiArc> random_arcs(std::mt19937_64& rng, std::size_t k, std::size_t l) {
  std::uniform_int_distribution<Vertex> pick(0, static_cast<Vertex>(k - 1));
  std::vector<DiArc> arcs;
  for (std::size_t i = 0; i < l; ++i) arcs.push_back({pick(rng), pick(rng)});
  return arcs;
}

CsrGraph<std::uint64_t> random_undirected(std::mt19937_64& rng, std::size_t k, std::size_t l) {
  std::uniform_int_distribution<Vertex> pick(0, static_cast<Vertex>(k - 1));
  std::vector<sssp::Edge<std::uint64_t>> e;
  for (std::size_t i = 0; i < l; ++i) e.push_back({pick(rng), pick(rng), 1});
  return CsrGraph<std::uint64_t>::from_edges(k, std::move(e));
}

}  // namespace

TEST_CASE("Digraph keeps each tail's arcs in input order with their input indices",
          "[paper][search][digraph]") {
  const std::vector<DiArc> arcs{{2, 0}, {0, 1}, {2, 1}, {0, 2}, {2, 2}};
  const auto g = Digraph::from_arcs(4, arcs);
  REQUIRE(g.num_vertices() == 4);
  REQUIRE(g.num_arcs() == 5);
  REQUIRE(std::vector<Vertex>(g.heads(0).begin(), g.heads(0).end()) == std::vector<Vertex>{1, 2});
  REQUIRE(std::vector<EdgeId>(g.arc_ids(0).begin(), g.arc_ids(0).end()) ==
          std::vector<EdgeId>{1, 3});
  REQUIRE(std::vector<Vertex>(g.heads(2).begin(), g.heads(2).end()) ==
          std::vector<Vertex>{0, 1, 2});
  REQUIRE(std::vector<EdgeId>(g.arc_ids(2).begin(), g.arc_ids(2).end()) ==
          std::vector<EdgeId>{0, 2, 4});
  REQUIRE(g.heads(1).empty());
  REQUIRE(g.heads(3).empty());
  REQUIRE_THROWS_AS(Digraph::from_arcs(2, std::vector<DiArc>{{0, 2}}), std::out_of_range);
  REQUIRE(Digraph::from_arcs(0, std::vector<DiArc>{}).num_vertices() == 0);
}

TEST_CASE("DFS and BFS discovery orders on a hand example", "[paper][search][dfs][bfs]") {
  // 0-1, 0-2, 1-3, 2-3, 3-4, and isolated 5.
  const auto g = CsrGraph<std::uint64_t>::from_edges(
      6, {{0, 1, 1}, {0, 2, 1}, {1, 3, 1}, {2, 3, 1}, {3, 4, 1}});
  const auto d = sssp::paper::dfs_order(g, 0);
  REQUIRE(d.order == std::vector<Vertex>{0, 1, 3, 2, 4});
  REQUIRE(d.parent == std::vector<Vertex>{kNoVertex, 0, 3, 1, 3, kNoVertex});
  REQUIRE(d.parent_edge == std::vector<EdgeId>{kNoEdge, 0, 3, 2, 4, kNoEdge});
  const auto b = sssp::paper::bfs_order(g, 0);
  REQUIRE(b.order == std::vector<Vertex>{0, 1, 2, 3, 4});
  REQUIRE(b.parent == std::vector<Vertex>{kNoVertex, 0, 0, 1, 3, kNoVertex});
  REQUIRE(b.parent_edge == std::vector<EdgeId>{kNoEdge, 0, 1, 2, 4, kNoEdge});
  const auto iso = sssp::paper::dfs_order(g, 5);
  REQUIRE(iso.order == std::vector<Vertex>{5});
}

TEST_CASE("DFS equals the recursive preorder on random graphs", "[paper][search][dfs]") {
  for (std::uint64_t seed = 0; seed < 40; ++seed) {
    std::mt19937_64 rng(seed);
    const std::size_t k = 1 + seed % 23;
    const std::size_t l = (seed * 7) % 60;
    const auto dg = Digraph::from_arcs(k, random_arcs(rng, k, l));
    const auto ug = random_undirected(rng, k, l);
    const auto check = [&](const auto& g, const Adj& adj) {
      for (Vertex root = 0; root < k; ++root) {
        std::vector<char> seen(k, 0);
        std::vector<Vertex> order;
        std::vector<Vertex> parent(k, kNoVertex);
        std::vector<EdgeId> parent_edge(k, kNoEdge);
        ref_dfs(adj, root, seen, order, parent, parent_edge);
        const auto t = sssp::paper::dfs_order(g, root);
        REQUIRE(t.order == order);
        REQUIRE(t.parent == parent);
        REQUIRE(t.parent_edge == parent_edge);
      }
    };
    check(dg, adjacency_of(dg));
    check(ug, adjacency_of(ug));
  }
}

TEST_CASE("BFS visits the reachable set in nondecreasing hop distance", "[paper][search][bfs]") {
  for (std::uint64_t seed = 0; seed < 40; ++seed) {
    std::mt19937_64 rng(1000 + seed);
    const std::size_t k = 1 + seed % 19;
    const std::size_t l = (seed * 5) % 50;
    const auto dg = Digraph::from_arcs(k, random_arcs(rng, k, l));
    const auto ug = random_undirected(rng, k, l);
    const auto check = [&](const auto& g, const Adj& adj) {
      const auto reach = closure(adj);
      for (Vertex root = 0; root < k; ++root) {
        const auto h = hop_distances(adj, root);
        const auto t = sssp::paper::bfs_order(g, root);
        REQUIRE(t.order.front() == root);
        std::vector<char> in_order(k, 0);
        for (std::size_t i = 0; i < t.order.size(); ++i) {
          const Vertex v = t.order[i];
          REQUIRE(in_order[v] == 0);  // each vertex discovered once
          in_order[v] = 1;
          if (i > 0) REQUIRE(h[t.order[i - 1]] <= h[v]);
        }
        for (Vertex v = 0; v < k; ++v) {
          REQUIRE(static_cast<bool>(in_order[v]) == static_cast<bool>(reach[root][v]));
          if (v == root || !in_order[v]) {
            REQUIRE(t.parent[v] == kNoVertex);
            continue;
          }
          // The parent arc exists and is one hop closer to the root.
          const Vertex p = t.parent[v];
          REQUIRE(h[p] + 1 == h[v]);
          const auto& out = adj[p];
          REQUIRE(std::find(out.begin(), out.end(), std::make_pair(v, t.parent_edge[v])) !=
                  out.end());
        }
      }
    };
    check(dg, adjacency_of(dg));
    check(ug, adjacency_of(ug));
  }
}

#if SSSP_DIAGNOSTICS  // counters are compiled out when SSSP_DIAGNOSTICS=0
TEST_CASE("searches examine each adjacency entry of the reached part once",
          "[paper][search][lemma3.1]") {
  std::mt19937_64 rng(7);
  // A connected undirected graph: a spanning path plus random extra edges.
  const std::size_t k = 300;
  std::vector<sssp::Edge<std::uint64_t>> e;
  for (Vertex v = 1; v < k; ++v) e.push_back({v - 1, v, 1});
  std::uniform_int_distribution<Vertex> pick(0, k - 1);
  for (int i = 0; i < 900; ++i) e.push_back({pick(rng), pick(rng), 1});
  const auto g = CsrGraph<std::uint64_t>::from_edges(k, std::move(e));
  for (int which = 0; which < 2; ++which) {
    sssp::OpCounters ops;
    const auto t =
        which == 0 ? sssp::paper::dfs_order(g, 0, &ops) : sssp::paper::bfs_order(g, 0, &ops);
    REQUIRE(t.order.size() == k);
    REQUIRE(ops.vertex_scans == k);
    REQUIRE(ops.edge_scans == g.num_arcs());
  }
}
#endif

TEST_CASE("DFS uses an explicit stack (long path does not overflow the call stack)",
          "[paper][search][dfs]") {
  const std::size_t k = 200000;
  std::vector<DiArc> arcs;
  for (Vertex v = 1; v < k; ++v) arcs.push_back({v - 1, v});
  const auto g = Digraph::from_arcs(k, arcs);
  const auto t = sssp::paper::dfs_order(g, 0);
  std::vector<Vertex> expect(k);
  std::iota(expect.begin(), expect.end(), Vertex{0});
  REQUIRE(t.order == expect);
}

TEST_CASE("searches reject an out-of-range root", "[paper][search]") {
  const auto g = CsrGraph<double>::from_edges(0, {});
  REQUIRE_THROWS_AS(sssp::paper::dfs_order(g, 0), std::out_of_range);
  REQUIRE_THROWS_AS(sssp::paper::bfs_order(g, 0), std::out_of_range);
  const auto d = Digraph::from_arcs(2, std::vector<DiArc>{{0, 1}});
  REQUIRE_THROWS_AS(sssp::paper::dfs_order(d, 2), std::out_of_range);
  REQUIRE_THROWS_AS(sssp::paper::bfs_order(d, 2), std::out_of_range);
}

TEST_CASE("connected components label exactly the mutually reachable classes",
          "[paper][search][components]") {
  SECTION("empty graph") {
    const auto c = sssp::paper::connected_components(CsrGraph<double>::from_edges(0, {}));
    REQUIRE(c.count == 0);
    REQUIRE(c.label.empty());
  }
  SECTION("single vertex with a self-loop") {
    const auto c =
        sssp::paper::connected_components(CsrGraph<std::uint64_t>::from_edges(1, {{0, 0, 3}}));
    REQUIRE(c.count == 1);
    REQUIRE(c.label == std::vector<Vertex>{0});
  }
  SECTION("hand example, labels numbered by smallest member") {
    const auto c = sssp::paper::connected_components(
        CsrGraph<std::uint64_t>::from_edges(7, {{5, 3, 1}, {1, 6, 0}, {3, 0, 2}, {6, 6, 1}}));
    REQUIRE(c.count == 4);
    REQUIRE(c.label == std::vector<Vertex>{0, 1, 2, 0, 3, 0, 1});
  }
  SECTION("random graphs versus transitive closure") {
    for (std::uint64_t seed = 0; seed < 40; ++seed) {
      std::mt19937_64 rng(2000 + seed);
      const std::size_t k = 1 + seed % 25;
      const auto g = random_undirected(rng, k, seed % 30);
      const auto reach = closure(adjacency_of(g));
      const auto c = sssp::paper::connected_components(g);
      REQUIRE(c.label.size() == k);
      Vertex next = 0;
      for (Vertex u = 0; u < k; ++u) {
        REQUIRE(c.label[u] <= next);  // first appearance order
        if (c.label[u] == next) ++next;
        for (Vertex v = 0; v < k; ++v) REQUIRE((c.label[u] == c.label[v]) == (reach[u][v] == 1));
      }
      REQUIRE(c.count == next);
    }
  }
}

TEST_CASE("topological sort orders every arc forward", "[paper][search][topo]") {
  SECTION("empty and single vertex") {
    REQUIRE(sssp::paper::topological_sort(Digraph::from_arcs(0, std::vector<DiArc>{})).empty());
    REQUIRE(sssp::paper::topological_sort(Digraph::from_arcs(1, std::vector<DiArc>{})) ==
            std::vector<Vertex>{0});
  }
  SECTION("hand example with parallel arcs") {
    const auto g = Digraph::from_arcs(4, std::vector<DiArc>{{3, 1}, {3, 1}, {1, 0}, {2, 0}});
    // Kahn with a FIFO of zero-indegree vertices, seeded in index order.
    REQUIRE(sssp::paper::topological_sort(g) == std::vector<Vertex>{2, 3, 1, 0});
  }
  SECTION("random DAGs") {
    for (std::uint64_t seed = 0; seed < 50; ++seed) {
      std::mt19937_64 rng(3000 + seed);
      const std::size_t k = 1 + seed % 40;
      std::vector<Vertex> rank(k);
      std::iota(rank.begin(), rank.end(), Vertex{0});
      std::shuffle(rank.begin(), rank.end(), rng);
      std::vector<DiArc> arcs;
      for (const auto& a : random_arcs(rng, k, 3 * k)) {
        if (rank[a.tail] < rank[a.head]) arcs.push_back(a);
      }
      const auto order = sssp::paper::topological_sort(Digraph::from_arcs(k, arcs));
      REQUIRE(order.size() == k);
      std::vector<std::size_t> pos(k, k);
      for (std::size_t i = 0; i < k; ++i) {
        REQUIRE(pos[order[i]] == k);  // a permutation
        pos[order[i]] = i;
      }
      for (const auto& a : arcs) REQUIRE(pos[a.tail] < pos[a.head]);
    }
  }
  SECTION("cycles are rejected") {
    REQUIRE_THROWS_AS(
        sssp::paper::topological_sort(Digraph::from_arcs(1, std::vector<DiArc>{{0, 0}})),
        std::invalid_argument);
    REQUIRE_THROWS_AS(sssp::paper::topological_sort(Digraph::from_arcs(
                          4, std::vector<DiArc>{{0, 1}, {1, 2}, {2, 3}, {3, 1}})),
                      std::invalid_argument);
  }
#if SSSP_DIAGNOSTICS
  SECTION("each arc is examined once") {
    std::vector<DiArc> arcs;
    for (Vertex v = 0; v + 1 < 50; ++v) {
      arcs.push_back({v, v + 1});
      arcs.push_back({v, v + 1});
    }
    sssp::OpCounters ops;
    const auto order = sssp::paper::topological_sort(Digraph::from_arcs(50, arcs), &ops);
    REQUIRE(order.size() == 50);
    REQUIRE(ops.edge_scans == arcs.size());
    REQUIRE(ops.vertex_scans == 50);
  }
#endif
}
