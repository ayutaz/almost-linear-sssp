// Tests for the base case of Algorithm 4 (l.2-5; Sec 2 l.286-290): one native state per
// vertex, both orientations of every nonloop edge, and the equal-weight-simultaneous Kruskal
// merge tree of the MST (after keeping a lightest edge per vertex pair) with MST edges as links.
// Distances are checked with naive Floyd-Warshall written here.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/base_case.hpp"
#include "paper/primitives/edge_pairs.hpp"
#include "paper/primitives/merge_events.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
using sssp::paper::kNoNode;
using sssp::paper::NodeId;
using sssp::paper::OrderedHierarchy;
using sssp::paper::StateId;

namespace {

template <class W>
constexpr W kInf = sssp::WeightTraits<W>::infinity();

// All-pairs shortest distances (sum) on the undirected graph.
template <class W>
std::vector<std::vector<W>> fw_sum(std::size_t n, const std::vector<Edge<W>>& edges) {
  std::vector<std::vector<W>> d(n, std::vector<W>(n, kInf<W>));
  for (std::size_t i = 0; i < n; ++i) d[i][i] = W{0};
  for (const auto& e : edges) {
    if (e.weight < d[e.u][e.v]) d[e.u][e.v] = d[e.v][e.u] = e.weight;
  }
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        if (d[i][k] != kInf<W> && d[k][j] != kInf<W> && d[i][k] + d[k][j] < d[i][j])
          d[i][j] = d[i][k] + d[k][j];
  return d;
}

// All-pairs shortest distances on a directed state graph.
template <class W>
std::vector<std::vector<W>> fw_states(const sssp::paper::StateGraph<W>& g) {
  const std::size_t n = g.num_states;
  std::vector<std::vector<W>> d(n, std::vector<W>(n, kInf<W>));
  for (std::size_t i = 0; i < n; ++i) d[i][i] = W{0};
  for (const auto& a : g.arcs) d[a.tail][a.head] = std::min(d[a.tail][a.head], a.length);
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        if (d[i][k] != kInf<W> && d[k][j] != kInf<W> && d[i][k] + d[k][j] < d[i][j])
          d[i][j] = d[i][k] + d[k][j];
  return d;
}

// All-pairs minimax (bottleneck) path weights on the undirected graph, ignoring self-loops.
template <class W>
std::vector<std::vector<W>> fw_bottleneck(std::size_t n, const std::vector<Edge<W>>& edges) {
  std::vector<std::vector<W>> b(n, std::vector<W>(n, kInf<W>));
  for (std::size_t i = 0; i < n; ++i) b[i][i] = W{0};
  for (const auto& e : edges) {
    if (e.u != e.v && e.weight < b[e.u][e.v]) b[e.u][e.v] = b[e.v][e.u] = e.weight;
  }
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j) {
        const W via = std::max(b[i][k], b[k][j]);
        if (via < b[i][j]) b[i][j] = via;
      }
  return b;
}

template <class W>
NodeId lca(const OrderedHierarchy<W>& h, NodeId a, NodeId b) {
  std::vector<NodeId> up;
  for (NodeId v = a; v != kNoNode; v = h.nodes[v].parent) up.push_back(v);
  for (NodeId v = b; v != kNoNode; v = h.nodes[v].parent) {
    if (std::find(up.begin(), up.end(), v) != up.end()) return v;
  }
  return kNoNode;
}

template <class W>
std::size_t internal_count(const OrderedHierarchy<W>& h) {
  return static_cast<std::size_t>(std::count_if(h.nodes.begin(), h.nodes.end(),
                                                [](const auto& v) { return !v.children.empty(); }));
}

}  // namespace

TEMPLATE_TEST_CASE("base case merges equal-weight MST edges simultaneously",
                   "[paper][local][base_case]", std::uint64_t, double) {
  using W = TestType;
  SECTION("path with two equal edges gives one ternary node") {
    const auto g = CsrGraph<W>::from_edges(3, {{0, 1, W{3}}, {1, 2, W{3}}});
    const auto r = sssp::paper::local_hierarchy_base_case(g);
    const auto& h = r.hierarchy;
    REQUIRE(internal_count(h) == 1);
    REQUIRE(h.roots.size() == 1);
    const auto& root = h.nodes[h.roots[0]];
    REQUIRE(root.height == W{3});
    REQUIRE(root.children.size() == 3);
    REQUIRE(root.links.size() == 2);
    for (const auto& lk : root.links) REQUIRE(lk.length == W{3});
  }
  SECTION("equal-weight triangle: one node, two of the three edges are links") {
    const auto g = CsrGraph<W>::from_edges(3, {{0, 1, W{2}}, {1, 2, W{2}}, {0, 2, W{2}}});
    const auto r = sssp::paper::local_hierarchy_base_case(g);
    REQUIRE(internal_count(r.hierarchy) == 1);
    const auto& root = r.hierarchy.nodes[r.hierarchy.roots[0]];
    REQUIRE(root.children.size() == 3);
    REQUIRE(root.links.size() == 2);
    REQUIRE(r.hierarchy.graph.arcs.size() == 6);  // all three edges stay as arcs
  }
  SECTION("two light pairs joined by a heavy edge give two nodes under a root") {
    const auto g = CsrGraph<W>::from_edges(4, {{0, 1, W{1}}, {2, 3, W{1}}, {1, 2, W{5}}});
    const auto r = sssp::paper::local_hierarchy_base_case(g);
    const auto& h = r.hierarchy;
    REQUIRE(internal_count(h) == 3);
    REQUIRE(h.roots.size() == 1);
    const auto& root = h.nodes[h.roots[0]];
    REQUIRE(root.height == W{5});
    REQUIRE(root.children.size() == 2);
    REQUIRE(root.links.size() == 1);
    for (const NodeId c : root.children) {
      REQUIRE(h.nodes[c].height == W{1});
      REQUIRE(h.nodes[c].children.size() == 2);
    }
  }
}

TEMPLATE_TEST_CASE("base case output satisfies Definition 4.1 with c0 = 1 and no light arcs",
                   "[paper][local][base_case][random]", std::uint64_t, double) {
  using W = TestType;
  std::size_t merged_ties = 0;  // internal nodes with >= 3 children (simultaneous merges)
  for (std::uint64_t seed = 1; seed <= 300; ++seed) {
    std::mt19937_64 rng(seed);
    const std::size_t n = 1 + rng() % 12;
    const std::size_t m = rng() % (3 * n + 1);
    const std::uint64_t max_w = 1 + rng() % 5;  // many ties
    std::vector<Edge<W>> edges;
    for (std::size_t i = 0; i < m; ++i) {
      edges.push_back({static_cast<Vertex>(rng() % n), static_cast<Vertex>(rng() % n),
                       static_cast<W>(1 + rng() % max_w)});
    }
    const auto g = CsrGraph<W>::from_edges(n, edges);
    const auto r = sssp::paper::local_hierarchy_base_case(g);
    const auto& h = r.hierarchy;
    INFO("seed " << seed << " n " << n << " m " << m);

    const auto rep = sssp::paper::validate_ordered_hierarchy(h, g);
    INFO(rep.summary());
    REQUIRE(rep.ok());
    REQUIRE(h.c0 == 1);

    // One native entry state per vertex, mapped to that vertex.
    REQUIRE(h.num_states() == n);
    REQUIRE(r.native_entry.size() == n);
    for (Vertex v = 0; v < n; ++v) REQUIRE(h.projection[r.native_entry[v]] == v);

    // Both orientations of every nonloop edge (parallel edges kept), no self-loops.
    std::size_t nonloop = 0;
    for (const auto& e : edges) nonloop += e.u != e.v ? 1 : 0;
    REQUIRE(h.graph.arcs.size() == 2 * nonloop);
    REQUIRE(r.arc_edge.size() == h.graph.arcs.size());
    for (std::size_t i = 0; i < h.graph.arcs.size(); ++i) {
      const auto& a = h.graph.arcs[i];
      const auto& e = edges[r.arc_edge[i]];
      REQUIRE(a.tail != a.head);
      REQUIRE(a.length == e.weight);
      const bool same = (h.projection[a.tail] == e.u && h.projection[a.head] == e.v) ||
                        (h.projection[a.tail] == e.v && h.projection[a.head] == e.u);
      REQUIRE(same);
    }

    // No light arcs (Lemma 3.5): every arc between distinct children of its LCA has length
    // >= the LCA's height.
    for (const auto& a : h.graph.arcs) {
      const NodeId l = lca(h, h.leaf_of_state[a.tail], h.leaf_of_state[a.head]);
      REQUIRE(l != kNoNode);
      REQUIRE(a.length >= h.nodes[l].height);
    }

    // Exact heights: h(LCA(x, y)) equals the minimax path weight between x and y; vertices in
    // different components lie in different trees. Heights strictly increase upwards.
    const auto b = fw_bottleneck<W>(n, edges);
    for (Vertex x = 0; x < n; ++x) {
      for (Vertex y = x + 1; y < n; ++y) {
        const NodeId l =
            lca(h, h.leaf_of_state[r.native_entry[x]], h.leaf_of_state[r.native_entry[y]]);
        if (b[x][y] == kInf<W>) {
          REQUIRE(l == kNoNode);
        } else {
          REQUIRE(l != kNoNode);
          REQUIRE(h.nodes[l].height == b[x][y]);
        }
      }
    }
    for (const auto& node : h.nodes) {
      if (!node.children.empty() && node.parent != kNoNode) {
        REQUIRE(node.height < h.nodes[node.parent].height);
      }
      if (node.children.size() >= 3) ++merged_ties;
      for (const auto& lk : node.links) REQUIRE(lk.length == node.height);  // MST edges, c0 = 1
    }

    // Every geodesic is realised by states from the native entries (Lemma 5.10).
    const auto dp = fw_sum<W>(n, edges);
    const auto ds = fw_states(h.graph);
    for (Vertex x = 0; x < n; ++x) {
      for (Vertex y = 0; y < n; ++y) REQUIRE(ds[r.native_entry[x]][r.native_entry[y]] == dp[x][y]);
    }
  }
  REQUIRE(merged_ties >= 30);  // the random family exercises simultaneous merges
}

TEST_CASE("base case on a single vertex is one leaf root", "[paper][local][base_case]") {
  const auto g = CsrGraph<std::uint64_t>::from_edges(1, {{0, 0, 4}});  // self-loop is dropped
  const auto r = sssp::paper::local_hierarchy_base_case(g);
  const auto& h = r.hierarchy;
  REQUIRE(h.num_states() == 1);
  REQUIRE(h.nodes.size() == 1);
  REQUIRE(h.roots == std::vector<NodeId>{0});
  REQUIRE(h.graph.arcs.empty());
  REQUIRE(r.native_entry == std::vector<StateId>{0});
  REQUIRE(sssp::paper::validate_ordered_hierarchy(h, g).ok());
}

TEST_CASE("base case on a disconnected graph gives one tree per component",
          "[paper][local][base_case]") {
  // {0,1} joined, {2,3} joined, 4 isolated.
  const auto g = CsrGraph<std::uint64_t>::from_edges(5, {{0, 1, 2}, {2, 3, 7}});
  const auto r = sssp::paper::local_hierarchy_base_case(g);
  const auto& h = r.hierarchy;
  REQUIRE(h.roots.size() == 3);
  REQUIRE(internal_count(h) == 2);
  REQUIRE(h.nodes[h.leaf_of_state[4]].parent == kNoNode);
  REQUIRE(sssp::paper::validate_ordered_hierarchy(h, g).ok());
}

TEST_CASE("base case keeps parallel arcs but links through a lightest edge per pair",
          "[paper][local][base_case]") {
  // Edge ids: 0 = {0,1} w5, 1 = {0,1} w2, 2 = loop at 1, 3 = {1,2} w2.
  const auto g =
      CsrGraph<std::uint64_t>::from_edges(3, {{0, 1, 5}, {1, 0, 2}, {1, 1, 1}, {1, 2, 2}});
  const auto r = sssp::paper::local_hierarchy_base_case(g);
  const auto& h = r.hierarchy;
  REQUIRE(h.graph.arcs.size() == 6);
  for (const auto& a : h.graph.arcs) REQUIRE(a.tail != a.head);
  REQUIRE(internal_count(h) == 1);
  const auto& root = h.nodes[h.roots[0]];
  REQUIRE(root.height == 2);
  REQUIRE(root.children.size() == 3);
  REQUIRE(root.links.size() == 2);
  for (const auto& lk : root.links) {
    REQUIRE(lk.length == 2);
    REQUIRE(lk.arc != sssp::paper::kNoArc);
    REQUIRE(r.arc_edge[lk.arc] != 0);  // never the heavy parallel edge
  }
  REQUIRE(sssp::paper::validate_ordered_hierarchy(h, g).ok());
}

TEST_CASE("base case rejects zero weights (zero components are contracted by the caller)",
          "[paper][local][base_case]") {
  const auto g = CsrGraph<std::uint64_t>::from_edges(2, {{0, 1, 0}});
  REQUIRE_THROWS_AS(sssp::paper::local_hierarchy_base_case(g), std::invalid_argument);
  const auto gl = CsrGraph<double>::from_edges(2, {{0, 1, 1.0}, {1, 1, 0.0}});
  REQUIRE_THROWS_AS(sssp::paper::local_hierarchy_base_case(gl), std::invalid_argument);
}

TEST_CASE("base case records path counters in Diagnostics", "[paper][local][base_case]") {
  const auto g = CsrGraph<std::uint64_t>::from_edges(4, {{0, 1, 1}, {2, 3, 1}, {1, 2, 5}});
  sssp::Diagnostics diag;
  const auto r = sssp::paper::local_hierarchy_base_case(g, &diag);
  REQUIRE(r.hierarchy.nodes.size() == 7);
#if SSSP_DIAGNOSTICS
  REQUIRE(diag.counters["local_base_case_calls"] == 1);
  REQUIRE(diag.counters["local_base_case_states"] == 4);
  REQUIRE(diag.counters["local_base_case_merge_nodes"] == 3);
#else
  REQUIRE(diag.counters.empty());  // counters are compiled out in timed builds
#endif
}

TEMPLATE_TEST_CASE("base case hierarchy is the merge forest of the lightest edges per pair",
                   "[paper][local][base_case][random]", std::uint64_t, double) {
  // Alg 4 l.4 = lightest edge per pair (P3-EDGE-PAIRS) + equal-weight simultaneous Kruskal
  // merges (P3-MERGE-EVENTS): node i of the hierarchy is node i of the merge forest, with the
  // same height, the same children in increasing id order, and its forest edges as links.
  using W = TestType;
  for (std::uint64_t seed = 1; seed <= 300; ++seed) {
    std::mt19937_64 rng(seed);
    const std::size_t n = 1 + rng() % 16;
    const std::size_t m = rng() % (4 * n + 1);
    const std::uint64_t max_w = 1 + rng() % 4;
    std::vector<Edge<W>> edges;
    for (std::size_t i = 0; i < m; ++i) {
      edges.push_back({static_cast<Vertex>(rng() % n), static_cast<Vertex>(rng() % n),
                       static_cast<W>(1 + rng() % max_w)});
    }
    const auto g = CsrGraph<W>::from_edges(n, edges);
    const auto r = sssp::paper::local_hierarchy_base_case(g);
    const auto& h = r.hierarchy;
    const auto kept = sssp::paper::lightest_edge_per_pair<W>(n, g.edges());
    const auto f = sssp::paper::build_merge_events<W>(n, g.edges(), kept);
    INFO("seed " << seed << " n " << n << " m " << m);
    REQUIRE(h.nodes.size() == f.num_nodes());
    REQUIRE(h.roots == f.roots);
    for (std::size_t i = 0; i < f.internal.size(); ++i) {
      const auto& node = h.nodes[n + i];
      const auto& want = f.internal[i];
      REQUIRE(node.height == want.height);
      REQUIRE(node.children == want.children);
      REQUIRE(node.links.size() == want.edges.size());
      for (std::size_t k = 0; k < want.edges.size(); ++k) {
        const auto& e = g.edges()[want.edges[k]];
        REQUIRE(node.links[k].x == e.u);
        REQUIRE(node.links[k].y == e.v);
        REQUIRE(node.links[k].length == e.weight);
        REQUIRE(r.arc_edge[node.links[k].arc] == want.edges[k]);
      }
    }
  }
}
