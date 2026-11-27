// Tests for the Kruskal merge forest with equal weights merged simultaneously (Sec 2 l.288-290,
// Sec 3.3 ordered merge events, Lemma 3.5; used by Alg 4 l.4 and the Y trees of Sec 5.4).
// The references below are written here independently of the library: components by naive
// relabelling for every weight threshold, and minimax distances by a Floyd-Warshall variant.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/primitives/merge_events.hpp"
#include "paper/primitives/mst.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::Edge;
using sssp::EdgeId;
using sssp::Vertex;
using sssp::paper::build_merge_events;
using sssp::paper::kNoMergeNode;
using sssp::paper::MergeForest;
using sssp::paper::MergeNodeId;

namespace {

template <class W>
std::vector<Edge<W>> random_multigraph(std::size_t n, std::size_t m, std::uint64_t max_weight,
                                       std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::vector<Edge<W>> edges;
  for (std::size_t i = 0; i < m; ++i) {
    const auto u = static_cast<Vertex>(rng() % n);
    const auto v = static_cast<Vertex>(rng() % n);
    edges.push_back({u, v, static_cast<W>(1 + rng() % max_weight)});
  }
  return edges;
}

// Vertices below every node (leaves are the vertices themselves), each sorted.
template <class W>
std::vector<std::vector<Vertex>> vertex_sets(const MergeForest<W>& f) {
  std::vector<std::vector<Vertex>> s(f.num_nodes());
  for (std::size_t v = 0; v < f.num_leaves; ++v) s[v] = {static_cast<Vertex>(v)};
  for (std::size_t i = 0; i < f.internal.size(); ++i) {
    auto& own = s[f.num_leaves + i];
    for (const MergeNodeId c : f.internal[i].children) {
      REQUIRE(c < f.num_leaves + i);  // children are created before their parent
      own.insert(own.end(), s[c].begin(), s[c].end());
    }
    std::sort(own.begin(), own.end());
  }
  return s;
}

// Components of the graph on `n` vertices using the edges with pred(weight) true; label per
// vertex = least vertex of its component. Naive relabelling, O(n m).
template <class W, class Pred>
std::vector<Vertex> components(std::size_t n, const std::vector<Edge<W>>& edges, Pred pred) {
  std::vector<Vertex> label(n);
  std::iota(label.begin(), label.end(), Vertex{0});
  bool changed = true;
  while (changed) {
    changed = false;
    for (const auto& e : edges) {
      if (!pred(e.weight)) continue;
      const Vertex m = std::min(label[e.u], label[e.v]);
      if (label[e.u] != m || label[e.v] != m) {
        label[e.u] = label[e.v] = m;
        changed = true;
      }
    }
  }
  return label;
}

template <class W>
std::vector<std::vector<Vertex>> groups_of(const std::vector<Vertex>& label) {
  std::vector<std::vector<Vertex>> g(label.size());
  for (std::size_t v = 0; v < label.size(); ++v) g[label[v]].push_back(static_cast<Vertex>(v));
  return g;
}

// Reference node set: (w, C) for every component C of G_{<=w} that contains at least two
// components of G_{<w}, over all distinct weights w. Returned as (w, C) pairs, sorted.
template <class W>
std::vector<std::pair<W, std::vector<Vertex>>> reference_nodes(std::size_t n,
                                                               const std::vector<Edge<W>>& edges) {
  std::set<W> weights;
  for (const auto& e : edges) {
    if (e.u != e.v) weights.insert(e.weight);
  }
  std::vector<std::pair<W, std::vector<Vertex>>> out;
  for (const W w : weights) {
    const auto below = components(n, edges, [&](W x) { return x < w; });
    const auto upto = components(n, edges, [&](W x) { return !(w < x); });
    for (const auto& comp : groups_of<W>(upto)) {
      if (comp.empty()) continue;
      std::set<Vertex> parts;
      for (const Vertex v : comp) parts.insert(below[v]);
      if (parts.size() >= 2) out.emplace_back(w, comp);
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

// Minimax (bottleneck) path value between all pairs; +infinity if disconnected.
template <class W>
std::vector<std::vector<W>> minimax(std::size_t n, const std::vector<Edge<W>>& edges) {
  const W inf = std::numeric_limits<W>::max();
  std::vector<std::vector<W>> d(n, std::vector<W>(n, inf));
  for (const auto& e : edges) {
    if (e.u != e.v && e.weight < d[e.u][e.v]) d[e.u][e.v] = d[e.v][e.u] = e.weight;
  }
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        const W via = std::max(d[i][k], d[k][j]);
        if (via < d[i][j]) d[i][j] = via;
      }
    }
  }
  return d;
}

template <class W>
MergeNodeId lca(const MergeForest<W>& f, MergeNodeId a, MergeNodeId b) {
  std::vector<bool> seen(f.num_nodes(), false);
  for (MergeNodeId x = a; x != kNoMergeNode; x = f.parent[x]) seen[x] = true;
  for (MergeNodeId x = b; x != kNoMergeNode; x = f.parent[x]) {
    if (seen[x]) return x;
  }
  return kNoMergeNode;
}

// Structural checks that hold for every output.
template <class W>
void check_structure(const MergeForest<W>& f, std::size_t n, const std::vector<Edge<W>>& edges) {
  REQUIRE(f.num_leaves == n);
  REQUIRE(f.parent.size() == f.num_nodes());
  std::size_t total_edges = 0;
  for (std::size_t i = 0; i < f.internal.size(); ++i) {
    const auto id = static_cast<MergeNodeId>(n + i);
    const auto& node = f.internal[i];
    REQUIRE(node.children.size() >= 2);
    REQUIRE(std::is_sorted(node.children.begin(), node.children.end()));
    REQUIRE(node.edges.size() + 1 == node.children.size());
    for (const EdgeId e : node.edges) REQUIRE(edges[e].weight == node.height);
    for (const MergeNodeId c : node.children) {
      REQUIRE(f.parent[c] == id);
      // Equal weights merge simultaneously: an internal child is strictly lower.
      if (c >= n) REQUIRE(f.internal[c - n].height < node.height);
    }
    if (i > 0) REQUIRE_FALSE(node.height < f.internal[i - 1].height);
    total_edges += node.edges.size();
  }
  REQUIRE(total_edges == f.forest_edges.size());
  std::vector<MergeNodeId> roots;
  for (std::size_t x = 0; x < f.num_nodes(); ++x) {
    if (f.parent[x] == kNoMergeNode) roots.push_back(static_cast<MergeNodeId>(x));
  }
  REQUIRE(f.roots == roots);
}

}  // namespace

TEMPLATE_TEST_CASE("merge events on hand examples", "[paper][merge_events]", std::uint64_t,
                   double) {
  using T = TestType;
  SECTION("path with two equal edges gives one ternary node") {
    const std::vector<Edge<T>> e{{0, 1, T(2)}, {1, 2, T(2)}};
    const auto f = build_merge_events<T>(3, e);
    check_structure(f, 3, e);
    REQUIRE(f.internal.size() == 1);
    REQUIRE(f.internal[0].height == T(2));
    REQUIRE(f.internal[0].children == std::vector<MergeNodeId>{0, 1, 2});
    REQUIRE(f.internal[0].edges == std::vector<EdgeId>{0, 1});
    REQUIRE(f.roots == std::vector<MergeNodeId>{3});
  }
  SECTION("equal-weight triangle: one node, the two least-index edges join it") {
    const std::vector<Edge<T>> e{{2, 0, T(1)}, {0, 1, T(1)}, {1, 2, T(1)}};
    const auto f = build_merge_events<T>(3, e);
    check_structure(f, 3, e);
    REQUIRE(f.internal.size() == 1);
    REQUIRE(f.internal[0].children == std::vector<MergeNodeId>{0, 1, 2});
    REQUIRE(f.internal[0].edges == std::vector<EdgeId>{0, 1});
    REQUIRE(f.forest_edges == std::vector<EdgeId>{0, 1});
  }
  SECTION("two light pairs joined by a heavy edge") {
    const std::vector<Edge<T>> e{{2, 3, T(1)}, {0, 1, T(1)}, {1, 2, T(5)}};
    const auto f = build_merge_events<T>(4, e);
    check_structure(f, 4, e);
    REQUIRE(f.internal.size() == 3);
    // Nodes of one weight class appear in the order of their first forest edge (id order).
    REQUIRE(f.internal[0].children == std::vector<MergeNodeId>{2, 3});
    REQUIRE(f.internal[1].children == std::vector<MergeNodeId>{0, 1});
    REQUIRE(f.internal[2].height == T(5));
    REQUIRE(f.internal[2].children == std::vector<MergeNodeId>{4, 5});
    REQUIRE(f.internal[2].edges == std::vector<EdgeId>{2});
    REQUIRE(f.forest_edges == std::vector<EdgeId>{0, 1, 2});
  }
  SECTION("a node may merge an internal child with leaves") {
    // 0-1 at 1, then 1-2 and 3-0 at 4: one node of height 4 over {node(0,1), 2, 3}.
    const std::vector<Edge<T>> e{{1, 2, T(4)}, {0, 1, T(1)}, {3, 0, T(4)}, {2, 3, T(7)}};
    const auto f = build_merge_events<T>(4, e);
    check_structure(f, 4, e);
    REQUIRE(f.internal.size() == 2);
    REQUIRE(f.internal[1].height == T(4));
    REQUIRE(f.internal[1].children == std::vector<MergeNodeId>{2, 3, 4});
    REQUIRE(f.internal[1].edges == std::vector<EdgeId>{0, 2});
  }
}

TEMPLATE_TEST_CASE("merge events boundary inputs", "[paper][merge_events]", std::uint64_t, double) {
  using T = TestType;
  const auto empty = build_merge_events<T>(0, std::vector<Edge<T>>{});
  REQUIRE(empty.num_nodes() == 0);
  REQUIRE(empty.roots.empty());
  // A single vertex with a loop: one leaf root, no internal node.
  const std::vector<Edge<T>> loop{{0, 0, T(3)}};
  const auto one = build_merge_events<T>(1, loop);
  REQUIRE(one.internal.empty());
  REQUIRE(one.roots == std::vector<MergeNodeId>{0});
  REQUIRE(one.forest_edges.empty());
  // Disconnected: one root per component (isolated vertex 2 is a leaf root).
  const std::vector<Edge<T>> two{{0, 1, T(1)}, {3, 4, T(2)}, {4, 3, T(1)}};
  const auto f = build_merge_events<T>(5, two);
  check_structure(f, 5, two);
  REQUIRE(f.roots == std::vector<MergeNodeId>{2, 5, 6});
  REQUIRE(f.internal[1].edges == std::vector<EdgeId>{2});
  // Invalid input.
  REQUIRE_THROWS_AS(build_merge_events<T>(1, std::vector<Edge<T>>{{0, 1, T(1)}}),
                    std::out_of_range);
  const std::vector<EdgeId> bad{3};
  REQUIRE_THROWS_AS(build_merge_events<T>(5, std::span<const Edge<T>>(two), bad),
                    std::out_of_range);
}

TEMPLATE_TEST_CASE("merge events use only the candidate edges", "[paper][merge_events]",
                   std::uint64_t, double) {
  using T = TestType;
  // Without the light edge 0 (not a candidate), vertices 0, 1, 2 merge at weight 3 at once.
  const std::vector<Edge<T>> e{{0, 1, T(1)}, {0, 1, T(3)}, {1, 2, T(3)}};
  const std::vector<EdgeId> cand{1, 2};
  const auto f = build_merge_events<T>(3, std::span<const Edge<T>>(e), cand);
  REQUIRE(f.internal.size() == 1);
  REQUIRE(f.internal[0].children == std::vector<MergeNodeId>{0, 1, 2});
  REQUIRE(f.internal[0].edges == std::vector<EdgeId>{1, 2});
  REQUIRE(f.forest_edges == std::vector<EdgeId>{1, 2});
}

TEMPLATE_TEST_CASE("merge events equal the threshold-component reference (Lemma 3.5)",
                   "[paper][merge_events]", std::uint64_t, double) {
  using T = TestType;
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    const std::size_t n = 1 + seed % 17;
    const std::size_t m = (seed * 5) % (3 * n + 1);
    const auto edges = random_multigraph<T>(n, m, 1 + seed % 5, seed);
    INFO("seed " << seed << " n " << n << " m " << m);
    const auto f = build_merge_events<T>(n, edges);
    check_structure(f, n, edges);

    // Exact node set: heights and vertex sets equal the threshold-component reference.
    const auto sets = vertex_sets(f);
    std::vector<std::pair<T, std::vector<Vertex>>> got;
    for (std::size_t i = 0; i < f.internal.size(); ++i) {
      got.emplace_back(f.internal[i].height, sets[n + i]);
    }
    std::sort(got.begin(), got.end());
    REQUIRE(got == reference_nodes(n, edges));

    // The forest edges are the (length, index) minimum spanning forest.
    auto fe = f.forest_edges;
    std::sort(fe.begin(), fe.end());
    REQUIRE(fe == sssp::paper::mst_kruskal_ref<T>(n, edges));

    // Lemma 3.5: no edge lighter than h(C) joins two children of C; in fact the height of the
    // LCA is the minimax value b_T(x, y), and equals +infinity iff x, y are disconnected.
    const auto mm = minimax(n, edges);
    for (std::size_t x = 0; x < n; ++x) {
      for (std::size_t y = x + 1; y < n; ++y) {
        const MergeNodeId a = lca(f, static_cast<MergeNodeId>(x), static_cast<MergeNodeId>(y));
        if (a == kNoMergeNode) {
          REQUIRE(mm[x][y] == std::numeric_limits<T>::max());
        } else {
          REQUIRE(f.internal[a - n].height == mm[x][y]);
        }
      }
    }
    for (const auto& e : edges) {
      if (e.u == e.v) continue;
      const MergeNodeId a = lca(f, e.u, e.v);
      REQUIRE_FALSE(e.weight < f.internal[a - n].height);
    }
  }
}

TEST_CASE("merge events record counters in Diagnostics", "[paper][merge_events]") {
  const std::vector<Edge<std::uint64_t>> e{{0, 1, 1}, {2, 3, 1}, {1, 2, 5}, {0, 3, 5}};
  sssp::Diagnostics diag;
  const auto f = build_merge_events<std::uint64_t>(4, e, &diag);
  REQUIRE(f.internal.size() == 3);
#if SSSP_DIAGNOSTICS
  REQUIRE(diag.counters["merge_events.calls"] == 1);
  REQUIRE(diag.counters["merge_events.nodes"] == 3);
  REQUIRE(diag.counters["merge_events.forest_edges"] == 3);
  REQUIRE(diag.counters["merge_events.weight_classes"] == 2);
  REQUIRE(diag.counters["dsu.unions"] == 3);
#else
  REQUIRE(diag.counters.empty());
#endif
}
