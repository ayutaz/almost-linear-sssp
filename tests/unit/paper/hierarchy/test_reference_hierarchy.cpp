// Tests for the threshold-SCC reference hierarchy (test support; a supplement to the paper,
// research notes 02 sec 2.11). It must satisfy Definition 4.1 with c0 = 1 on every directed
// state graph with positive lengths, and its heights must be exact: the LCA of two states has
// height equal to the least W at which they are mutually reachable using arcs of length <= W.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "sssp/graph.hpp"

using sssp::Vertex;
using sssp::paper::kNoNode;
using sssp::paper::NodeId;
using sssp::paper::OrderedHierarchy;
using sssp::paper::StateArc;
using sssp::paper::StateGraph;
using sssp::paper::StateId;

namespace {

template <class W>
StateGraph<W> random_state_graph(std::mt19937_64& rng, std::size_t n, std::size_t m,
                                 std::uint64_t max_len) {
  StateGraph<W> g;
  g.num_states = n;
  for (std::size_t i = 0; i < m; ++i) {
    const auto u = static_cast<StateId>(rng() % n);
    const auto v = static_cast<StateId>(rng() % n);
    g.arcs.push_back({u, v, static_cast<W>(1 + rng() % max_len)});
  }
  return g;
}

// reach[x][y]: y reachable from x using arcs of length <= w (naive closure, test oracle).
template <class W>
std::vector<std::vector<bool>> closure(const StateGraph<W>& g, W w) {
  const std::size_t n = g.num_states;
  std::vector<std::vector<bool>> r(n, std::vector<bool>(n, false));
  for (std::size_t i = 0; i < n; ++i) r[i][i] = true;
  for (const auto& a : g.arcs) {
    if (a.length <= w) r[a.tail][a.head] = true;
  }
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t i = 0; i < n; ++i)
      if (r[i][k])
        for (std::size_t j = 0; j < n; ++j)
          if (r[k][j]) r[i][j] = true;
  return r;
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

}  // namespace

TEMPLATE_TEST_CASE("reference hierarchy of a hand-built cycle has one node at the cycle's max",
                   "[paper][hierarchy][reference]", std::uint64_t, double) {
  using W = TestType;
  StateGraph<W> g;
  g.num_states = 4;
  // 0 -> 1 (1), 1 -> 2 (1), 2 -> 0 (2) form an SCC at W = 2; 3 is reached but never returns.
  g.arcs = {{0, 1, W{1}}, {1, 2, W{1}}, {2, 0, W{2}}, {2, 3, W{1}}};
  const auto h = sssp::paper::build_threshold_scc_hierarchy<W>(g, {0, 1, 2, 3});
  const auto phys = sssp::paper::physical_graph_of_arcs<W>(g, {0, 1, 2, 3}, 4);
  const auto rep = sssp::paper::validate_ordered_hierarchy(h, phys);
  INFO(rep.summary());
  REQUIRE(rep.ok());
  REQUIRE(h.c0 == 1);
  REQUIRE(h.nodes.size() == 5);  // four leaves + one internal node
  REQUIRE(h.roots.size() == 2);  // {0,1,2} and {3}
  const auto& c = h.nodes[4];
  REQUIRE(c.height == W{2});
  REQUIRE(c.children.size() == 3);
  REQUIRE(c.links.size() == 2);
  for (const auto& lk : c.links) REQUIRE(lk.length <= c.height);
  REQUIRE(h.nodes[3].parent == kNoNode);
}

TEST_CASE("reference hierarchy of a state graph without arcs is a forest of leaves",
          "[paper][hierarchy][reference]") {
  StateGraph<std::uint64_t> g;
  g.num_states = 3;
  const auto h = sssp::paper::build_threshold_scc_hierarchy<std::uint64_t>(g, {0, 0, 1});
  REQUIRE(h.nodes.size() == 3);
  REQUIRE(h.roots.size() == 3);
  const auto phys = sssp::paper::physical_graph_of_arcs<std::uint64_t>(g, {0, 0, 1}, 2);
  REQUIRE(sssp::paper::validate_ordered_hierarchy(h, phys).ok());
}

TEST_CASE("reference hierarchy rejects lengths that are not positive and finite",
          "[paper][hierarchy][reference]") {
  StateGraph<std::uint64_t> g;
  g.num_states = 2;
  g.arcs = {{0, 1, 0}, {1, 0, 3}};
  REQUIRE_THROWS_AS(sssp::paper::build_threshold_scc_hierarchy<std::uint64_t>(g, {0, 1}),
                    std::invalid_argument);
  StateGraph<double> gd;
  gd.num_states = 2;
  gd.arcs = {{0, 1, -1.0}};
  REQUIRE_THROWS_AS(sssp::paper::build_threshold_scc_hierarchy<double>(gd, {0, 1}),
                    std::invalid_argument);
  gd.arcs = {{0, 1, std::numeric_limits<double>::quiet_NaN()}};
  REQUIRE_THROWS_AS(sssp::paper::build_threshold_scc_hierarchy<double>(gd, {0, 1}),
                    std::invalid_argument);
  StateGraph<std::uint64_t> bad;
  bad.num_states = 2;
  bad.arcs = {{0, 5, 1}};
  REQUIRE_THROWS_AS(sssp::paper::build_threshold_scc_hierarchy<std::uint64_t>(bad, {0, 1}),
                    std::out_of_range);
  REQUIRE_THROWS_AS(sssp::paper::build_threshold_scc_hierarchy<std::uint64_t>(g, {0}),
                    std::invalid_argument);
}

TEMPLATE_TEST_CASE("reference hierarchy is valid, strictly monotone and exact on random graphs",
                   "[paper][hierarchy][reference][random]", std::uint64_t, double) {
  using W = TestType;
  std::size_t internal_nodes = 0;
  std::size_t light_arcs = 0;
  std::size_t multi_root = 0;
  for (std::uint64_t seed = 1; seed <= 300; ++seed) {
    std::mt19937_64 rng(seed);
    const std::size_t n = 1 + rng() % 12;
    const std::size_t m = n + rng() % (3 * n + 1);
    const std::uint64_t max_len = 1 + rng() % 6;  // small range: many equal lengths
    const auto g = random_state_graph<W>(rng, n, m, max_len);
    const std::size_t n_phys = 1 + rng() % n;
    std::vector<Vertex> pi(n);
    for (auto& p : pi) p = static_cast<Vertex>(rng() % n_phys);
    const auto h = sssp::paper::build_threshold_scc_hierarchy<W>(g, pi);
    const auto phys = sssp::paper::physical_graph_of_arcs<W>(g, pi, n_phys);
    INFO("seed " << seed << " n " << n << " m " << m);

    const auto rep = sssp::paper::validate_ordered_hierarchy(h, phys);
    INFO(rep.summary());
    REQUIRE(rep.ok());
    REQUIRE(h.c0 == 1);

    // Heights strictly increase towards the root (one node per merge threshold).
    for (const auto& node : h.nodes) {
      if (node.children.empty() || node.parent == kNoNode) continue;
      REQUIRE(node.height < h.nodes[node.parent].height);
    }
    for (const auto& node : h.nodes) internal_nodes += node.children.empty() ? 0 : 1;
    if (h.roots.size() > 1) ++multi_root;

    // Exactness against the naive closure: h(LCA(x, y)) is the least arc length W at which x
    // and y are mutually reachable with arcs <= W; different trees iff never mutually reachable.
    std::vector<W> lengths;
    for (const auto& a : g.arcs) lengths.push_back(a.length);
    std::sort(lengths.begin(), lengths.end());
    lengths.erase(std::unique(lengths.begin(), lengths.end()), lengths.end());
    std::vector<std::vector<std::vector<bool>>> reach;
    for (const W w : lengths) reach.push_back(closure(g, w));
    for (StateId x = 0; x < n; ++x) {
      for (StateId y = x + 1; y < n; ++y) {
        std::size_t first = lengths.size();
        for (std::size_t k = 0; k < lengths.size(); ++k) {
          if (reach[k][x][y] && reach[k][y][x]) {
            first = k;
            break;
          }
        }
        const NodeId l = lca(h, h.leaf_of_state[x], h.leaf_of_state[y]);
        if (first == lengths.size()) {
          REQUIRE(l == kNoNode);
        } else {
          REQUIRE(l != kNoNode);
          REQUIRE(h.nodes[l].height == lengths[first]);
        }
      }
    }
    // Light arcs (length < h(LCA)) actually occur, so Evaluate tests on these are non-trivial.
    for (const auto& a : g.arcs) {
      if (a.tail == a.head) continue;
      const NodeId l = lca(h, h.leaf_of_state[a.tail], h.leaf_of_state[a.head]);
      if (l != kNoNode && a.length < h.nodes[l].height) ++light_arcs;
    }
  }
  // The random family is not vacuous.
  REQUIRE(internal_nodes >= 300);
  REQUIRE(light_arcs >= 100);
  REQUIRE(multi_root >= 50);
}
