// Changing the heights of the threshold-SCC reference hierarchy (research notes 02 sec 2.11
// item 1, "scaled heights") against Definition 4.1 (Sec 4, l.646-658).
//
// Expected behaviour, derived from Def 4.1 and not from the checker:
//  * A node C of the reference hierarchy with merge threshold W(C) is a strongly connected
//    component of the arcs of length <= W(C); its >= 2 children are the components of the arcs
//    of length < W(C) inside it. All paths that make C strongly connected stay inside C.
//  * Raising h(C) to any h' > W(C) makes every arc of length <= W(C) between distinct children
//    light (length strictly below h', l.650), so the light quotient on the children is strongly
//    connected on >= 2 vertices: H3 fails at C. In particular the "h(C) = 2W" variant fails H3
//    at EVERY internal node; it is not a valid ordered hierarchy.
//  * Lowering heights keeps H3 (the light arcs only shrink), H2 needs h(child) <= h(C) and
//    h > 0 (l.649), and H5 needs every link <= c0 * h(C) (l.652). So any heights with
//    0 < h'(child) <= h'(C) <= W(C) are valid with c0 = max over C of ceil(max link / h'(C));
//    halving every height is valid with c0 = 2 (heights scaled by lambda <= 1 need c0/lambda).
// The H3 expectation is also recomputed here with a naive quotient-closure oracle.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <set>
#include <type_traits>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "sssp/graph.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
using sssp::paper::Condition;
using sssp::paper::kNoNode;
using sssp::paper::NodeId;
using sssp::paper::OrderedHierarchy;
using sssp::paper::StateGraph;
using sssp::paper::StateId;
using sssp::paper::ViolationKind;

namespace {

struct Instance {
  std::size_t n = 0;
  std::size_t n_phys = 0;
};

// Random directed state graph with lengths step * (1 .. max_len) (step = 2 makes halving
// exact for integer weights), and random images.
template <class W>
StateGraph<W> random_graph(std::mt19937_64& rng, std::uint64_t step, std::vector<Vertex>& pi,
                           Instance& inst) {
  inst.n = 2 + rng() % 11;
  const std::size_t m = inst.n + rng() % (3 * inst.n + 1);
  const std::uint64_t max_len = 1 + rng() % 6;
  StateGraph<W> g;
  g.num_states = inst.n;
  for (std::size_t i = 0; i < m; ++i) {
    g.arcs.push_back({static_cast<StateId>(rng() % inst.n), static_cast<StateId>(rng() % inst.n),
                      static_cast<W>(step * (1 + rng() % max_len))});
  }
  inst.n_phys = 1 + rng() % inst.n;
  pi.assign(inst.n, 0);
  for (auto& p : pi) p = static_cast<Vertex>(rng() % inst.n_phys);
  return g;
}

template <class W>
CsrGraph<W> own_physical_graph(const StateGraph<W>& g, const std::vector<Vertex>& pi,
                               std::size_t n_phys) {
  std::vector<Edge<W>> edges;
  for (const auto& a : g.arcs) edges.push_back({pi[a.tail], pi[a.head], a.length});
  return CsrGraph<W>::from_edges(n_phys, edges);
}

// Naive H3 oracle at node c with height `height`: quotient of the children of c on the arcs
// between distinct children of length < height; a cycle exists iff some child reaches itself.
template <class W>
bool naive_light_cycle(const OrderedHierarchy<W>& h, NodeId c, W height) {
  const auto& kids = h.nodes[c].children;
  const std::size_t k = kids.size();
  auto child_of = [&](StateId x) -> std::size_t {
    NodeId prev = kNoNode;
    for (NodeId v = h.leaf_of_state[x]; v != kNoNode; prev = v, v = h.nodes[v].parent) {
      if (v == c)
        return static_cast<std::size_t>(std::find(kids.begin(), kids.end(), prev) - kids.begin());
    }
    return k;  // not below c
  };
  std::vector<std::vector<bool>> r(k, std::vector<bool>(k, false));
  for (const auto& a : h.graph.arcs) {
    const std::size_t u = child_of(a.tail);
    const std::size_t v = child_of(a.head);
    if (u < k && v < k && u != v && a.length < height) r[u][v] = true;
  }
  for (std::size_t m = 0; m < k; ++m)
    for (std::size_t i = 0; i < k; ++i)
      if (r[i][m])
        for (std::size_t j = 0; j < k; ++j)
          if (r[m][j]) r[i][j] = true;
  for (std::size_t i = 0; i < k; ++i) {
    if (r[i][i]) return true;
  }
  return false;
}

template <class W>
std::vector<NodeId> internal_nodes(const OrderedHierarchy<W>& h) {
  std::vector<NodeId> out;
  for (NodeId v = 0; v < h.nodes.size(); ++v) {
    if (!h.nodes[v].children.empty()) out.push_back(v);
  }
  return out;
}

// Smallest value strictly above w (exact for the integral values used here).
template <class W>
W just_above(W w) {
  if constexpr (std::is_floating_point_v<W>) {
    return std::nextafter(w, std::numeric_limits<W>::infinity());
  } else {
    return w + 1;
  }
}

}  // namespace

TEMPLATE_TEST_CASE("reference heights: raising one height above its threshold violates H3 there",
                   "[paper][hierarchy][reference][heights][def41]", std::uint64_t, double) {
  using W = TestType;
  std::size_t checked = 0;
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    std::mt19937_64 rng(seed);
    std::vector<Vertex> pi;
    Instance inst;
    const auto g = random_graph<W>(rng, 1, pi, inst);
    const auto base = sssp::paper::build_threshold_scc_hierarchy<W>(g, pi);
    const auto phys = own_physical_graph(g, pi, inst.n_phys);
    INFO("seed " << seed << " n " << inst.n);
    REQUIRE(sssp::paper::validate_ordered_hierarchy(base, phys).ok());
    for (const NodeId c : internal_nodes(base)) {
      const W w = base.nodes[c].height;
      // The unchanged height is the largest valid one: no light cycle at W itself.
      REQUIRE_FALSE(naive_light_cycle(base, c, w));
      for (const W raised : {just_above(w), static_cast<W>(2 * w)}) {
        const NodeId p = base.nodes[c].parent;
        if (p != kNoNode && base.nodes[p].height < raised) continue;  // would also break H2
        auto h = base;
        h.nodes[c].height = raised;
        REQUIRE(naive_light_cycle(h, c, raised));  // the paper-level claim, independently
        const auto rep = sssp::paper::validate_ordered_hierarchy(h, phys);
        INFO("node " << c << " W " << w << " raised to " << raised << "\n" << rep.summary());
        REQUIRE(rep.violations.size() == 1);
        REQUIRE(rep.violations[0].kind == ViolationKind::kLightArcCycle);
        REQUIRE(rep.violations[0].node == c);
        ++checked;
      }
    }
  }
  REQUIRE(checked >= 300);
}

TEMPLATE_TEST_CASE("reference heights: the doubled-height (2W) variant violates H3 at every node",
                   "[paper][hierarchy][reference][heights][def41]", std::uint64_t, double) {
  using W = TestType;
  std::size_t nonempty = 0;
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    std::mt19937_64 rng(seed);
    std::vector<Vertex> pi;
    Instance inst;
    const auto g = random_graph<W>(rng, 1, pi, inst);
    auto h = sssp::paper::build_threshold_scc_hierarchy<W>(g, pi);
    const auto phys = own_physical_graph(g, pi, inst.n_phys);
    const auto nodes = internal_nodes(h);
    for (const NodeId c : nodes) h.nodes[c].height = static_cast<W>(2 * h.nodes[c].height);
    const auto rep = sssp::paper::validate_ordered_hierarchy(h, phys);
    INFO("seed " << seed << " n " << inst.n << "\n" << rep.summary());
    std::set<NodeId> reported;
    for (const auto& v : rep.violations) {
      REQUIRE(v.kind == ViolationKind::kLightArcCycle);
      reported.insert(v.node);
    }
    REQUIRE(reported == std::set<NodeId>(nodes.begin(), nodes.end()));
    for (const NodeId c : nodes) REQUIRE(naive_light_cycle(h, c, h.nodes[c].height));
    if (!nodes.empty()) ++nonempty;
  }
  REQUIRE(nonempty >= 150);
}

TEMPLATE_TEST_CASE("reference heights: halving every height keeps Def 4.1 with c0 = 2",
                   "[paper][hierarchy][reference][heights][def41]", std::uint64_t, double) {
  using W = TestType;
  std::size_t nonempty = 0;
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    std::mt19937_64 rng(seed);
    std::vector<Vertex> pi;
    Instance inst;
    const auto g = random_graph<W>(rng, 2, pi, inst);  // even lengths: W / 2 is exact
    auto h = sssp::paper::build_threshold_scc_hierarchy<W>(g, pi);
    const auto phys = own_physical_graph(g, pi, inst.n_phys);
    REQUIRE(h.c0 == 1);
    const auto nodes = internal_nodes(h);
    for (const NodeId c : nodes) h.nodes[c].height = h.nodes[c].height / 2;
    h.c0 = 2;
    const auto rep = sssp::paper::validate_ordered_hierarchy(h, phys);
    INFO("seed " << seed << " n " << inst.n << "\n" << rep.summary());
    REQUIRE(rep.ok());
    for (const NodeId c : nodes) REQUIRE_FALSE(naive_light_cycle(h, c, h.nodes[c].height));
    if (!nodes.empty()) ++nonempty;
  }
  REQUIRE(nonempty >= 150);
}

TEMPLATE_TEST_CASE("reference heights: any monotone heights at most the thresholds are valid",
                   "[paper][hierarchy][reference][heights][def41]", std::uint64_t, double) {
  using W = TestType;
  std::size_t lowered = 0;
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    std::mt19937_64 rng(seed);
    std::vector<Vertex> pi;
    Instance inst;
    const auto g = random_graph<W>(rng, 1, pi, inst);
    auto h = sssp::paper::build_threshold_scc_hierarchy<W>(g, pi);
    const auto phys = own_physical_graph(g, pi, inst.n_phys);
    // Bottom-up: a node is processed after all its children (children have smaller heights,
    // so sorting by the original height is a valid order).
    auto nodes = internal_nodes(h);
    std::sort(nodes.begin(), nodes.end(),
              [&](NodeId a, NodeId b) { return h.nodes[a].height < h.nodes[b].height; });
    std::uint64_t c0 = 1;
    for (const NodeId c : nodes) {
      std::uint64_t lo = 1;
      for (const NodeId ch : h.nodes[c].children) {
        if (!h.nodes[ch].children.empty()) {
          lo = std::max(lo, static_cast<std::uint64_t>(h.nodes[ch].height));
        }
      }
      const auto hi = static_cast<std::uint64_t>(h.nodes[c].height);
      REQUIRE(lo <= hi);
      const std::uint64_t pick = std::uniform_int_distribution<std::uint64_t>(lo, hi)(rng);
      if (pick < hi) ++lowered;
      h.nodes[c].height = static_cast<W>(pick);
      for (const auto& lk : h.nodes[c].links) {
        const auto len = static_cast<std::uint64_t>(lk.length);
        c0 = std::max(c0, (len + pick - 1) / pick);
      }
    }
    h.c0 = static_cast<std::uint32_t>(c0);
    const auto rep = sssp::paper::validate_ordered_hierarchy(h, phys);
    INFO("seed " << seed << " n " << inst.n << " c0 " << c0 << "\n" << rep.summary());
    REQUIRE(rep.ok());
  }
  REQUIRE(lowered >= 100);
}
