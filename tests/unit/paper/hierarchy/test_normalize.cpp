// Tests for normalisation (Sec 4.1, l.674-723, Lemma 4.3).
//
// Every expectation is derived from the paper's text, recomputed here by independent, naive
// oracles (recursive masses, a brute-force power-of-two table per link set, naive LCA,
// Floyd-Warshall physical distances); nothing is copied from an implementation run.
//
// Inputs: the threshold-SCC reference hierarchy on random directed state graphs with random
// many-to-one projections, and the Alg 4 base case on generated physical graphs; several
// seeds, for std::uint64_t and double.
//
// Tolerance (double only, MassBoundsPhysicalDiameter): S(C) and the Floyd-Warshall distance
// are both rounded sums of at most ~2N nonnegative terms in different orders, so each carries a
// relative error below 2N * 2^-53 < 1e-14 for N <= 64. The check allows a relative slack of
// 1e-12. Integer weights and every structural property are compared exactly.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "paper/local/base_case.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
using sssp::paper::ArcId;
using sssp::paper::kNoNode;
using sssp::paper::NodeId;
using sssp::paper::NormalizedHierarchy;
using sssp::paper::OrderedHierarchy;
using sssp::paper::PhysicalLink;
using sssp::paper::StateArc;
using sssp::paper::StateGraph;
using sssp::paper::StateId;

namespace {

template <class W>
constexpr W kInf = sssp::WeightTraits<W>::infinity();

// ------------------------------------------------------------------ input generators

// Random directed state graph with positive lengths (self-loops and parallel arcs occur). For
// double, lengths are random reals in [0.5, 16.5) so that sums are genuinely rounded.
template <class W>
StateGraph<W> random_state_graph(std::mt19937_64& rng, std::size_t n, std::size_t m) {
  StateGraph<W> g;
  g.num_states = n;
  for (std::size_t i = 0; i < m; ++i) {
    const auto u = static_cast<StateId>(rng() % n);
    const auto v = static_cast<StateId>(rng() % n);
    W len{};
    if constexpr (std::is_integral_v<W>) {
      len = static_cast<W>(1 + rng() % 12);
    } else {
      len = 0.5 + std::uniform_real_distribution<double>(0.0, 16.0)(rng);
    }
    g.arcs.push_back({u, v, len});
  }
  return g;
}

struct RefInput {
  std::size_t states;
  std::size_t arcs;
  std::size_t physical;
};

template <class W>
struct Instance {
  OrderedHierarchy<W> h;
  CsrGraph<W> phys;
};

template <class W>
Instance<W> reference_instance(std::uint64_t seed, RefInput in) {
  std::mt19937_64 rng(seed);
  auto g = random_state_graph<W>(rng, in.states, in.arcs);
  std::vector<Vertex> pi(in.states);
  for (auto& p : pi) p = static_cast<Vertex>(rng() % in.physical);
  auto phys = sssp::paper::physical_graph_of_arcs<W>(g, pi, in.physical);
  auto h = sssp::paper::build_threshold_scc_hierarchy<W>(std::move(g), std::move(pi));
  return {std::move(h), std::move(phys)};
}

template <class W>
Instance<W> base_case_instance(std::uint64_t seed, int family) {
  sssp::gen::GenOptions opt;
  opt.seed = seed;
  opt.connected = (seed % 3) != 0;  // some inputs stay disconnected (several roots)
  if constexpr (std::is_integral_v<W>) {
    opt.weights = {sssp::gen::WeightDist::kUniformInt, 1.0, 9.0};
  } else {
    opt.weights = {sssp::gen::WeightDist::kUniformReal, 0.5, 9.0};
  }
  sssp::EdgeList<W> el;
  switch (family) {
    case 0:
      el = sssp::gen::erdos_renyi<W>(10 + (seed % 25), 3.0, opt);
      break;
    case 1:
      el = sssp::gen::random_tree<W>(5 + (seed % 30), opt);
      break;
    default:
      el = sssp::gen::grid_2d<W>(2 + (seed % 4), 2 + (seed % 5), opt);
      break;
  }
  auto g = sssp::to_csr(std::move(el));
  auto r = sssp::paper::local_hierarchy_base_case(g);
  return {std::move(r.hierarchy), std::move(g)};
}

// ------------------------------------------------------------------ oracles

template <class W>
std::vector<std::vector<W>> floyd_warshall(const CsrGraph<W>& g) {
  const std::size_t n = g.num_vertices();
  std::vector<std::vector<W>> d(n, std::vector<W>(n, kInf<W>));
  for (std::size_t i = 0; i < n; ++i) d[i][i] = W{0};
  for (const auto& e : g.edges()) {
    if (e.weight < d[e.u][e.v]) d[e.u][e.v] = d[e.v][e.u] = e.weight;
  }
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        if (d[i][k] != kInf<W> && d[k][j] != kInf<W> && d[i][k] + d[k][j] < d[i][j])
          d[i][j] = d[i][k] + d[k][j];
  return d;
}

// Sec 4.1 l.686-688: S(C) = total link lengths + heights of internal nodes in C's subtree;
// leaves have mass zero. Plain recursion over the INPUT hierarchy.
template <class W>
W oracle_mass(const OrderedHierarchy<W>& h, NodeId v) {
  const auto& nd = h.nodes[v];
  if (nd.children.empty()) return W{0};
  W s = nd.height;
  for (const auto& lk : nd.links) s += lk.length;
  for (const NodeId c : nd.children) s += oracle_mass(h, c);
  return s;
}

template <class W>
std::vector<W> oracle_masses(const OrderedHierarchy<W>& h) {
  std::vector<W> s(h.nodes.size());
  for (NodeId v = 0; v < h.nodes.size(); ++v) s[v] = oracle_mass(h, v);
  return s;
}

// Sec 4.1 l.689: an internal child C is joined to its parent P iff 4 S(C) >= h(P).
template <class W>
std::vector<bool> oracle_linked(const OrderedHierarchy<W>& h, const std::vector<W>& s) {
  std::vector<bool> linked(h.nodes.size(), false);
  for (NodeId v = 0; v < h.nodes.size(); ++v) {
    const auto& nd = h.nodes[v];
    if (nd.children.empty() || nd.parent == kNoNode) continue;
    linked[v] = !(s[v] + s[v] + s[v] + s[v] < h.nodes[nd.parent].height);
  }
  return linked;
}

// Sec 4.1 l.690-692: within each connected set of normalisation links, a table of powers of two
// starts at the least height of the set, and nu(C) is h(C) rounded down on that table. A node
// without links is a one-node set (Q-24), so its table starts at h(C) and nu(C) = h(C).
template <class W>
std::vector<W> oracle_units(const OrderedHierarchy<W>& h, const std::vector<bool>& linked) {
  const std::size_t n = h.nodes.size();
  std::vector<NodeId> top(n, kNoNode);
  for (NodeId v = 0; v < n; ++v) {
    if (h.nodes[v].children.empty()) continue;
    NodeId t = v;
    while (linked[t]) t = h.nodes[t].parent;
    top[v] = t;
  }
  std::vector<W> unit(n, W{0});
  for (NodeId v = 0; v < n; ++v) {
    if (top[v] == kNoNode) continue;
    W hmin = h.nodes[v].height;
    for (NodeId u = 0; u < n; ++u) {
      if (top[u] == top[v]) hmin = std::min(hmin, h.nodes[u].height);
    }
    W t = hmin;
    while (!(h.nodes[v].height < t + t)) t = t + t;
    unit[v] = t;
  }
  return unit;
}

template <class W>
NodeId naive_lca(const OrderedHierarchy<W>& h, NodeId a, NodeId b) {
  std::vector<NodeId> up;
  for (NodeId v = a; v != kNoNode; v = h.nodes[v].parent) up.push_back(v);
  for (NodeId v = b; v != kNoNode; v = h.nodes[v].parent) {
    if (std::find(up.begin(), up.end(), v) != up.end()) return v;
  }
  return kNoNode;
}

template <class W>
NodeId child_toward(const OrderedHierarchy<W>& h, NodeId anc, NodeId v) {
  while (h.nodes[v].parent != anc) v = h.nodes[v].parent;
  return v;
}

template <class W>
NodeId root_of(const OrderedHierarchy<W>& h, NodeId v) {
  while (h.nodes[v].parent != kNoNode) v = h.nodes[v].parent;
  return v;
}

template <class W>
void states_below(const OrderedHierarchy<W>& h, NodeId v, std::vector<StateId>& out) {
  if (h.nodes[v].children.empty()) {
    out.push_back(h.nodes[v].state);
    return;
  }
  for (const NodeId c : h.nodes[v].children) states_below(h, c, out);
}

// True iff a <= b * 2^k for some k >= 0 with a == b * 2^k (doublings only).
template <class W>
bool is_power_of_two_multiple(W a, W b) {
  W t = b;
  while (t < a) t = t + t;
  return t == a;
}

// Kahn's algorithm on a small quotient: true iff acyclic.
bool acyclic(std::size_t k, const std::vector<std::pair<std::size_t, std::size_t>>& arcs) {
  std::vector<std::size_t> indeg(k, 0);
  std::vector<std::vector<std::size_t>> out(k);
  for (const auto& [a, b] : arcs) {
    out[a].push_back(b);
    ++indeg[b];
  }
  std::vector<std::size_t> q;
  for (std::size_t i = 0; i < k; ++i)
    if (indeg[i] == 0) q.push_back(i);
  std::size_t seen = 0;
  while (!q.empty()) {
    const std::size_t x = q.back();
    q.pop_back();
    ++seen;
    for (const std::size_t y : out[x])
      if (--indeg[y] == 0) q.push_back(y);
  }
  return seen == k;
}

template <class W>
std::size_t internal_count(const OrderedHierarchy<W>& h) {
  return static_cast<std::size_t>(std::count_if(h.nodes.begin(), h.nodes.end(),
                                                [](const auto& v) { return !v.children.empty(); }));
}

template <class W>
std::size_t link_count(const OrderedHierarchy<W>& h) {
  std::size_t k = 0;
  for (const auto& v : h.nodes) k += v.links.size();
  return k;
}

// ------------------------------------------------------------------ the property checks

// Checks every claim of Sec 4.1 / Lemma 4.3 on `nh` and returns sum_C S(C)/nu(C) / N.
template <class W>
double check_normalized(const OrderedHierarchy<W>& in, const CsrGraph<W>& phys,
                        const NormalizedHierarchy<W>& nh) {
  const auto& f = nh.forest;
  const std::size_t nf = f.nodes.size();
  REQUIRE(nh.unit.size() == nf);
  REQUIRE(nh.mass.size() == nf);
  REQUIRE(nh.original_height.size() == nf);
  REQUIRE(nh.original_node.size() == nf);
  REQUIRE(nh.node_of_input.size() == in.nodes.size());
  REQUIRE(nh.light_begin.size() == nf + 1);
  REQUIRE(f.num_states() == in.num_states());

  const auto s = oracle_masses(in);
  const auto linked = oracle_linked(in, s);
  const auto nu = oracle_units(in, linked);

  // Contraction (l.710-713): exactly the normalisation links with equal units disappear.
  std::size_t contracted = 0;
  for (NodeId v = 0; v < in.nodes.size(); ++v) {
    if (linked[v] && nu[v] == nu[in.nodes[v].parent]) ++contracted;
  }
  REQUIRE(nf == in.nodes.size() - contracted);

  for (NodeId d = 0; d < in.nodes.size(); ++d) {
    const NodeId x = nh.node_of_input[d];
    REQUIRE(x < nf);
    if (in.nodes[d].children.empty()) {
      REQUIRE(f.nodes[x].children.empty());
      REQUIRE(f.nodes[x].state == in.nodes[d].state);
      REQUIRE(nh.mass[x] == W{0});
      continue;
    }
    // Units follow the power-of-two table of the link set and lie in (h/2, h] (l.690-691),
    // with the pre-contraction height of every original node.
    REQUIRE(nh.unit[x] == nu[d]);
    REQUIRE(nh.unit[x] <= in.nodes[d].height);
    REQUIRE(in.nodes[d].height < nh.unit[x] + nh.unit[x]);
    const bool any_link =
        linked[d] || std::any_of(in.nodes[d].children.begin(), in.nodes[d].children.end(),
                                 [&](NodeId c) { return linked[c]; });
    if (!any_link) REQUIRE(nh.unit[x] == in.nodes[d].height);  // Q-24
    // A node is absorbed exactly when its normalisation link joins equal units.
    const bool absorbed = linked[d] && nu[d] == nu[in.nodes[d].parent];
    REQUIRE((nh.original_node[x] == d) == !absorbed);
    if (absorbed) REQUIRE(x == nh.node_of_input[in.nodes[d].parent]);
  }
  for (NodeId x = 0; x < nf; ++x) {
    const NodeId d = nh.original_node[x];
    REQUIRE(nh.node_of_input[d] == x);
    if (f.nodes[x].children.empty()) continue;
    // The topmost original mass and height are kept (l.713).
    REQUIRE(nh.mass[x] == s[d]);
    REQUIRE(nh.original_height[x] == in.nodes[d].height);
    // The forest's own height is the unit (header contract).
    REQUIRE(f.nodes[x].height == nh.unit[x]);
  }

  // Lemma 4.3: units at least double; dyadic alignment or 2 S(C) < nu(P).
  for (NodeId x = 0; x < nf; ++x) {
    const auto& nd = f.nodes[x];
    if (nd.children.empty() || nd.parent == kNoNode) continue;
    const W cu = nh.unit[x];
    const W pu = nh.unit[nd.parent];
    REQUIRE(cu + cu <= pu);  // UnitsAtLeastDouble
    const bool dyadic = is_power_of_two_multiple(pu, cu);
    const bool separated = nh.mass[x] + nh.mass[x] < pu;
    REQUIRE((dyadic || separated));  // DyadicOrMassSeparated
  }

  // MassBoundsPhysicalDiameter: S(C) >= every physical distance between states below C.
  const auto dist = floyd_warshall(phys);
  for (NodeId x = 0; x < nf; ++x) {
    if (f.nodes[x].children.empty()) continue;
    std::vector<StateId> st;
    states_below(f, x, st);
    for (const StateId a : st) {
      for (const StateId b : st) {
        const W dab = dist[f.projection[a]][f.projection[b]];
        REQUIRE(dab != kInf<W>);
        if constexpr (std::is_integral_v<W>) {
          REQUIRE(dab <= nh.mass[x]);
        } else {
          REQUIRE(dab <= nh.mass[x] * (1.0 + 1e-12));
        }
      }
    }
  }

  // MassOverUnitLinear: each mass item z (a height h(D) or a link of D, z <= max(1,c0) h(D))
  // is charged to the surviving ancestors of D, whose units at least double from nu(D) >
  // h(D)/2 upwards, so it contributes < 4 max(1,c0). Hence sum S/nu < 4 max(1,c0) (#internal
  // + #links) of the input.
  double sum = 0.0;
  for (NodeId x = 0; x < nf; ++x) {
    if (f.nodes[x].children.empty()) continue;
    sum += static_cast<double>(nh.mass[x]) / static_cast<double>(nh.unit[x]);
  }
  const double c0 = std::max<double>(1.0, in.c0);
  REQUIRE(sum < 4.0 * c0 * static_cast<double>(internal_count(in) + link_count(in)) + 1e-9);

  // Arc assignment (l.719-722) against a naive LCA on the normalised forest.
  std::vector<std::vector<std::pair<ArcId, std::pair<NodeId, NodeId>>>> expect(nf);
  std::vector<ArcId> expect_root;
  for (ArcId a = 0; a < f.graph.arcs.size(); ++a) {
    const auto& arc = f.graph.arcs[a];
    if (arc.tail == arc.head) continue;  // self-arcs are discarded
    const NodeId lx = f.leaf_of_state[arc.tail];
    const NodeId ly = f.leaf_of_state[arc.head];
    if (root_of(f, lx) != root_of(f, ly)) {
      expect_root.push_back(a);
      continue;
    }
    const NodeId l = naive_lca(f, lx, ly);
    if (arc.length < nh.unit[l]) {
      expect[l].push_back({a, {child_toward(f, l, lx), child_toward(f, l, ly)}});
    }
  }
  std::vector<ArcId> got_root = nh.root_arcs;
  std::sort(got_root.begin(), got_root.end());
  REQUIRE(got_root == expect_root);
  for (NodeId x = 0; x < nf; ++x) {
    std::vector<std::pair<ArcId, std::pair<NodeId, NodeId>>> got;
    for (const auto& la : nh.light(x)) got.push_back({la.arc, {la.tail_child, la.head_child}});
    std::sort(got.begin(), got.end());
    REQUIRE(got == expect[x]);
    // LightQuotientAcyclicAfterContraction (proof of Lemma 4.3, l.711-712).
    const auto& ch = f.nodes[x].children;
    std::vector<std::pair<std::size_t, std::size_t>> q;
    for (const auto& [arc, ends] : got) {
      const auto pos = [&](NodeId c) {
        return static_cast<std::size_t>(std::find(ch.begin(), ch.end(), c) - ch.begin());
      };
      REQUIRE(pos(ends.first) < ch.size());
      REQUIRE(pos(ends.second) < ch.size());
      q.emplace_back(pos(ends.first), pos(ends.second));
    }
    REQUIRE(acyclic(ch.size(), q));
  }

  // Root order: a permutation of the roots in which every inter-root arc goes forward.
  std::vector<NodeId> ro = nh.root_order;
  std::vector<NodeId> roots = f.roots;
  std::sort(ro.begin(), ro.end());
  std::sort(roots.begin(), roots.end());
  REQUIRE(ro == roots);
  std::vector<std::size_t> pos(nf, 0);
  for (std::size_t i = 0; i < nh.root_order.size(); ++i) pos[nh.root_order[i]] = i;
  for (const ArcId a : nh.root_arcs) {
    const auto& arc = f.graph.arcs[a];
    REQUIRE(pos[root_of(f, f.leaf_of_state[arc.tail])] <
            pos[root_of(f, f.leaf_of_state[arc.head])]);
  }

  // The normalised forest is itself an ordered hierarchy (Def 4.1) with heights nu and
  // connection constant 2 c0 (header contract): links <= c0 h(D) < 2 c0 nu(D).
  REQUIRE(f.c0 == 2 * in.c0);
  const auto rep = sssp::paper::validate_ordered_hierarchy(f, phys);
  INFO(rep.summary());
  REQUIRE(rep.ok());
  return sum / static_cast<double>(f.num_states());
}

}  // namespace

TEMPLATE_TEST_CASE("normalize: hand examples of contraction, dyadic alignment and separation",
                   "[paper][hierarchy][normalize]", std::uint64_t, double) {
  using W = TestType;
  // Three states on a physical path 0 - 1 - 2 with edge lengths hc and hp; the state arcs are
  // both orientations of these edges. C = {0, 1} has height hc and link (0, 1, hc); its parent
  // P = {C, 2} has height hp and link (1, 2, hp). S(C) = 2 hc. No arc is lighter than the LCA
  // height, so the input satisfies Def 4.1 with c0 = 1 (checked).
  NodeId c = kNoNode;
  NodeId p = kNoNode;
  const auto make = [&](W hc, W hp) {
    StateGraph<W> g;
    g.num_states = 3;
    g.arcs = {{0, 1, hc}, {1, 0, hc}, {1, 2, hp}, {2, 1, hp}};
    auto h = sssp::paper::make_leaf_forest<W>(g, {0, 1, 2});
    c = sssp::paper::add_internal_node<W>(h, hc, {0, 1}, {{0, 1, hc}});
    p = sssp::paper::add_internal_node<W>(h, hp, {c, 2}, {{1, 2, hp}});
    sssp::paper::collect_roots(h);
    return h;
  };
  const auto phys_of = [](W hc, W hp) {
    return CsrGraph<W>::from_edges(3, {{0, 1, hc}, {1, 2, hp}});
  };

  SECTION("4S(C) >= h(P) and equal units: the child is contracted into the parent") {
    // hc = 3, hp = 5: 4 S(C) = 24 >= 5, linked. The table of the set starts at 3: {3} (6 > 5),
    // so nu(C) = nu(P) = 3 and C is contracted into P.
    const auto h = make(W{3}, W{5});
    const auto phys = phys_of(W{3}, W{5});
    REQUIRE(sssp::paper::validate_ordered_hierarchy(h, phys).ok());
    const auto nh = sssp::paper::normalize(h);
    REQUIRE(nh.forest.nodes.size() == 4);
    REQUIRE(nh.forest.roots.size() == 1);
    const NodeId r = nh.forest.roots[0];
    REQUIRE(nh.forest.nodes[r].children.size() == 3);
    REQUIRE(nh.forest.nodes[r].links.size() == 2);
    REQUIRE(nh.unit[r] == W{3});
    REQUIRE(nh.mass[r] == W{16});  // h(P) + link 5 + S(C) = 5 + 5 + 6 (topmost original mass)
    REQUIRE(nh.original_height[r] == W{5});
    REQUIRE(nh.node_of_input[c] == r);
    REQUIRE(nh.node_of_input[p] == r);
    check_normalized(h, phys, nh);
  }
  SECTION("linked with different units: dyadically aligned, not contracted") {
    // hc = 3, hp = 7: linked (24 >= 7); table {3, 6}: nu(C) = 3, nu(P) = 6.
    const auto h = make(W{3}, W{7});
    const auto phys = phys_of(W{3}, W{7});
    REQUIRE(sssp::paper::validate_ordered_hierarchy(h, phys).ok());
    const auto nh = sssp::paper::normalize(h);
    REQUIRE(nh.forest.nodes.size() == 5);
    REQUIRE(nh.unit[nh.node_of_input[c]] == W{3});
    REQUIRE(nh.unit[nh.node_of_input[p]] == W{6});
    check_normalized(h, phys, nh);
  }
  SECTION("4S(C) < h(P): unlinked, each node keeps nu = h and 2S(C) < nu(P)") {
    // hc = 1, hp = 9: 4 S(C) = 8 < 9, so no link; both are one-node sets (Q-24).
    const auto h = make(W{1}, W{9});
    const auto phys = phys_of(W{1}, W{9});
    REQUIRE(sssp::paper::validate_ordered_hierarchy(h, phys).ok());
    const auto nh = sssp::paper::normalize(h);
    REQUIRE(nh.unit[nh.node_of_input[c]] == W{1});
    REQUIRE(nh.unit[nh.node_of_input[p]] == W{9});
    check_normalized(h, phys, nh);
  }
}

TEMPLATE_TEST_CASE("normalize: light arcs are assigned to their LCA with length below its unit",
                   "[paper][hierarchy][normalize]", std::uint64_t, double) {
  using W = TestType;
  // C = {0,1} (h = 2) and leaf 2 under P (h = 8, nu = 8 since 4 S(C) = 16 >= 8 links them:
  // table {2, 4, 8}). Arc 0 -> 2 of length 5 < nu(P) = 8 is light at P; arc 2 -> 1 of length 8
  // is not; the self-arc 1 -> 1 is discarded; arc 0 -> 1 of length 1 < nu(C) = 2 is light at C.
  StateGraph<W> g;
  g.num_states = 3;
  g.arcs = {{0, 1, W{1}}, {1, 0, W{2}}, {0, 2, W{5}}, {2, 1, W{8}}, {1, 1, W{1}}};
  auto h = sssp::paper::make_leaf_forest<W>(g, {0, 1, 2});
  const NodeId c = sssp::paper::add_internal_node<W>(h, W{2}, {0, 1}, {{0, 1, W{2}}});
  const NodeId p = sssp::paper::add_internal_node<W>(h, W{8}, {c, 2}, {{0, 2, W{5}}});
  sssp::paper::collect_roots(h);
  const auto phys = sssp::paper::physical_graph_of_arcs<W>(g, {0, 1, 2}, 3);
  REQUIRE(sssp::paper::validate_ordered_hierarchy(h, phys).ok());
  const auto nh = sssp::paper::normalize(h);
  const NodeId nc = nh.node_of_input[c];
  const NodeId np = nh.node_of_input[p];
  REQUIRE(nh.unit[nc] == W{2});
  REQUIRE(nh.unit[np] == W{8});
  REQUIRE(nh.light(np).size() == 1);
  REQUIRE(nh.light(np)[0].arc == ArcId{2});
  REQUIRE(nh.light(np)[0].tail_child == nc);
  REQUIRE(nh.light(np)[0].head_child == nh.node_of_input[2]);
  REQUIRE(nh.light(nc).size() == 1);
  REQUIRE(nh.light(nc)[0].arc == ArcId{0});
  REQUIRE(nh.root_arcs.empty());
  check_normalized(h, phys, nh);
}

TEMPLATE_TEST_CASE("normalize: Lemma 4.3 on threshold-SCC reference hierarchies",
                   "[paper][hierarchy][normalize][random]", std::uint64_t, double) {
  using W = TestType;
  std::size_t contracted_total = 0;
  double kappa = 0.0;
  for (std::uint64_t seed = 1; seed <= 150; ++seed) {
    const RefInput in{2 + seed % 30, 1 + (seed * 7) % 90, 1 + seed % 20};
    const auto inst = reference_instance<W>(seed, in);
    REQUIRE(sssp::paper::validate_ordered_hierarchy(inst.h, inst.phys).ok());
    const auto nh = sssp::paper::normalize(inst.h);
    INFO("seed " << seed);
    kappa = std::max(kappa, check_normalized(inst.h, inst.phys, nh));
    contracted_total += inst.h.nodes.size() - nh.forest.nodes.size();
  }
  WARN("Lemma 4.3 reference hierarchies: max sum S/nu / N = " << kappa);
  // The generator family must actually exercise contraction (otherwise the test is vacuous).
  REQUIRE(contracted_total > 0);
}

TEMPLATE_TEST_CASE("normalize: Lemma 4.3 on Alg 4 base-case hierarchies",
                   "[paper][hierarchy][normalize][random]", std::uint64_t, double) {
  using W = TestType;
  double kappa = 0.0;
  for (int family = 0; family < 3; ++family) {
    for (std::uint64_t seed = 1; seed <= 40; ++seed) {
      const auto inst = base_case_instance<W>(seed, family);
      const auto nh = sssp::paper::normalize(inst.h);
      INFO("family " << family << " seed " << seed);
      kappa = std::max(kappa, check_normalized(inst.h, inst.phys, nh));
      // The base case has no light arcs (Lemma 3.5); normalisation cannot create any.
      for (NodeId x = 0; x < nh.forest.nodes.size(); ++x) REQUIRE(nh.light(x).empty());
    }
  }
  WARN("Lemma 4.3 base-case hierarchies: max sum S/nu / N = " << kappa);
}

TEST_CASE("normalize: inter-root arcs go to the root quotient in topological order",
          "[paper][hierarchy][normalize]") {
  using W = std::uint64_t;
  // Three single-state roots with 2 -> 0 and 0 -> 1: the only topological order is 2, 0, 1.
  StateGraph<W> g;
  g.num_states = 3;
  g.arcs = {{2, 0, 4}, {0, 1, 1}, {1, 1, 3}};
  auto h = sssp::paper::make_leaf_forest<W>(g, {0, 1, 2});
  sssp::paper::collect_roots(h);
  const auto nh = sssp::paper::normalize(h);
  REQUIRE(nh.root_arcs.size() == 2);
  REQUIRE(nh.root_order == std::vector<NodeId>{2, 0, 1});
  REQUIRE(nh.forest.c0 == 2);
}

TEST_CASE("normalize: a cyclic root quotient violates H4 and is rejected",
          "[paper][hierarchy][normalize]") {
  using W = std::uint64_t;
  StateGraph<W> g;
  g.num_states = 2;
  g.arcs = {{0, 1, 1}, {1, 0, 1}};
  auto h = sssp::paper::make_leaf_forest<W>(g, {0, 1});
  sssp::paper::collect_roots(h);
  REQUIRE_THROWS_AS(sssp::paper::normalize(h), std::invalid_argument);
}

TEST_CASE("normalize records path counters in Diagnostics", "[paper][hierarchy][normalize]") {
  using W = std::uint64_t;
  const auto inst = reference_instance<W>(7, {20, 60, 10});
  sssp::Diagnostics diag;
  const auto nh = sssp::paper::normalize(inst.h, &diag);
#if SSSP_DIAGNOSTICS
  REQUIRE(diag.counters["normalize_calls"] == 1);
  REQUIRE(diag.counters["normalize_contracted_nodes"] ==
          inst.h.nodes.size() - nh.forest.nodes.size());
  std::size_t light = 0;
  for (NodeId x = 0; x < nh.forest.nodes.size(); ++x) light += nh.light(x).size();
  REQUIRE(diag.counters["normalize_light_arcs"] == light);
  REQUIRE(diag.counters["normalize_root_arcs"] == nh.root_arcs.size());
  REQUIRE(diag.ops.additions > 0);
  REQUIRE(diag.ops.comparisons > 0);
#else
  static_cast<void>(nh);
#endif
}

TEST_CASE("normalize: a mass beyond the uint64 range saturates at infinity (D-13 / D-23)",
          "[paper][hierarchy][normalize]") {
  using W = std::uint64_t;
  // The hand example above with hc = 2^63 and hp = 3 * 2^62: S(C) = 2 hc = 2^64 and S(P) =
  // hp + hp + S(C) do not fit. Both are upper bounds (Sec 4.1 l.686-688), so the saturated value
  // +infinity is still one; 4 S(C) >= h(P) holds for the true value and for infinity, so C is
  // linked to P. The table of the set starts at 2^63 and stops there (2^64 does not fit), so
  // nu(C) = nu(P) = 2^63 and C is contracted into P, which keeps its own (saturated) mass.
  constexpr W hc = W{1} << 63;
  constexpr W hp = W{3} << 62;
  StateGraph<W> g;
  g.num_states = 3;
  g.arcs = {{0, 1, hc}, {1, 0, hc}, {1, 2, hp}, {2, 1, hp}};
  auto h = sssp::paper::make_leaf_forest<W>(g, {0, 1, 2});
  const NodeId c = sssp::paper::add_internal_node<W>(h, hc, {0, 1}, {{0, 1, hc}});
  const NodeId p = sssp::paper::add_internal_node<W>(h, hp, {c, 2}, {{1, 2, hp}});
  sssp::paper::collect_roots(h);
  const auto phys = CsrGraph<W>::from_edges(3, {{0, 1, hc}, {1, 2, hp}});
  REQUIRE(sssp::paper::validate_ordered_hierarchy(h, phys).ok());
  NormalizedHierarchy<W> nh;
  REQUIRE_NOTHROW(nh = sssp::paper::normalize(h));
  REQUIRE(nh.forest.nodes.size() == 4);  // three leaves and P (C contracted)
  const NodeId x = nh.node_of_input[p];
  CHECK(nh.node_of_input[c] == x);
  CHECK(nh.unit[x] == hc);
  CHECK(nh.mass[x] == kInf<W>);
  CHECK(nh.original_height[x] == hp);
}
