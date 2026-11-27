// Tests for the child order of Sec 4.2 (l.725-779, Lemma 4.4, Lemma 4.5).
//
// Expectations come from the paper's text: hand examples are traced by hand through the DFS of
// Sec 4.2 (l.726-730) and the charging of the proof of Lemma 4.4 (l.737-751); random inputs are
// checked against independent oracles (Floyd-Warshall physical distances, a direct traversal of
// the link forest F, the bounds of Lemma 4.4 and 4.5 recomputed here).
//
// Inputs: normalize() of the threshold-SCC reference hierarchy on random directed state graphs
// and of the Alg 4 base case on generated graphs; several seeds, std::uint64_t and double.
//
// Tolerance (double only): gaps are sums of forest distances rd(a) + rd(b) - 2 rd(lca) and arc
// lengths. Each root distance is a rounded sum of at most N terms, each bounded by
// `scale` = (total link length + total arc length), so it errs by <= N eps scale; a forest
// distance by <= 4 N eps scale; a gap (at most 2N such terms) by <= 8 N^2 eps scale, which is
// < 1e-11 scale for N <= 64 (eps = 2^-53). Comparisons with Floyd-Warshall distances and with
// the Lemma 4.4 bound therefore allow an absolute slack of 1e-10 * scale. Integer weights are
// compared exactly.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "paper/local/base_case.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Vertex;
using sssp::paper::kNoNode;
using sssp::paper::kNoState;
using sssp::paper::NodeId;
using sssp::paper::NormalizedHierarchy;
using sssp::paper::OrderedHierarchy;
using sssp::paper::StateGraph;
using sssp::paper::StateId;

namespace {

template <class W>
constexpr W kInf = sssp::WeightTraits<W>::infinity();

// ------------------------------------------------------------------ input generators

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

template <class W>
struct Instance {
  OrderedHierarchy<W> h;
  CsrGraph<W> phys;
};

template <class W>
Instance<W> reference_instance(std::uint64_t seed, std::size_t n, std::size_t m, std::size_t p) {
  std::mt19937_64 rng(seed);
  auto g = random_state_graph<W>(rng, n, m);
  std::vector<Vertex> pi(n);
  for (auto& x : pi) x = static_cast<Vertex>(rng() % p);
  auto phys = sssp::paper::physical_graph_of_arcs<W>(g, pi, p);
  auto h = sssp::paper::build_threshold_scc_hierarchy<W>(std::move(g), std::move(pi));
  return {std::move(h), std::move(phys)};
}

template <class W>
Instance<W> base_case_instance(std::uint64_t seed, int family) {
  sssp::gen::GenOptions opt;
  opt.seed = seed;
  opt.connected = (seed % 3) != 0;
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

// Distances in the link forest F (Sec 4.1 l.676-679) by a direct search over the links:
// fdist[a][b] = length of the unique F path, kInf if a and b lie in different trees.
template <class W>
std::vector<std::vector<W>> forest_distances(const OrderedHierarchy<W>& f) {
  const std::size_t n = f.num_states();
  std::vector<std::vector<std::pair<StateId, W>>> adj(n);
  for (const auto& nd : f.nodes) {
    for (const auto& lk : nd.links) {
      adj[lk.x].push_back({lk.y, lk.length});
      adj[lk.y].push_back({lk.x, lk.length});
    }
  }
  std::vector<std::vector<W>> d(n, std::vector<W>(n, kInf<W>));
  for (StateId s = 0; s < n; ++s) {
    std::vector<StateId> stack{s};
    d[s][s] = W{0};
    while (!stack.empty()) {
      const StateId u = stack.back();
      stack.pop_back();
      for (const auto& [v, w] : adj[u]) {
        if (d[s][v] == kInf<W>) {
          d[s][v] = d[s][u] + w;
          stack.push_back(v);
        }
      }
    }
  }
  return d;
}

template <class W>
void states_below(const OrderedHierarchy<W>& h, NodeId v, std::vector<StateId>& out) {
  if (h.nodes[v].children.empty()) {
    out.push_back(h.nodes[v].state);
    return;
  }
  for (const NodeId c : h.nodes[v].children) states_below(h, c, out);
}

template <class W>
W total_scale(const OrderedHierarchy<W>& f) {
  W s{0};
  for (const auto& nd : f.nodes)
    for (const auto& lk : nd.links) s += lk.length;
  for (const auto& a : f.graph.arcs) s += a.length;
  return s;
}

template <class W>
bool le_tol(W a, W b, W scale) {
  if constexpr (std::is_integral_v<W>) {
    static_cast<void>(scale);
    return a <= b;
  } else {
    return a <= b + 1e-10 * scale;
  }
}

struct Ratios {
  double gap_over_unit = 0.0;  ///< sum_X sum_k g / nu(X), divided by N
};

// Checks every claim of Sec 4.2 / Lemma 4.4 / Lemma 4.5 on an ordered NormalizedHierarchy.
// `before` is the same hierarchy before order_children (child sets, light arcs).
template <class W>
Ratios check_ordered(const NormalizedHierarchy<W>& before, const NormalizedHierarchy<W>& nh,
                     const CsrGraph<W>& phys) {
  const auto& f = nh.forest;
  const std::size_t nf = f.nodes.size();
  REQUIRE(nh.children_ordered);
  REQUIRE(nh.rank.size() == nf);
  REQUIRE(nh.rep.size() == nf);
  REQUIRE(nh.gap_next.size() == nf);

  const auto lf = sssp::paper::build_link_forest(nh);
  const auto dist = floyd_warshall(phys);
  const W scale = total_scale(f);

  // The link forest: one edge per link, parents precede children in the preorder, and the root
  // distances add up the parent links.
  std::size_t links = 0;
  for (const auto& nd : f.nodes) links += nd.links.size();
  std::size_t f_roots = 0;
  for (StateId s = 0; s < f.num_states(); ++s) {
    REQUIRE(lf.order[lf.preorder[s]] == s);
    if (lf.parent[s] == kNoState) {
      ++f_roots;
      REQUIRE(lf.root_distance[s] == W{0});
    } else {
      REQUIRE(lf.preorder[lf.parent[s]] < lf.preorder[s]);
      REQUIRE(lf.root_distance[s] == lf.root_distance[lf.parent[s]] + lf.parent_length[s]);
    }
  }
  REQUIRE(f_roots + links == f.num_states());
  const auto fdist = forest_distances(f);

  // Roots: rank = position in the root order.
  for (std::size_t i = 0; i < nh.root_order.size(); ++i) {
    REQUIRE(nh.rank[nh.root_order[i]] == i);
    REQUIRE(nh.rep[nh.root_order[i]] == kNoState);
  }

  Ratios out;
  for (NodeId x = 0; x < nf; ++x) {
    const auto& ch = f.nodes[x].children;
    if (ch.empty()) continue;
    // Same children as before, now in rank order.
    auto a = ch;
    auto b = before.forest.nodes[x].children;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    REQUIRE(a == b);
    for (std::size_t k = 0; k < ch.size(); ++k) {
      REQUIRE(nh.rank[ch[k]] == k);
      REQUIRE(f.nodes[ch[k]].parent == x);
    }

    // LightArcsForward (Lemma 4.5: a topological order of the light child DAG).
    for (const auto& la : nh.light(x)) REQUIRE(nh.rank[la.tail_child] < nh.rank[la.head_child]);

    // Representatives: a state below the child; either its F-preorder-minimum state (DFS root,
    // provisional representative, Q-25) or the head of a light arc from a lower-ranked child
    // (DFS tree arc, l.729).
    W mass_children{0};
    for (const NodeId c : ch) {
      std::vector<StateId> st;
      states_below(f, c, st);
      REQUIRE(std::find(st.begin(), st.end(), nh.rep[c]) != st.end());
      const StateId prov = *std::min_element(st.begin(), st.end(), [&](StateId p, StateId q) {
        return lf.preorder[p] < lf.preorder[q];
      });
      bool entered = false;
      for (const auto& la : nh.light(x)) {
        entered = entered || (la.head_child == c && f.graph.arcs[la.arc].head == nh.rep[c] &&
                              nh.rank[la.tail_child] < nh.rank[c]);
      }
      REQUIRE((nh.rep[c] == prov || entered));
      mass_children += nh.mass[c];
    }

    // GapIsPhysicalUpperBound and TourBoundLemma44.
    W gaps{0};
    double gap_over_unit = 0.0;
    for (std::size_t k = 0; k + 1 < ch.size(); ++k) {
      const W g = nh.gap_next[ch[k]];
      const W d = dist[f.projection[nh.rep[ch[k]]]][f.projection[nh.rep[ch[k + 1]]]];
      REQUIRE(d != kInf<W>);
      REQUIRE(le_tol(d, g, scale));
      gaps += g;
      gap_over_unit += static_cast<double>(g) / static_cast<double>(nh.unit[x]);
    }
    REQUIRE(nh.gap_next[ch.back()] == W{0});
    W d_nu{0};
    for (std::size_t k = 0; k < ch.size(); ++k) d_nu += nh.unit[x];
    const W bound = mass_children + mass_children + d_nu + d_nu + nh.mass[x] + nh.mass[x];
    REQUIRE(le_tol(gaps, bound, scale));
    // Lemma 4.5 proof: the child masses sum to at most S(X).
    REQUIRE(mass_children <= nh.mass[x]);
    // The F path between consecutive representatives never exceeds the Lemma 4.4 total either;
    // in particular F is connected below X.
    for (std::size_t k = 0; k + 1 < ch.size(); ++k) {
      REQUIRE(fdist[nh.rep[ch[k]]][nh.rep[ch[k + 1]]] != kInf<W>);
    }
    out.gap_over_unit += gap_over_unit;
  }
  out.gap_over_unit /= static_cast<double>(f.num_states());
  return out;
}

// GapOverUnitLinear (Lemma 4.5): sum_X sum g/nu(X) <= sum_X (2 sum_i S(i) + 2 d nu + 2 S(X)) /
// nu(X) <= 4 sum_X S(X)/nu(X) + 2 sum_X d(X), using sum_i S(i) <= S(X).
template <class W>
double lemma45_bound(const NormalizedHierarchy<W>& nh) {
  double sum = 0.0;
  for (NodeId x = 0; x < nh.forest.nodes.size(); ++x) {
    const auto& ch = nh.forest.nodes[x].children;
    if (ch.empty()) continue;
    sum += 4.0 * static_cast<double>(nh.mass[x]) / static_cast<double>(nh.unit[x]) +
           2.0 * static_cast<double>(ch.size());
  }
  return sum / static_cast<double>(nh.forest.num_states());
}

}  // namespace

TEMPLATE_TEST_CASE("child order: hand example with a cyclic start inside a child",
                   "[paper][hierarchy][child_order]", std::uint64_t, double) {
  using W = TestType;
  // States 0..5 on the link path 0 - 1 - 2 - 3 - 4 - 5 (F rooted at state 0, so the F preorder
  // is 0, 1, 2, 3, 4, 5). A = {1, 2, 3} has height 1 and links (1,2,1), (2,3,1); X has children
  // Z = leaf 0, A, leaf 4, leaf 5, height 8 and links (0,1,2), (3,4,2), (4,5,2). Light arcs at X
  // (all < nu(X) = 8; S(A) = 3 links A to X, table {1, 2, 4, 8}): 0 -> 2 (4), 3 -> 4 (3),
  // 1 -> 5 (5).
  //
  // Hand trace (Sec 4.2): DFS roots in F preorder of the provisional representatives: Z first.
  // Z: arc 0 -> 2 enters A at state 2. A: outgoing tails in cyclic preorder from 2: 3, then 1.
  //   3 -> 4: move F(2,3) = 1, arc 3; leaf 4 finishes [mark 1]; return 3.
  //   1 -> 5: move F(3,1) = 2, arc 5; leaf 5 finishes [mark 2]; g(4,5) = 3 + 2 + 5 = 10.
  //   return 5; move F(1,2) = 1; A finishes [mark 3]; g(5,A) = 5 + 1 = 6.
  // back in Z: return 4; Z finishes [mark 4]; g(A,Z) = 4.
  // Rank order = reverse finish = Z, A, 5, 4; representatives 0, 2, 5, 4.
  StateGraph<W> g;
  g.num_states = 6;
  g.arcs = {{0, 2, W{4}}, {3, 4, W{3}}, {1, 5, W{5}}};
  auto h = sssp::paper::make_leaf_forest<W>(g, {0, 1, 2, 3, 4, 5});
  const NodeId a =
      sssp::paper::add_internal_node<W>(h, W{1}, {1, 2, 3}, {{1, 2, W{1}}, {2, 3, W{1}}});
  sssp::paper::add_internal_node<W>(h, W{8}, {0, a, 4, 5},
                                    {{0, 1, W{2}}, {3, 4, W{2}}, {4, 5, W{2}}});
  sssp::paper::collect_roots(h);
  const auto phys = CsrGraph<W>::from_edges(6, {{0, 1, W{2}},
                                                {1, 2, W{1}},
                                                {2, 3, W{1}},
                                                {3, 4, W{2}},
                                                {4, 5, W{2}},
                                                {0, 2, W{4}},
                                                {3, 4, W{3}},
                                                {1, 5, W{5}}});
  REQUIRE(sssp::paper::validate_ordered_hierarchy(h, phys).ok());

  auto nh = sssp::paper::normalize(h);
  const auto before = nh;
  sssp::paper::order_children(nh);
  const NodeId x = nh.forest.roots[0];
  const NodeId na = nh.node_of_input[a];
  const NodeId z = nh.node_of_input[0];
  const NodeId l4 = nh.node_of_input[4];
  const NodeId l5 = nh.node_of_input[5];
  REQUIRE(nh.unit[x] == W{8});
  REQUIRE(nh.forest.nodes[x].children == std::vector<NodeId>{z, na, l5, l4});
  REQUIRE(nh.rep[z] == StateId{0});
  REQUIRE(nh.rep[na] == StateId{2});
  REQUIRE(nh.rep[l5] == StateId{5});
  REQUIRE(nh.rep[l4] == StateId{4});
  REQUIRE(nh.gap_next[z] == W{4});
  REQUIRE(nh.gap_next[na] == W{6});
  REQUIRE(nh.gap_next[l5] == W{10});
  REQUIRE(nh.gap_next[l4] == W{0});
  // Inside A there are no light arcs: its leaves keep the F preorder reversed (finish order of
  // three DFS roots 1, 2, 3 is 1, 2, 3), gaps are F distances 1 and 1.
  REQUIRE(nh.forest.nodes[na].children ==
          std::vector<NodeId>{nh.node_of_input[3], nh.node_of_input[2], nh.node_of_input[1]});
  REQUIRE(nh.gap_next[nh.node_of_input[3]] == W{1});
  REQUIRE(nh.gap_next[nh.node_of_input[2]] == W{1});
  check_ordered(before, nh, phys);
}

TEMPLATE_TEST_CASE("child order: children without outgoing light arcs finish immediately",
                   "[paper][hierarchy][child_order]", std::uint64_t, double) {
  using W = TestType;
  // Base-case hierarchies have no light arcs (Lemma 3.5): every child is a DFS root entered at
  // its provisional representative (the F-preorder minimum below it, Q-25) and finishes at
  // once. Hence the rank order is the reverse F preorder of those representatives and every gap
  // is exactly the F distance between consecutive representatives (l.747: moves between
  // successive DFS roots along F).
  for (int family = 0; family < 3; ++family) {
    for (std::uint64_t seed = 1; seed <= 20; ++seed) {
      const auto inst = base_case_instance<W>(seed, family);
      auto nh = sssp::paper::normalize(inst.h);
      const auto before = nh;
      sssp::paper::order_children(nh);
      const auto lf = sssp::paper::build_link_forest(nh);
      const auto fdist = forest_distances(nh.forest);
      INFO("family " << family << " seed " << seed);
      for (NodeId x = 0; x < nh.forest.nodes.size(); ++x) {
        const auto& ch = nh.forest.nodes[x].children;
        if (ch.empty()) continue;
        REQUIRE(nh.light(x).empty());
        for (const NodeId c : ch) {
          std::vector<StateId> st;
          states_below(nh.forest, c, st);
          const StateId prov = *std::min_element(st.begin(), st.end(), [&](StateId p, StateId q) {
            return lf.preorder[p] < lf.preorder[q];
          });
          REQUIRE(nh.rep[c] == prov);
        }
        for (std::size_t k = 0; k + 1 < ch.size(); ++k) {
          REQUIRE(lf.preorder[nh.rep[ch[k]]] > lf.preorder[nh.rep[ch[k + 1]]]);
          const W fd = fdist[nh.rep[ch[k]]][nh.rep[ch[k + 1]]];
          if constexpr (std::is_integral_v<W>) {
            REQUIRE(nh.gap_next[ch[k]] == fd);
          } else {
            REQUIRE(le_tol(fd, nh.gap_next[ch[k]], total_scale(nh.forest)));
            REQUIRE(le_tol(nh.gap_next[ch[k]], fd, total_scale(nh.forest)));
          }
        }
      }
      check_ordered(before, nh, inst.phys);
    }
  }
}

TEMPLATE_TEST_CASE("child order: Lemma 4.4 and 4.5 on threshold-SCC reference hierarchies",
                   "[paper][hierarchy][child_order][random]", std::uint64_t, double) {
  using W = TestType;
  double kappa = 0.0;
  std::size_t light_total = 0;
  for (std::uint64_t seed = 1; seed <= 150; ++seed) {
    const auto inst =
        reference_instance<W>(seed, 2 + seed % 30, 1 + (seed * 7) % 90, 1 + seed % 20);
    auto nh = sssp::paper::normalize(inst.h);
    const auto before = nh;
    sssp::paper::order_children(nh);
    INFO("seed " << seed);
    const auto r = check_ordered(before, nh, inst.phys);
    REQUIRE(r.gap_over_unit <= lemma45_bound(nh) + 1e-9);  // GapOverUnitLinear
    kappa = std::max(kappa, r.gap_over_unit);
    light_total += nh.light_arcs.size();
  }
  REQUIRE(light_total > 0);  // the family must exercise light arcs
  WARN("Lemma 4.5 reference hierarchies: max sum g/nu / N = " << kappa);
}

TEMPLATE_TEST_CASE("child order: Lemma 4.4 and 4.5 on Alg 4 base-case hierarchies",
                   "[paper][hierarchy][child_order][random]", std::uint64_t, double) {
  using W = TestType;
  double kappa = 0.0;
  for (int family = 0; family < 3; ++family) {
    for (std::uint64_t seed = 21; seed <= 50; ++seed) {
      const auto inst = base_case_instance<W>(seed, family);
      auto nh = sssp::paper::normalize(inst.h);
      const auto before = nh;
      sssp::paper::order_children(nh);
      INFO("family " << family << " seed " << seed);
      const auto r = check_ordered(before, nh, inst.phys);
      REQUIRE(r.gap_over_unit <= lemma45_bound(nh) + 1e-9);
      kappa = std::max(kappa, r.gap_over_unit);
    }
  }
  WARN("Lemma 4.5 base-case hierarchies: max sum g/nu / N = " << kappa);
}

TEST_CASE("child order: a cyclic light quotient violates H3 and is rejected",
          "[paper][hierarchy][child_order]") {
  using W = std::uint64_t;
  // Two leaves under X (h = 4) with light arcs both ways: Def 4.1 H3 fails at X.
  StateGraph<W> g;
  g.num_states = 2;
  g.arcs = {{0, 1, 1}, {1, 0, 1}};
  auto h = sssp::paper::make_leaf_forest<W>(g, {0, 1});
  sssp::paper::add_internal_node<W>(h, 4, {0, 1}, {{0, 1, 4}});
  sssp::paper::collect_roots(h);
  auto nh = sssp::paper::normalize(h);
  REQUIRE(nh.light_arcs.size() == 2);
  REQUIRE_THROWS_AS(sssp::paper::order_children(nh), std::invalid_argument);
}

TEST_CASE("child order records path counters in Diagnostics", "[paper][hierarchy][child_order]") {
  using W = std::uint64_t;
  const auto inst = reference_instance<W>(11, 24, 70, 12);
  auto nh = sssp::paper::normalize(inst.h);
  sssp::Diagnostics diag;
  sssp::paper::order_children(nh, &diag);
#if SSSP_DIAGNOSTICS
  std::size_t internal = 0;
  for (const auto& nd : nh.forest.nodes) internal += nd.children.empty() ? 0 : 1;
  REQUIRE(diag.counters["child_order_calls"] == 1);
  REQUIRE(diag.counters["child_order_nodes"] == internal);
  REQUIRE(diag.counters["child_order_light_arcs"] == nh.light_arcs.size());
  REQUIRE(diag.counters["child_order_dfs_roots"] >= internal);
  REQUIRE(diag.ops.additions > 0);
#endif
}
