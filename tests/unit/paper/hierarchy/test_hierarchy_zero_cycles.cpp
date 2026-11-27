// Zero-length directed cycles through two or more states (open question Q-21) against the
// Def 4.1 checker (Sec 4, l.646-658).
//
// Expected behaviour, derived from Def 4.1 and not from the checker:
//  * Def 4.1 allows nonnegative arc lengths (l.646) but requires a positive height at every
//    internal node (l.649). An arc of length 0 between distinct children of C is therefore
//    always light (length strictly below h(C), l.650).
//  * Let Z be a zero-length cycle through >= 2 distinct states.
//    - If the states of Z lie in one tree, let C be the lowest node containing all of them. C is
//      internal and Z meets >= 2 of its children; every arc of Z between distinct children is
//      light, so the light quotient at C contains a closed walk through >= 2 children, i.e. a
//      cycle: H3 is violated at C (l.650).
//    - Otherwise Z meets >= 2 roots and the root quotient on all inter-root arcs has a cycle:
//      H4 is violated (l.651).
//    Hence no ordered hierarchy exists for such a state graph, whatever the forest and heights;
//    the caller must contract zero-length components first (Sec 5.1 l.961-962; Q-21 decision).
//  * Zero-length arcs that close no cycle, and zero-length self-arcs (discarded by
//    normalisation, Sec 4.1 l.719-720), are allowed.
// Lowest common nodes, H1 walks and the physical graph are computed here independently.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <type_traits>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "sssp/graph.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
using sssp::paper::Condition;
using sssp::paper::kNoArc;
using sssp::paper::kNoNode;
using sssp::paper::NodeId;
using sssp::paper::OrderedHierarchy;
using sssp::paper::StateArc;
using sssp::paper::StateGraph;
using sssp::paper::StateId;
using sssp::paper::ValidationReport;
using sssp::paper::ViolationKind;

namespace {

// Every violation is an H3 or H4 violation (nothing else is wrong with the inputs below).
bool only_h3_h4(const ValidationReport& rep) {
  return std::all_of(rep.violations.begin(), rep.violations.end(), [](const auto& v) {
    return v.condition() == Condition::kH3 || v.condition() == Condition::kH4;
  });
}

bool has_light_cycle_at(const ValidationReport& rep, NodeId node) {
  return std::any_of(rep.violations.begin(), rep.violations.end(), [&](const auto& v) {
    return v.kind == ViolationKind::kLightArcCycle && v.node == node;
  });
}

// Lowest node whose subtree contains every state in `states` (kNoNode if they span >= 2
// trees). Naive parent walks, independent of the checker.
template <class W>
NodeId lowest_common_node(const OrderedHierarchy<W>& h, const std::vector<StateId>& states) {
  std::vector<NodeId> common;
  for (NodeId v = h.leaf_of_state[states[0]]; v != kNoNode; v = h.nodes[v].parent) {
    common.push_back(v);
  }
  for (std::size_t i = 1; i < states.size(); ++i) {
    std::vector<NodeId> up;
    for (NodeId v = h.leaf_of_state[states[i]]; v != kNoNode; v = h.nodes[v].parent) {
      up.push_back(v);
    }
    std::erase_if(common, [&](NodeId c) { return std::find(up.begin(), up.end(), c) == up.end(); });
  }
  return common.empty() ? kNoNode : common.front();  // ancestors are listed bottom-up
}

// One physical edge {pi(x), pi(y)} of length l per state arc: every arc is its own walk (H1).
template <class W>
CsrGraph<W> own_physical_graph(const StateGraph<W>& g, const std::vector<Vertex>& pi,
                               std::size_t n_phys) {
  std::vector<Edge<W>> edges;
  for (const auto& a : g.arcs) edges.push_back({pi[a.tail], pi[a.head], a.length});
  return CsrGraph<W>::from_edges(n_phys, edges);
}

// All states map to the single physical vertex 0, so every arc and link is a walk of length 0.
template <class W>
CsrGraph<W> one_point() {
  return CsrGraph<W>::from_edges(1, std::vector<Edge<W>>{});
}

template <class W>
OrderedHierarchy<W> leaves_at_one_point(std::size_t n, std::vector<StateArc<W>> arcs) {
  StateGraph<W> d;
  d.num_states = n;
  d.arcs = std::move(arcs);
  return sssp::paper::make_leaf_forest<W>(std::move(d), std::vector<Vertex>(n, 0));
}

template <class W>
std::vector<W> positive_heights() {
  if constexpr (std::is_floating_point_v<W>) {
    return {std::numeric_limits<W>::denorm_min(), W{1}, W{1e300}};
  } else {
    return {W{1}, W{1000}, std::numeric_limits<W>::max() / 2};
  }
}

}  // namespace

TEMPLATE_TEST_CASE("Q-21: a zero-length 2-cycle below one node violates H3 there, at any height",
                   "[paper][hierarchy][def41][h3][q21]", std::uint64_t, double) {
  using W = TestType;
  for (const W height : positive_heights<W>()) {
    INFO("height " << height);
    auto h = leaves_at_one_point<W>(2, {{0, 1, W{0}}, {1, 0, W{0}}});
    const NodeId c = sssp::paper::add_internal_node<W>(h, height, {0, 1}, {{0, 1, W{0}, kNoArc}});
    sssp::paper::collect_roots(h);
    const auto rep = sssp::paper::validate_ordered_hierarchy(h, one_point<W>());
    INFO(rep.summary());
    REQUIRE_FALSE(rep.ok());
    REQUIRE(rep.violations.size() == 1);
    REQUIRE(rep.violations[0].kind == ViolationKind::kLightArcCycle);
    REQUIRE(rep.violations[0].node == c);
  }
}

TEMPLATE_TEST_CASE("Q-21: a zero-length 2-cycle between two roots violates H4",
                   "[paper][hierarchy][def41][h4][q21]", std::uint64_t, double) {
  using W = TestType;
  auto h = leaves_at_one_point<W>(2, {{0, 1, W{0}}, {1, 0, W{0}}});
  sssp::paper::collect_roots(h);
  const auto rep = sssp::paper::validate_ordered_hierarchy(h, one_point<W>());
  INFO(rep.summary());
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kRootArcCycle);
}

TEMPLATE_TEST_CASE("Q-21: a zero-length 3-cycle is reported at the lowest node containing it",
                   "[paper][hierarchy][def41][h3][q21]", std::uint64_t, double) {
  using W = TestType;
  // 0 -> 1 -> 2 -> 0, all of length 0. A = {0, 1} (h = 1) sees only 0 -> 1 between its
  // children (acyclic); R = {A, 2} (h = 2) sees A -> 2 and 2 -> A: a light cycle at R only.
  auto h = leaves_at_one_point<W>(3, {{0, 1, W{0}}, {1, 2, W{0}}, {2, 0, W{0}}});
  const NodeId a = sssp::paper::add_internal_node<W>(h, W{1}, {0, 1}, {{0, 1, W{0}, kNoArc}});
  const NodeId r = sssp::paper::add_internal_node<W>(h, W{2}, {a, 2}, {{1, 2, W{0}, kNoArc}});
  sssp::paper::collect_roots(h);
  REQUIRE(lowest_common_node(h, {0, 1, 2}) == r);
  const auto rep = sssp::paper::validate_ordered_hierarchy(h, one_point<W>());
  INFO(rep.summary());
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kLightArcCycle);
  REQUIRE(rep.violations[0].node == r);
}

TEMPLATE_TEST_CASE("Q-21: a zero-length cycle inside a child is reported at that child only",
                   "[paper][hierarchy][def41][h3][q21]", std::uint64_t, double) {
  using W = TestType;
  // Physical edge {0, 1} of length 5. States 0, 1 at vertex 0, state 2 at vertex 1.
  // Zero cycle 0 <-> 1 inside A = {0, 1}; arcs 1 -> 2 and 2 -> 0 of length 5 are not light at
  // R = {A, 2} (h = 5, not strictly below), so only A is violated.
  StateGraph<W> d;
  d.num_states = 3;
  d.arcs = {{0, 1, W{0}}, {1, 0, W{0}}, {1, 2, W{5}}, {2, 0, W{5}}};
  auto h = sssp::paper::make_leaf_forest<W>(std::move(d), {0, 0, 1});
  const NodeId a = sssp::paper::add_internal_node<W>(h, W{1}, {0, 1}, {{0, 1, W{0}, kNoArc}});
  const NodeId r = sssp::paper::add_internal_node<W>(h, W{5}, {a, 2}, {{1, 2, W{5}, 2}});
  sssp::paper::collect_roots(h);
  const auto phys = CsrGraph<W>::from_edges(2, {{0, 1, W{5}}});
  REQUIRE(lowest_common_node(h, {0, 1}) == a);
  REQUIRE(lowest_common_node(h, {0, 1, 2}) == r);
  const auto rep = sssp::paper::validate_ordered_hierarchy(h, phys);
  INFO(rep.summary());
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kLightArcCycle);
  REQUIRE(rep.violations[0].node == a);
}

TEMPLATE_TEST_CASE("Q-21: a zero-length cycle across two trees violates H4, not H3",
                   "[paper][hierarchy][def41][h4][q21]", std::uint64_t, double) {
  using W = TestType;
  // A = {0, 1} is a root, state 2 is another root. The zero cycle 0 -> 1 -> 2 -> 0 meets both.
  auto h = leaves_at_one_point<W>(3, {{0, 1, W{0}}, {1, 2, W{0}}, {2, 0, W{0}}});
  sssp::paper::add_internal_node<W>(h, W{1}, {0, 1}, {{0, 1, W{0}, kNoArc}});
  sssp::paper::collect_roots(h);
  REQUIRE(lowest_common_node(h, {0, 1, 2}) == kNoNode);
  const auto rep = sssp::paper::validate_ordered_hierarchy(h, one_point<W>());
  INFO(rep.summary());
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kRootArcCycle);
}

TEMPLATE_TEST_CASE("Q-21: zero-length arcs on no cycle and zero-length self-arcs are accepted",
                   "[paper][hierarchy][def41][q21]", std::uint64_t, double) {
  using W = TestType;
  // 0 -> 1 has length 0 (light at C), 1 -> 0 has length 1 = h(C) (not light): acyclic.
  auto h = leaves_at_one_point<W>(
      2, {{0, 1, W{0}}, {1, 0, W{1}}, {0, 0, W{0}}, {1, 1, W{0}}, {1, 1, W{3}}});
  sssp::paper::add_internal_node<W>(h, W{1}, {0, 1}, {{0, 1, W{0}, 0}});
  sssp::paper::collect_roots(h);
  const auto rep = sssp::paper::validate_ordered_hierarchy(h, one_point<W>());
  INFO(rep.summary());
  REQUIRE(rep.ok());
  // Two roots with only a zero self-arc and a one-way zero arc: acyclic root quotient.
  auto roots = leaves_at_one_point<W>(2, {{0, 1, W{0}}, {0, 0, W{0}}});
  sssp::paper::collect_roots(roots);
  const auto rep2 = sssp::paper::validate_ordered_hierarchy(roots, one_point<W>());
  INFO(rep2.summary());
  REQUIRE(rep2.ok());
}

TEMPLATE_TEST_CASE("Q-21: no forest on three states admits a zero-length cycle through >= 2 states",
                   "[paper][hierarchy][def41][q21][exhaustive]", std::uint64_t, double) {
  using W = TestType;
  // Cycles through >= 2 distinct states of {0, 1, 2}.
  const std::vector<std::vector<StateId>> cycles = {{0, 1}, {0, 2}, {1, 2}, {0, 1, 2}, {0, 2, 1}};
  // Forest shapes: pair = the two states joined first (or none), top = whether a node joins
  // everything. Heights range over {1, 2, 3} with h(pair) <= h(top).
  struct Shape {
    int pair;  // -1: no pair node; otherwise the state left out of the pair
    bool top;  // a node over everything (the pair node and the remaining leaf, or all three)
  };
  const std::vector<Shape> shapes = {{-1, false}, {-1, true}, {0, false}, {1, false},
                                     {2, false},  {0, true},  {1, true},  {2, true}};
  for (const W len : {W{0}, W{5}}) {
    for (const auto& cyc : cycles) {
      std::vector<StateArc<W>> arcs;
      for (std::size_t i = 0; i < cyc.size(); ++i) {
        arcs.push_back({cyc[i], cyc[(i + 1) % cyc.size()], len});
      }
      std::size_t accepted = 0;
      for (const auto& shape : shapes) {
        for (std::uint64_t hp = 1; hp <= 3; ++hp) {
          for (std::uint64_t ht = hp; ht <= 3; ++ht) {
            auto h = leaves_at_one_point<W>(3, arcs);
            if (shape.pair < 0) {
              // Heights of absent nodes do not vary: all-roots runs once, one top node runs
              // once per height.
              if (!shape.top && (hp != 1 || ht != 1)) continue;
              if (shape.top && hp != ht) continue;
              if (shape.top) {
                sssp::paper::add_internal_node<W>(h, static_cast<W>(ht), {0, 1, 2},
                                                  {{0, 1, W{0}, kNoArc}, {1, 2, W{0}, kNoArc}});
              }
            } else {
              if (!shape.top && hp != ht) continue;
              const auto out = static_cast<StateId>(shape.pair);
              std::vector<NodeId> in;
              for (StateId x = 0; x < 3; ++x) {
                if (x != out) in.push_back(x);
              }
              const NodeId p = sssp::paper::add_internal_node<W>(h, static_cast<W>(hp), in,
                                                                 {{in[0], in[1], W{0}, kNoArc}});
              if (shape.top) {
                sssp::paper::add_internal_node<W>(h, static_cast<W>(ht), {p, out},
                                                  {{in[0], out, W{0}, kNoArc}});
              }
            }
            sssp::paper::collect_roots(h);
            const auto rep = sssp::paper::validate_ordered_hierarchy(h, one_point<W>());
            INFO("len " << len << " cycle size " << cyc.size() << " pair " << shape.pair << " top "
                        << shape.top << " h " << hp << "/" << ht << "\n"
                        << rep.summary());
            if (rep.ok()) ++accepted;
            if (len != W{0}) continue;
            REQUIRE_FALSE(rep.ok());
            REQUIRE(only_h3_h4(rep));
            const NodeId lcn = lowest_common_node(h, cyc);
            if (lcn == kNoNode) {
              REQUIRE(rep.has(ViolationKind::kRootArcCycle));
            } else {
              REQUIRE(has_light_cycle_at(rep, lcn));
            }
          }
        }
      }
      // Control: the same cycle with positive length 5 > every height admits hierarchies, so
      // the rejections above are caused by the zero length and not by the enumeration.
      INFO("len " << len << " cycle size " << cyc.size());
      if (len == W{0}) {
        REQUIRE(accepted == 0);
      } else {
        REQUIRE(accepted > 0);
      }
    }
  }
}

TEMPLATE_TEST_CASE("Q-21: injecting a zero-length cycle into a valid reference hierarchy breaks it",
                   "[paper][hierarchy][def41][q21][random]", std::uint64_t, double) {
  using W = TestType;
  std::size_t same_tree = 0;
  std::size_t across_trees = 0;
  for (std::uint64_t seed = 1; seed <= 300; ++seed) {
    std::mt19937_64 rng(seed);
    const std::size_t n = 2 + rng() % 11;
    const std::size_t m = n + rng() % (3 * n + 1);
    StateGraph<W> g;
    g.num_states = n;
    for (std::size_t i = 0; i < m; ++i) {
      g.arcs.push_back({static_cast<StateId>(rng() % n), static_cast<StateId>(rng() % n),
                        static_cast<W>(1 + rng() % 6)});
    }
    const std::size_t n_phys = 1 + rng() % n;
    std::vector<Vertex> pi(n);
    for (auto& p : pi) p = static_cast<Vertex>(rng() % n_phys);
    auto h = sssp::paper::build_threshold_scc_hierarchy<W>(g, pi);
    INFO("seed " << seed << " n " << n << " m " << m);
    {
      const auto rep =
          sssp::paper::validate_ordered_hierarchy(h, own_physical_graph(g, pi, n_phys));
      INFO(rep.summary());
      REQUIRE(rep.ok());  // precondition: the uncorrupted hierarchy is valid
    }
    // A zero-length cycle through k >= 2 distinct random states. Its arcs become zero-length
    // physical edges, so H1 holds for them and every other walk only gets shorter (H1, H5 keep
    // holding).
    std::vector<StateId> perm(n);
    for (StateId x = 0; x < n; ++x) perm[x] = x;
    std::shuffle(perm.begin(), perm.end(), rng);
    const std::size_t k = 2 + rng() % std::min<std::size_t>(n - 1, 3);
    const std::vector<StateId> cyc(perm.begin(), perm.begin() + static_cast<std::ptrdiff_t>(k));
    for (std::size_t i = 0; i < k; ++i) h.graph.arcs.push_back({cyc[i], cyc[(i + 1) % k], W{0}});
    const auto rep =
        sssp::paper::validate_ordered_hierarchy(h, own_physical_graph(h.graph, pi, n_phys));
    INFO("cycle length " << k << "\n" << rep.summary());
    REQUIRE_FALSE(rep.ok());
    REQUIRE(only_h3_h4(rep));
    const NodeId lcn = lowest_common_node(h, cyc);
    if (lcn == kNoNode) {
      ++across_trees;
      REQUIRE(rep.has(ViolationKind::kRootArcCycle));
    } else {
      ++same_tree;
      REQUIRE(has_light_cycle_at(rep, lcn));
    }
  }
  // Both branches of the expectation are exercised.
  REQUIRE(same_tree >= 50);
  REQUIRE(across_trees >= 50);
}
