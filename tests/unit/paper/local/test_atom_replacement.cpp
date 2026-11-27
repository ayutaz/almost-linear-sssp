// Unit tests for Sec 5.5 (txt l.1302-1364): capping a local hierarchy at the cutoff of a
// replaceable atom (Def 5.11, Lemma 5.13) and substituting the capped hierarchies into a
// partition hierarchy (Def 5.12, the replacement procedure, Lemma 5.14-5.15).
//
// Expected values come from the paper's procedure applied by hand: Alg 4's base case on a path
// is the Kruskal merge tree (Sec 2 l.286-290); capping keeps the maximal subtrees of height
// strictly below c and gives several of them a common root of height c whose links come from
// zero-length copy certificates and the atom's spanning tree (proof of Lemma 5.13).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/atom_replacement.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/paper/eval_oracle.hpp"
#include "unit/paper/local/atom_test_support.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::EdgeId;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::cap_local_hierarchy;
using sssp::paper::kNoArc;
using sssp::paper::kNoAtom;
using sssp::paper::kNoNode;
using sssp::paper::LocalHierarchyParams;
using sssp::paper::NodeId;
using sssp::paper::OrderedHierarchy;
using sssp::paper::PartitionHierarchy;
using sssp::paper::StateId;
using sssp::paper::substitute_atoms;
namespace at = sssp::test::atom;
namespace lt = sssp::test::local;
namespace tp = sssp::test::paper;

namespace {

template <class W>
CsrGraph<W> graph_of(std::size_t n, std::vector<Edge<W>> e) {
  return CsrGraph<W>::from_edges(n, std::move(e));
}

/// The projections (H vertices) of the states below node v, as a set.
template <class W>
std::set<Vertex> images_below(const OrderedHierarchy<W>& h, NodeId v) {
  std::set<Vertex> out;
  std::vector<NodeId> stack{v};
  while (!stack.empty()) {
    const NodeId u = stack.back();
    stack.pop_back();
    if (h.is_leaf(u)) {
      out.insert(h.projection[h.nodes[u].state]);
    } else {
      for (const NodeId c : h.nodes[u].children) stack.push_back(c);
    }
  }
  return out;
}

/// The states below node v, as a set.
template <class W>
std::set<StateId> states_set(const OrderedHierarchy<W>& h, NodeId v) {
  std::set<StateId> out;
  std::vector<NodeId> stack{v};
  while (!stack.empty()) {
    const NodeId u = stack.back();
    stack.pop_back();
    if (h.is_leaf(u)) {
      out.insert(h.nodes[u].state);
    } else {
      for (const NodeId c : h.nodes[u].children) stack.push_back(c);
    }
  }
  return out;
}

template <class W>
W max_internal_height(const OrderedHierarchy<W>& h) {
  W m{0};
  for (NodeId v = 0; v < h.nodes.size(); ++v) {
    if (!h.is_leaf(v)) m = std::max(m, h.nodes[v].height);
  }
  return m;
}

}  // namespace

// ------------------------------------------------------------------ Lemma 5.13

TEMPLATE_TEST_CASE("cap: maximal subtrees below c get a common root of height c", "", std::uint64_t,
                   double) {
  using W = TestType;
  // Path 0-1-2-3-4 with lengths 1, 5, 1, 2. Base case (Kruskal, equal weights merged together):
  // A = {0,1} at 1, B = {2,3} at 1, C = {B,4} at 2, R = {A,C} at 5. At c = 3 the maximal
  // subtrees below c are A and C; R (height 5 >= 3) is removed and a new root of height 3 joins
  // A and C with the tree edge 1-2 (length 5 < 2c = 6) as its link.
  const auto g = graph_of<W>(5, {{0, 1, W{1}}, {1, 2, W{5}}, {2, 3, W{1}}, {3, 4, W{2}}});
  const std::vector<EdgeId> tree{0, 1, 2, 3};
  Diagnostics diag;
  const auto capped = cap_local_hierarchy<W>(g, tree, W{3}, &diag);
  const auto& h = capped.local.hierarchy;

  REQUIRE(h.roots.size() == 1);
  REQUIRE(h.roots[0] == capped.root);
  CHECK(capped.added_root);
  CHECK(capped.removed_nodes == 1);
  const auto& root = h.nodes[capped.root];
  CHECK(root.height == W{3});
  REQUIRE(root.children.size() == 2);
  std::set<std::set<Vertex>> parts;
  for (const NodeId c : root.children) parts.insert(images_below(h, c));
  CHECK(parts == std::set<std::set<Vertex>>{{0, 1}, {2, 3, 4}});
  REQUIRE(root.links.size() == 1);
  const auto& link = root.links[0];
  CHECK(link.length == W{5});
  CHECK(link.arc == kNoArc);
  CHECK(std::set<Vertex>{h.projection[link.x], h.projection[link.y]} == std::set<Vertex>{1, 2});
  CHECK(max_internal_height(h) == W{3});
  CHECK(h.c0 == 2);
  // Native entries are kept: one state per vertex in the base case.
  REQUIRE(capped.local.native_entry.size() == 5);
  for (Vertex v = 0; v < 5; ++v) CHECK(h.projection[capped.local.native_entry[v]] == v);
  const auto report = sssp::paper::validate_ordered_hierarchy(h, g);
  INFO(report.summary());
  CHECK(report.ok());
  // Alg 4 took the paper's base case on k = 5 <= 2 L_5^4, which is recorded.
  REQUIRE(diag.fallbacks.size() == 1);
  CHECK(diag.fallbacks[0].reason ==
        sssp::paper::local_base_case_reason(sssp::paper::kPaperLocalParams));
}

TEMPLATE_TEST_CASE("cap: a cutoff above every height keeps Alg 4's tree and adds no root", "",
                   std::uint64_t, double) {
  using W = TestType;
  const auto g = graph_of<W>(5, {{0, 1, W{1}}, {1, 2, W{5}}, {2, 3, W{1}}, {3, 4, W{2}}});
  const std::vector<EdgeId> tree{0, 1, 2, 3};
  const auto local = sssp::paper::local_hierarchy(g);
  const auto capped = cap_local_hierarchy<W>(g, tree, W{6});
  const auto& h = capped.local.hierarchy;
  REQUIRE(h.roots.size() == 1);
  CHECK_FALSE(capped.added_root);
  CHECK(capped.removed_nodes == 0);
  CHECK(h.nodes[capped.root].height == W{5});
  CHECK(h.nodes.size() == local.hierarchy.nodes.size());
  CHECK(h.c0 == local.hierarchy.c0);
  CHECK(sssp::paper::validate_ordered_hierarchy(h, g).ok());
}

TEMPLATE_TEST_CASE("cap: a cutoff at the least height keeps only the leaves under the new root", "",
                   std::uint64_t, double) {
  using W = TestType;
  // All lengths 1: Alg 4 has one node of height 1, which is not strictly below c = 1.
  const auto g = graph_of<W>(4, {{0, 1, W{1}}, {1, 2, W{1}}, {2, 3, W{1}}});
  const std::vector<EdgeId> tree{0, 1, 2};
  const auto capped = cap_local_hierarchy<W>(g, tree, W{1});
  const auto& h = capped.local.hierarchy;
  REQUIRE(h.roots.size() == 1);
  CHECK(capped.added_root);
  CHECK(capped.removed_nodes == 1);
  const auto& root = h.nodes[capped.root];
  CHECK(root.height == W{1});
  CHECK(root.children.size() == 4);
  REQUIRE(root.links.size() == 3);
  for (const auto& l : root.links) CHECK(l.length < W{2});
  CHECK(sssp::paper::validate_ordered_hierarchy(h, g).ok());
}

TEMPLATE_TEST_CASE("cap: a one-state atom is its own root", "", std::uint64_t, double) {
  using W = TestType;
  const auto g = graph_of<W>(1, {});
  const auto capped = cap_local_hierarchy<W>(g, {}, W{5});
  const auto& h = capped.local.hierarchy;
  REQUIRE(h.num_states() == 1);
  REQUIRE(h.roots.size() == 1);
  CHECK(h.is_leaf(capped.root));
  CHECK_FALSE(capped.added_root);
  CHECK(capped.local.native_entry == std::vector<StateId>{0});
}

TEMPLATE_TEST_CASE("cap: rejects inputs outside Def 5.11", "", std::uint64_t, double) {
  using W = TestType;
  const auto g = graph_of<W>(5, {{0, 1, W{1}}, {1, 2, W{5}}, {2, 3, W{1}}, {3, 4, W{2}}});
  const std::vector<EdgeId> tree{0, 1, 2, 3};
  // c > 0.
  CHECK_THROWS_AS(cap_local_hierarchy<W>(g, tree, W{0}), std::invalid_argument);
  CHECK_THROWS_AS(cap_local_hierarchy<W>(g, tree, WeightTraits<W>::infinity()),
                  std::invalid_argument);
  if constexpr (std::is_floating_point_v<W>) {
    CHECK_THROWS_AS(cap_local_hierarchy<W>(g, tree, std::numeric_limits<W>::quiet_NaN()),
                    std::invalid_argument);
    CHECK_THROWS_AS(cap_local_hierarchy<W>(g, tree, W{-1}), std::invalid_argument);
  }
  // Tree edge 1-2 has length 5, not below 2c = 4.
  CHECK_THROWS_AS(cap_local_hierarchy<W>(g, tree, W{2}), std::invalid_argument);
  // Not a spanning tree: too few edges, a repeated edge, an id out of range.
  CHECK_THROWS_AS(cap_local_hierarchy<W>(g, std::vector<EdgeId>{0, 1, 2}, W{3}),
                  std::invalid_argument);
  CHECK_THROWS_AS(cap_local_hierarchy<W>(g, std::vector<EdgeId>{0, 1, 2, 2}, W{3}),
                  std::invalid_argument);
  CHECK_THROWS_AS(cap_local_hierarchy<W>(g, std::vector<EdgeId>{0, 1, 2, 9}, W{3}),
                  std::invalid_argument);
  // Disconnected internal graph, empty atom.
  CHECK_THROWS_AS(cap_local_hierarchy<W>(graph_of<W>(2, {}), {}, W{3}), std::invalid_argument);
  CHECK_THROWS_AS(cap_local_hierarchy<W>(graph_of<W>(0, {}), {}, W{3}), std::invalid_argument);
}

TEMPLATE_TEST_CASE("cap: on the radial branch the capped hierarchy is rooted, below c, valid", "",
                   std::uint64_t, double) {
  using W = TestType;
  std::size_t added = 0;
  std::size_t kept = 0;
  for (const auto f : {lt::Family::kRandomTree, lt::Family::kSparse, lt::Family::kGrid}) {
    for (const auto wd : {lt::Weights::kTies, lt::Weights::kUniform, lt::Weights::kExponential}) {
      for (const std::size_t b : {std::size_t{2}, std::size_t{6}}) {
        const std::size_t n = 40;
        const auto edges = lt::make_graph<W>(f, wd, n, 7 * n + b);
        const auto g = CsrGraph<W>::from_edges(n, edges);
        const LocalHierarchyParams p{1, 1, b};
        const auto tree = at::kruskal_tree<W>(n, edges);
        W tree_max{0};
        for (const EdgeId id : tree) tree_max = std::max(tree_max, edges[id].weight);
        // Independent reference: the uncapped Alg 4 result (deterministic). Candidate cutoffs:
        // a random admissible value, every height h, and h + 1.
        const auto local = sssp::paper::local_hierarchy(g, nullptr, p);
        const auto& u = local.hierarchy;
        sssp::test::SplitMix64 rng(b);
        std::vector<W> cuts{at::draw_cutoff<W>(at::half_of(tree_max), tree_max + W{1}, rng)};
        for (NodeId v = 0; v < u.nodes.size(); ++v) {
          if (!u.is_leaf(v) && at::half_of(tree_max) < u.nodes[v].height) {
            // c = h removes the nodes of height h; c = h + 1 keeps them.
            cuts.push_back(u.nodes[v].height);
            cuts.push_back(u.nodes[v].height + W{1});
          }
        }
        for (const W c : cuts) {
          if (!(tree_max < c + c)) continue;
          INFO(lt::name(f) << " " << lt::name(wd) << " b=" << b << " c=" << c);
          Diagnostics diag;
          const auto capped = cap_local_hierarchy<W>(g, tree, c, &diag, p);
          CHECK(diag.fallbacks.empty());
          const auto& h = capped.local.hierarchy;
          // Capping keeps Alg 4's states, arcs, projection and native entries.
          REQUIRE(h.num_states() == u.num_states());
          REQUIRE(h.graph.arcs.size() == u.graph.arcs.size());
          CHECK(h.projection == u.projection);
          CHECK(capped.local.native_entry == local.native_entry);
          REQUIRE(h.roots.size() == 1);
          CHECK(h.roots[0] == capped.root);
          for (NodeId v = 0; v < h.nodes.size(); ++v) {
            if (!h.is_leaf(v)) CHECK(!(c < h.nodes[v].height));
          }
          for (const auto& l : h.nodes[capped.root].links) CHECK(l.length < c + c);
          const auto report = sssp::paper::validate_ordered_hierarchy(h, g);
          INFO(report.summary());
          CHECK(report.ok());
          // The children of an added root are exactly the maximal subtrees of the uncapped
          // hierarchy of height strictly below c (leaves included).
          std::set<std::set<StateId>> want;
          for (NodeId v = 0; v < u.nodes.size(); ++v) {
            const bool low = u.is_leaf(v) || u.nodes[v].height < c;
            const NodeId par = u.nodes[v].parent;
            const bool parent_high = par == kNoNode || !(u.nodes[par].height < c);
            if (low && parent_high) want.insert(states_set(u, v));
          }
          if (want.size() == 1) {
            CHECK_FALSE(capped.added_root);
            CHECK(states_set(h, capped.root) == *want.begin());
            ++kept;
          } else {
            REQUIRE(capped.added_root);
            CHECK(h.nodes[capped.root].height == c);
            std::set<std::set<StateId>> got;
            for (const NodeId ch : h.nodes[capped.root].children) got.insert(states_set(h, ch));
            CHECK(got == want);
            ++added;
          }
        }
      }
    }
  }
  // Alg 4's radial result has one root per SCC of the local state graph (Sec 5.4 l.1280), and
  // on these inputs there are always several, so every cutoff adds a common root here; the
  // single-retained case is covered by the base-case tests above.
  CHECK(added > 0);
  CHECK(kept == 0);
}

// ------------------------------------------------------------------ Lemma 5.14-5.15

TEMPLATE_TEST_CASE("substitute: hand example keeps the outer forest and maps arcs and links", "",
                   std::uint64_t, double) {
  using W = TestType;
  // States 0..4, image = id. Atom {1, 2, 3} with internal path 1-2 (1), 2-3 (2), cutoff 3.
  // D: 0<->1 (5), 3<->4 (5), one-way 4->1 (9). One root R of height 5 over {0}, atom, {4}.
  PartitionHierarchy<W> ph;
  ph.graph.num_states = 5;
  ph.graph.arcs = {{0, 1, W{5}}, {1, 0, W{5}}, {3, 4, W{5}}, {4, 3, W{5}}, {4, 1, W{9}}};
  ph.projection = {0, 1, 2, 3, 4};
  ph.atoms.push_back({{1, 2, 3}, {{0, 1, W{1}}, {1, 2, W{2}}}, {0, 1}, W{3}});
  ph.nodes.resize(4);
  ph.leaf_atom = {kNoAtom, 0, kNoAtom, kNoAtom};
  ph.nodes[0].state = 0;
  ph.nodes[2].state = 4;
  for (const NodeId v : {NodeId{0}, NodeId{1}, NodeId{2}}) ph.nodes[v].parent = 3;
  ph.nodes[3].children = {0, 1, 2};
  ph.nodes[3].height = W{5};
  ph.nodes[3].links = {{0, 1, W{5}, 0}, {3, 4, W{5}, 2}};
  ph.c0 = 1;
  const auto physical =
      graph_of<W>(5, {{0, 1, W{5}}, {1, 2, W{1}}, {2, 3, W{2}}, {3, 4, W{5}}, {4, 1, W{9}}});

  Diagnostics diag;
  const auto sub = substitute_atoms(ph, &diag);
  const auto& h = sub.hierarchy;
  // Base case: one copy per state, so 5 states; 5 copied D arcs + 4 internal arcs.
  REQUIRE(h.num_states() == 5);
  REQUIRE(h.graph.arcs.size() == 9);
  REQUIRE(sub.entry.size() == 5);
  REQUIRE(sub.old_state.size() == 5);
  for (StateId x = 0; x < 5; ++x) {
    CHECK(sub.old_state[sub.entry[x]] == x);
    CHECK(h.projection[sub.entry[x]] == x);
  }
  // Every D arc is copied once (one tail copy each), with its head at the entry of its head.
  for (sssp::paper::ArcId a = 0; a < ph.graph.arcs.size(); ++a) {
    std::size_t copies = 0;
    for (sssp::paper::ArcId b = 0; b < h.graph.arcs.size(); ++b) {
      if (sub.old_arc[b] != a) continue;
      ++copies;
      CHECK(sub.old_state[h.graph.arcs[b].tail] == ph.graph.arcs[a].tail);
      CHECK(h.graph.arcs[b].head == sub.entry[ph.graph.arcs[a].head]);
      CHECK(h.graph.arcs[b].length == ph.graph.arcs[a].length);
    }
    CHECK(copies == 1);
  }
  // Atom {1,2,3}: base-case tree {1,2} at 1, then {.,3} at 2 < c = 3: kept, no new root.
  REQUIRE(sub.atom_root.size() == 1);
  CHECK(h.nodes[sub.atom_root[0]].height == W{2});
  // The outer root keeps its height, its three children and its links on the entries.
  REQUIRE(h.roots.size() == 1);
  const auto& r = h.nodes[h.roots[0]];
  CHECK(r.height == W{5});
  REQUIRE(r.children.size() == 3);
  CHECK(std::find(r.children.begin(), r.children.end(), sub.atom_root[0]) != r.children.end());
  REQUIRE(r.links.size() == 2);
  CHECK(r.links[0].x == sub.entry[0]);
  CHECK(r.links[0].y == sub.entry[1]);
  CHECK(r.links[1].x == sub.entry[3]);
  CHECK(r.links[1].y == sub.entry[4]);
  for (const auto& l : r.links) {
    REQUIRE(l.arc != kNoArc);
    const auto& arc = h.graph.arcs[l.arc];
    CHECK(std::set<StateId>{arc.tail, arc.head} == std::set<StateId>{l.x, l.y});
  }
  REQUIRE(sub.partition_node.size() == ph.nodes.size());
  CHECK(sub.partition_node[3] == h.roots[0]);
  CHECK(sub.partition_node[1] == sub.atom_root[0]);
  const auto report = sssp::paper::validate_ordered_hierarchy(h, physical);
  INFO(report.summary());
  CHECK(report.ok());
  // One Alg 4 call, in the base case, recorded.
  REQUIRE(diag.fallbacks.size() == 1);
  CHECK(diag.fallbacks[0].vertices == 3);

  // Seeds go to the native copy; Evaluate from state 0 gives 0, 5, 6, 8, 13 (path 0-1-2-3-4;
  // the one-way 4->1 and the physical edge 1-4 of length 9 are longer).
  const std::vector<sssp::paper::StateOffer<W>> offers{{0, W{0}}};
  const auto moved = sssp::paper::redirect_offers(sub, std::span{offers});
  REQUIRE(moved.size() == 1);
  CHECK(moved[0].state == sub.entry[0]);
  const sssp::paper::SeededOrderedHierarchy<W> seeded{h, {{0, W{0}}}, moved};
  REQUIRE(sssp::paper::validate_seeded_ordered_hierarchy(seeded, physical).ok());
  const auto ev = sssp::paper::evaluate(seeded);
  const auto got = tp::physical_minima<W>(5, h.projection, ev.emissions);
  CHECK(got == std::vector<W>{W{0}, W{5}, W{6}, W{8}, W{13}});
  const std::vector<sssp::paper::StateOffer<W>> bad{{7, W{0}}};
  CHECK_THROWS_AS(sssp::paper::redirect_offers(sub, std::span{bad}), std::out_of_range);
}

TEMPLATE_TEST_CASE("substitute: every copy of an old tail carries the old exits", "", std::uint64_t,
                   double) {
  using W = TestType;
  // Groups of 37..44 states take the radial branch with threshold k > L_k^2 (k >= 37); smaller
  // atoms take the base case. Radial atoms have several copies of some states.
  std::size_t multi_copy = 0;
  for (std::uint64_t seed = 1; seed <= 4; ++seed) {
    at::PartitionParams pp;
    pp.num_states = 90;
    pp.min_group = 37;
    pp.max_group = 44;
    pp.internal_max = 1000;
    pp.outer_min = 2000;
    pp.outer_max = 5000;
    pp.one_way = 6;
    pp.largest_first = true;
    pp.atom_percent = 100;
    const auto inst = at::make_partition_instance<W>(pp, seed);
    const auto& ph = inst.partition;
    const LocalHierarchyParams p{1, 2, 3};
    const auto sub = substitute_atoms(ph, nullptr, p);
    const auto& h = sub.hierarchy;
    INFO("seed " << seed);
    std::vector<std::size_t> copies(ph.graph.num_states, 0);
    for (StateId t = 0; t < h.num_states(); ++t) {
      ++copies[sub.old_state[t]];
      CHECK(h.projection[t] == ph.projection[sub.old_state[t]]);
    }
    for (StateId x = 0; x < ph.graph.num_states; ++x) {
      REQUIRE(copies[x] >= 1);
      CHECK(sub.old_state[sub.entry[x]] == x);
      if (copies[x] > 1) ++multi_copy;
    }
    // Old arc a = x -> y: exactly one copy per copy of x, each entering the entry of y.
    std::vector<std::vector<StateId>> tails(ph.graph.arcs.size());
    for (sssp::paper::ArcId b = 0; b < h.graph.arcs.size(); ++b) {
      const auto a = sub.old_arc[b];
      if (a == kNoArc) continue;
      tails[a].push_back(h.graph.arcs[b].tail);
      CHECK(h.graph.arcs[b].head == sub.entry[ph.graph.arcs[a].head]);
      CHECK(h.graph.arcs[b].length == ph.graph.arcs[a].length);
    }
    for (sssp::paper::ArcId a = 0; a < ph.graph.arcs.size(); ++a) {
      auto& ts = tails[a];
      std::sort(ts.begin(), ts.end());
      CHECK(std::adjacent_find(ts.begin(), ts.end()) == ts.end());
      CHECK(ts.size() == copies[ph.graph.arcs[a].tail]);
      for (const StateId t : ts) CHECK(sub.old_state[t] == ph.graph.arcs[a].tail);
    }
    const auto report = sssp::paper::validate_ordered_hierarchy(h, inst.physical);
    INFO(report.summary());
    CHECK(report.ok());
  }
  CHECK(multi_copy > 0);
}

TEMPLATE_TEST_CASE("substitute: rejects partition hierarchies outside Def 5.12", "", std::uint64_t,
                   double) {
  using W = TestType;
  const auto base = [] {
    PartitionHierarchy<W> ph;
    ph.graph.num_states = 3;
    ph.graph.arcs = {{0, 1, W{4}}, {1, 0, W{4}}};
    ph.projection = {0, 1, 2};
    ph.atoms.push_back({{1, 2}, {{0, 1, W{1}}}, {0}, W{2}});
    ph.nodes.resize(3);
    ph.leaf_atom = {kNoAtom, 0, kNoAtom};
    ph.nodes[0].state = 0;
    ph.nodes[0].parent = 2;
    ph.nodes[1].parent = 2;
    ph.nodes[2].children = {0, 1};
    ph.nodes[2].height = W{4};
    ph.nodes[2].links = {{0, 1, W{4}, 0}};
    return ph;
  };
  CHECK_NOTHROW(substitute_atoms(base()));
  {
    auto ph = base();  // an arc inside an atom belongs to its internal graph H
    ph.graph.arcs.push_back({1, 2, W{1}});
    CHECK_THROWS_AS(substitute_atoms(ph), std::invalid_argument);
  }
  {
    auto ph = base();  // cutoff above the parent's height (Def 5.11, Def 5.12)
    ph.atoms[0].cutoff = W{5};
    CHECK_THROWS_AS(substitute_atoms(ph), std::invalid_argument);
  }
  {
    auto ph = base();  // state 2 in two parts
    ph.atoms[0].states = {1, 2};
    ph.nodes[0].state = 2;
    CHECK_THROWS_AS(substitute_atoms(ph), std::invalid_argument);
  }
  {
    auto ph = base();  // state 0 in no part
    ph.graph.num_states = 4;
    ph.projection.push_back(3);
    CHECK_THROWS_AS(substitute_atoms(ph), std::invalid_argument);
  }
  {
    auto ph = base();  // atom without a leaf
    ph.atoms.push_back({{}, {}, {}, W{1}});
    CHECK_THROWS_AS(substitute_atoms(ph), std::invalid_argument);
  }
  {
    auto ph = base();  // leaf_atom has the wrong size
    ph.leaf_atom.pop_back();
    CHECK_THROWS_AS(substitute_atoms(ph), std::invalid_argument);
  }
  {
    auto ph = base();  // parent pointer and child list disagree
    ph.nodes[0].parent = kNoNode;
    CHECK_THROWS_AS(substitute_atoms(ph), std::invalid_argument);
  }
}
