// Property tests for Sec 5.5 (txt l.1302-1364) on random partition hierarchies (Def 5.12) built
// by tests/unit/paper/local/atom_test_support.hpp from threshold-SCC reference hierarchies:
//   * Lemma 5.15: the substituted forest validates Def 4.1 (c0 <= 2) on the physical graph;
//   * Lemma 5.13: below every atom's root lie exactly the copies of the atom's states, every
//     height is at most the cutoff c, and every link at the atom's root is shorter than 2c;
//   * replacement procedure: every copy keeps its old state's image, seeds and entering arcs go
//     to one copy of their state, and every copy of an old tail carries each old exit once;
//   * Lemma 5.14 at the level of the state graphs (independent of Evaluate): for ANY initial
//     values, the least state-graph distance over the copies of an old state equals its distance
//     in the underlying state graph (D plus both orientations of every internal edge);
//   * the outer forest is kept: each partition node maps to a node with the same height, child
//     count and link count; atom leaves map to the atom roots;
//   * fallback policy: one FallbackRecord per atom that takes Alg 4's base case, none otherwise.
// Weights: std::uint64_t and quarter-integer doubles (exact sums).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/atom_replacement.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/paper/eval_oracle.hpp"
#include "support/prng.hpp"
#include "unit/paper/local/atom_test_support.hpp"

using sssp::Diagnostics;
using sssp::WeightTraits;
using sssp::paper::ArcId;
using sssp::paper::kNoArc;
using sssp::paper::kNoAtom;
using sssp::paper::LocalHierarchyParams;
using sssp::paper::NodeId;
using sssp::paper::StateId;
namespace at = sssp::test::atom;
namespace tp = sssp::test::paper;

namespace {

struct Family {
  const char* name;
  at::PartitionParams params;
  LocalHierarchyParams local;
  std::uint64_t seeds;
};

std::vector<Family> families() {
  std::vector<Family> f;
  {
    at::PartitionParams p;  // small groups, undirected, identity projection
    f.push_back({"small-undirected", p, sssp::paper::kPaperLocalParams, 25});
  }
  {
    at::PartitionParams p;  // one-way arcs, folded projection
    p.num_states = 36;
    p.one_way = 12;
    p.fold_projection = true;
    p.singleton_percent = 25;
    f.push_back({"one-way-folded", p, sssp::paper::kPaperLocalParams, 25});
  }
  {
    at::PartitionParams p;  // interleaved lengths: atoms span several groups
    p.num_states = 40;
    p.internal_max = 30;
    p.outer_min = 1;
    p.outer_max = 30;
    p.cross_edges = 15;
    p.one_way = 5;
    p.chords_per_group = 3;
    f.push_back({"interleaved", p, sssp::paper::kPaperLocalParams, 25});
  }
  {
    at::PartitionParams p;  // disconnected groups: several roots
    p.num_states = 30;
    p.connect = false;
    p.cross_edges = 2;
    p.one_way = 3;
    f.push_back({"forest", p, sssp::paper::kPaperLocalParams, 15});
  }
  for (const bool exponential : {false, true}) {
    for (const std::size_t b : {std::size_t{2}, std::size_t{5}}) {
      at::PartitionParams p;  // groups of >= 37 states: the radial branch with k > L_k^2
      p.num_states = 100;
      p.min_group = 37;
      p.max_group = 50;
      p.internal_max = 1000;
      p.internal_exponential = exponential;
      p.outer_min = exponential ? (std::uint64_t{1} << 42) : 2000;
      p.outer_max = exponential ? (std::uint64_t{1} << 43) : 5000;
      p.one_way = 8;
      p.largest_first = true;
      p.atom_percent = 100;
      f.push_back({exponential ? "radial-exponential" : "radial-uniform", p,
                   LocalHierarchyParams{1, 2, b}, 3});
    }
  }
  return f;
}

template <class W>
void check_instance(const at::PartitionInstance<W>& inst, const LocalHierarchyParams& lp,
                    std::uint64_t seed, std::size_t& radial_atoms, std::size_t& added_roots) {
  const auto& ph = inst.partition;
  const std::size_t n = ph.graph.num_states;
  Diagnostics diag;
  const auto sub = sssp::paper::substitute_atoms(ph, &diag, lp);
  const auto& h = sub.hierarchy;

  // Lemma 5.15: Def 4.1 on the substituted forest.
  CHECK(h.c0 <= 2);
  const auto report = sssp::paper::validate_ordered_hierarchy(h, inst.physical);
  INFO(report.summary());
  REQUIRE(report.ok());

  // Copies, images, entries.
  REQUIRE(sub.old_state.size() == h.num_states());
  REQUIRE(sub.entry.size() == n);
  std::vector<std::size_t> copies(n, 0);
  for (StateId t = 0; t < h.num_states(); ++t) {
    REQUIRE(sub.old_state[t] < n);
    ++copies[sub.old_state[t]];
    CHECK(h.projection[t] == ph.projection[sub.old_state[t]]);
  }
  std::vector<std::uint32_t> atom_of(n, kNoAtom);
  for (std::uint32_t a = 0; a < ph.atoms.size(); ++a) {
    for (const StateId x : ph.atoms[a].states) atom_of[x] = a;
  }
  for (StateId x = 0; x < n; ++x) {
    REQUIRE(copies[x] >= 1);
    if (atom_of[x] == kNoAtom) CHECK(copies[x] == 1);
    REQUIRE(sub.entry[x] < h.num_states());
    CHECK(sub.old_state[sub.entry[x]] == x);
  }

  // Exits on every copy.
  REQUIRE(sub.old_arc.size() == h.graph.arcs.size());
  std::vector<std::vector<StateId>> tails(ph.graph.arcs.size());
  for (ArcId b = 0; b < h.graph.arcs.size(); ++b) {
    const auto& arc = h.graph.arcs[b];
    const ArcId a = sub.old_arc[b];
    if (a == kNoArc) {
      // An internal arc of an atom's local hierarchy: both ends are copies of that atom.
      const std::uint32_t at_tail = atom_of[sub.old_state[arc.tail]];
      REQUIRE(at_tail != kNoAtom);
      CHECK(atom_of[sub.old_state[arc.head]] == at_tail);
      continue;
    }
    REQUIRE(a < ph.graph.arcs.size());
    tails[a].push_back(arc.tail);
    CHECK(arc.head == sub.entry[ph.graph.arcs[a].head]);
    CHECK(arc.length == ph.graph.arcs[a].length);
  }
  for (ArcId a = 0; a < ph.graph.arcs.size(); ++a) {
    auto& ts = tails[a];
    std::sort(ts.begin(), ts.end());
    CHECK(std::adjacent_find(ts.begin(), ts.end()) == ts.end());
    CHECK(ts.size() == copies[ph.graph.arcs[a].tail]);
    for (const StateId t : ts) CHECK(sub.old_state[t] == ph.graph.arcs[a].tail);
  }

  // Lemma 5.13 per atom.
  REQUIRE(sub.atom_root.size() == ph.atoms.size());
  for (std::uint32_t a = 0; a < ph.atoms.size(); ++a) {
    const auto& atom = ph.atoms[a];
    const W c = atom.cutoff;
    std::vector<NodeId> stack{sub.atom_root[a]};
    std::size_t leaves = 0;
    while (!stack.empty()) {
      const NodeId v = stack.back();
      stack.pop_back();
      if (h.is_leaf(v)) {
        ++leaves;
        CHECK(atom_of[sub.old_state[h.nodes[v].state]] == a);
      } else {
        CHECK(!(c < h.nodes[v].height));
        for (const NodeId ch : h.nodes[v].children) stack.push_back(ch);
      }
    }
    std::size_t want = 0;
    for (const StateId x : atom.states) want += copies[x];
    CHECK(leaves == want);
    for (const auto& l : h.nodes[sub.atom_root[a]].links) CHECK(l.length < c + c);
    if (!h.is_leaf(sub.atom_root[a]) && h.nodes[sub.atom_root[a]].height == c) ++added_roots;
    if (sssp::paper::local_construction_applies(atom.states.size(), lp)) ++radial_atoms;
  }

  // The outer forest is kept.
  REQUIRE(sub.partition_node.size() == ph.nodes.size());
  for (NodeId v = 0; v < ph.nodes.size(); ++v) {
    const NodeId u = sub.partition_node[v];
    REQUIRE(u < h.nodes.size());
    if (ph.leaf_atom[v] != kNoAtom) {
      CHECK(u == sub.atom_root[ph.leaf_atom[v]]);
    } else if (ph.nodes[v].children.empty()) {
      REQUIRE(h.is_leaf(u));
      CHECK(h.nodes[u].state == sub.entry[ph.nodes[v].state]);
    } else {
      CHECK(h.nodes[u].height == ph.nodes[v].height);
      CHECK(h.nodes[u].children.size() == ph.nodes[v].children.size());
      CHECK(h.nodes[u].links.size() == ph.nodes[v].links.size());
    }
  }
  CHECK(h.roots.size() == [&] {
    std::size_t r = 0;
    for (const auto& node : ph.nodes) r += node.parent == sssp::paper::kNoNode ? 1 : 0;
    return r;
  }());

  // Fallback policy: one record per base-case atom.
  std::size_t base_atoms = 0;
  for (const auto& atom : ph.atoms) {
    base_atoms += sssp::paper::local_construction_applies(atom.states.size(), lp) ? 0 : 1;
  }
  CHECK(diag.fallbacks.size() == base_atoms);

  // Lemma 5.14 on the state graphs, for arbitrary initial values.
  sssp::test::SplitMix64 rng(seed * 977 + 3);
  for (int round = 0; round < 3; ++round) {
    std::vector<sssp::paper::StateOffer<W>> offers;
    const std::size_t k = 1 + rng.below(3);
    for (std::size_t i = 0; i < k; ++i) {
      offers.push_back({static_cast<StateId>(rng.below(n)), static_cast<W>(rng.between(0, 20))});
    }
    const auto before = tp::state_graph_distances<W>(
        inst.full_graph, tp::best_initial_offers<W>(n, std::span{offers}));
    const auto moved = sssp::paper::redirect_offers<W>(sub, offers);
    const auto after = tp::state_graph_distances<W>(
        h.graph, tp::best_initial_offers<W>(h.num_states(), std::span{moved}));
    std::vector<W> best(n, WeightTraits<W>::infinity());
    for (StateId t = 0; t < h.num_states(); ++t) {
      best[sub.old_state[t]] = std::min(best[sub.old_state[t]], after[t]);
    }
    for (StateId x = 0; x < n; ++x) {
      INFO("state " << x << " round " << round);
      CHECK(best[x] == before[x]);
    }
  }
}

}  // namespace

TEMPLATE_TEST_CASE("atom substitution: Def 4.1, Lemma 5.13 and Lemma 5.14 on random partitions", "",
                   std::uint64_t, double) {
  using W = TestType;
  std::size_t radial_atoms = 0;
  std::size_t added_roots = 0;
  std::size_t atoms = 0;
  std::size_t multi_atoms = 0;
  for (const auto& fam : families()) {
    for (std::uint64_t seed = 1; seed <= fam.seeds; ++seed) {
      INFO(fam.name << " seed " << seed);
      const auto inst = at::make_partition_instance<W>(fam.params, seed);
      atoms += inst.partition.atoms.size();
      for (const auto& a : inst.partition.atoms) multi_atoms += a.states.size() >= 2 ? 1 : 0;
      check_instance(inst, fam.local, seed, radial_atoms, added_roots);
    }
  }
  CHECK(multi_atoms > 50);
  CHECK(radial_atoms > 0);
  CHECK(added_roots > 0);
  CHECK(atoms > multi_atoms);  // one-state atoms occur as well
}
