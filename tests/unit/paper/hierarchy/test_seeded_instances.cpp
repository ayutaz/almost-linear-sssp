// Self-tests of the seeded ordered hierarchy generators (tests/support/paper/seeded_instances.hpp):
// every instance satisfies Def 4.1 / 4.2 (independent validator), generation is reproducible
// from the seed, and the families actually contain what they promise (light arcs, several
// roots, unprotected states, zero-length arcs without zero cycles, ties, several sources).

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "support/paper/eval_invariants.hpp"
#include "support/paper/eval_oracle.hpp"
#include "support/paper/seeded_instances.hpp"

using sssp::paper::kNoNode;
using sssp::paper::NodeId;
using sssp::paper::OrderedHierarchy;
using sssp::paper::StateId;
using namespace sssp::test::paper;

namespace {

template <class W>
NodeId child_towards(const OrderedHierarchy<W>& h, NodeId anc, NodeId v) {
  while (h.nodes[v].parent != anc) v = h.nodes[v].parent;
  return v;
}

template <class W>
std::vector<NodeId> ancestors(const OrderedHierarchy<W>& h, NodeId v) {
  std::vector<NodeId> a;
  for (; v != kNoNode; v = h.nodes[v].parent) a.push_back(v);
  return a;
}

/// Number of light arcs (length < h(C) between distinct children of their LCA C).
template <class W>
std::size_t count_light_arcs(const OrderedHierarchy<W>& h) {
  std::size_t k = 0;
  for (const auto& a : h.graph.arcs) {
    if (a.tail == a.head) continue;
    const auto up = ancestors(h, h.leaf_of_state[a.tail]);
    NodeId lca = kNoNode;
    for (NodeId v = h.leaf_of_state[a.head]; v != kNoNode && lca == kNoNode;
         v = h.nodes[v].parent) {
      for (const NodeId u : up) {
        if (u == v) lca = v;
      }
    }
    if (lca == kNoNode) continue;  // inter-root arc
    if (child_towards(h, lca, h.leaf_of_state[a.tail]) ==
        child_towards(h, lca, h.leaf_of_state[a.head])) {
      continue;
    }
    if (a.length < h.nodes[lca].height) ++k;
  }
  return k;
}

template <class W>
bool same_instance(const SeededInstance<W>& x, const SeededInstance<W>& y) {
  const auto& a = x.seeded;
  const auto& b = y.seeded;
  if (a.hierarchy.graph.arcs.size() != b.hierarchy.graph.arcs.size()) return false;
  for (std::size_t i = 0; i < a.hierarchy.graph.arcs.size(); ++i) {
    const auto& p = a.hierarchy.graph.arcs[i];
    const auto& q = b.hierarchy.graph.arcs[i];
    if (p.tail != q.tail || p.head != q.head || p.length != q.length) return false;
  }
  if (a.hierarchy.projection != b.hierarchy.projection) return false;
  if (a.offers.size() != b.offers.size() || a.sources.size() != b.sources.size()) return false;
  for (std::size_t i = 0; i < a.offers.size(); ++i) {
    if (a.offers[i].state != b.offers[i].state || a.offers[i].value != b.offers[i].value)
      return false;
  }
  for (std::size_t i = 0; i < a.sources.size(); ++i) {
    if (a.sources[i].vertex != b.sources[i].vertex || a.sources[i].offset != b.sources[i].offset)
      return false;
  }
  return x.physical.num_edges() == y.physical.num_edges();
}

/// True iff the length-0 arcs between distinct states contain a directed cycle.
template <class W>
bool zero_cycle(const sssp::paper::StateGraph<W>& g) {
  const std::size_t n = g.num_states;
  std::vector<std::vector<StateId>> out(n);
  for (const auto& a : g.arcs) {
    if (a.length == W{0} && a.tail != a.head) out[a.tail].push_back(a.head);
  }
  std::vector<int> color(n, 0);  // 0 new, 1 on stack, 2 done
  std::vector<std::pair<StateId, std::size_t>> st;
  for (StateId s = 0; s < n; ++s) {
    if (color[s] != 0) continue;
    st.emplace_back(s, 0);
    color[s] = 1;
    while (!st.empty()) {
      auto& [v, i] = st.back();
      if (i < out[v].size()) {
        const StateId w = out[v][i++];
        if (color[w] == 1) return true;
        if (color[w] == 0) {
          color[w] = 1;
          st.emplace_back(w, 0);
        }
      } else {
        color[v] = 2;
        st.pop_back();
      }
    }
  }
  return false;
}

}  // namespace

TEMPLATE_TEST_CASE("seeded instances: reference catalogue satisfies Def 4.1 and Def 4.2",
                   "[seeded-instances]", std::uint64_t, double) {
  using W = TestType;
  std::size_t light = 0, multi_root = 0, multi_source = 0, unprotected = 0, protected_ = 0;
  for (const auto& inst : reference_catalogue<W>(1, 300)) {
    INFO(inst.name);
    const auto rep = sssp::paper::validate_seeded_ordered_hierarchy(inst.seeded, inst.physical);
    INFO(rep.summary());
    REQUIRE(rep.ok());
    const auto ref = compute_eval_reference<W>(inst.seeded, inst.physical);
    light += count_light_arcs(inst.seeded.hierarchy);
    multi_root += inst.seeded.hierarchy.roots.size() > 1 ? 1 : 0;
    multi_source += inst.seeded.sources.size() > 1 ? 1 : 0;
    protected_ += ref.num_protected();
    unprotected += inst.seeded.hierarchy.num_states() - ref.num_protected();
  }
  CHECK(light > 100);
  CHECK(multi_root > 50);
  CHECK(multi_source > 50);
  CHECK(protected_ > 300);
  CHECK(unprotected > 300);
}

TEST_CASE("seeded instances: generation is reproducible from the seed", "[seeded-instances]") {
  ReferenceParams p;
  p.num_sources = 3;
  p.extra_physical_edges = 4;
  p.zero_arcs = 6;
  CHECK(same_instance(reference_instance<std::uint64_t>(p, 7),
                      reference_instance<std::uint64_t>(p, 7)));
  CHECK_FALSE(same_instance(reference_instance<std::uint64_t>(p, 7),
                            reference_instance<std::uint64_t>(p, 8)));
  CHECK(same_instance(zero_arc_instance<std::uint64_t>(p, 9),
                      zero_arc_instance<std::uint64_t>(p, 9)));
  const auto a = base_case_catalogue<std::uint64_t>(5);
  const auto b = base_case_catalogue<std::uint64_t>(5);
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) CHECK(same_instance(a[i], b[i]));
}

TEMPLATE_TEST_CASE("seeded instances: zero-length arcs are present but never form a cycle",
                   "[seeded-instances]", std::uint64_t, double) {
  using W = TestType;
  std::size_t zero_nonloop = 0;
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    sssp::test::SplitMix64 rng(seed);
    ReferenceParams p;
    p.num_states = 2 + rng.below(20);
    p.num_physical = 1 + rng.below(p.num_states);
    p.num_arcs = rng.below(3 * p.num_states + 1);
    p.zero_arcs = 1 + rng.below(2 * p.num_states);
    p.num_sources = 1 + rng.below(3);
    p.extra_physical_edges = rng.below(3);
    const auto inst = zero_arc_instance<W>(p, seed);
    INFO("seed " << seed);
    const auto& g = inst.seeded.hierarchy.graph;
    for (const auto& a : g.arcs) zero_nonloop += (a.length == W{0} && a.tail != a.head) ? 1 : 0;
    CHECK_FALSE(zero_cycle(g));
    const auto rep = sssp::paper::validate_seeded_ordered_hierarchy(inst.seeded, inst.physical);
    INFO(rep.summary());
    REQUIRE(rep.ok());
    CHECK_NOTHROW(compute_eval_reference<W>(inst.seeded, inst.physical));
  }
  CHECK(zero_nonloop > 200);
}

TEMPLATE_TEST_CASE("seeded instances: base-case catalogue serves every reachable vertex",
                   "[seeded-instances]", std::uint64_t, double) {
  using W = TestType;
  for (std::uint64_t seed : {1u, 2u, 3u}) {
    for (const auto& inst : base_case_catalogue<W>(seed)) {
      INFO(inst.name << " seed " << seed);
      const auto rep = sssp::paper::validate_seeded_ordered_hierarchy(inst.seeded, inst.physical);
      INFO(rep.summary());
      REQUIRE(rep.ok());
      const auto ref = compute_eval_reference<W>(inst.seeded, inst.physical);
      const auto served =
          served_targets<W>(inst.num_physical, inst.seeded.hierarchy.projection, ref);
      // Base case: states = vertices with both orientations of every edge, so a state-graph
      // path exists along every physical path and every reachable vertex is protected.
      std::vector<std::pair<StateId, W>> em;
      for (StateId x = 0; x < ref.state_distance.size(); ++x) {
        if (ref.state_distance[x] != sssp::WeightTraits<W>::infinity()) {
          em.emplace_back(x, ref.state_distance[x]);
        }
      }
      const auto minima =
          physical_minima<W>(inst.num_physical, inst.seeded.hierarchy.projection, em);
      for (sssp::Vertex v = 0; v < inst.num_physical; ++v) {
        INFO("vertex " << v);
        CHECK(served[v] == (ref.target[v] != sssp::WeightTraits<W>::infinity()));
        CHECK(minima[v] == ref.target[v]);
      }
      if (inst.name.find("/multi") != std::string::npos) CHECK(inst.seeded.sources.size() == 3);
    }
  }
}

TEST_CASE("eval invariants: Dijkstra traces on generated instances raise no violation",
          "[eval-invariants]") {
  // The flat state-graph Dijkstra emits dist_D, which is >= delta everywhere and = delta at
  // protected states: I1 and I2 must hold on every generated instance.
  using W = std::uint64_t;
  std::vector<SeededInstance<W>> all = reference_catalogue<W>(1000, 150);
  for (auto& x : base_case_catalogue<W>(11)) all.push_back(std::move(x));
  for (std::uint64_t seed = 1; seed <= 50; ++seed) {
    ReferenceParams p;
    p.zero_arcs = 8;
    p.num_sources = 1 + seed % 3;
    all.push_back(zero_arc_instance<W>(p, seed));
  }
  CheckerOptions opt;
  opt.allow_unscoped_emits = true;
  for (const auto& inst : all) {
    INFO(inst.name << " seed " << inst.seed);
    const auto ref = compute_eval_reference<W>(inst.seeded, inst.physical);
    InvariantChecker<W> c(inst.seeded, ref, opt);
    replay<W>(state_dijkstra_trace<W>(inst.seeded, ref), c);
    INFO(c.report().summary());
    CHECK(c.report().ok());
  }
}
