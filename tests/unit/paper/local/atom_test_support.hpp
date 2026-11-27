#pragma once
// Test support for Sec 5.5 (replacing atoms by capped local hierarchies): random partition
// hierarchies (Def 5.12) built from the definitions only, together with the fine ordered
// hierarchy they were collapsed from, so that Evaluate can be run before and after the
// replacement on the same underlying state graph. Nothing here calls
// src/paper/local/atom_replacement.
//
// Construction (every random choice from SplitMix64 with a fixed seed):
//   1. The states 0..n-1 are cut into consecutive groups of random sizes. Inside a group: a
//      random spanning tree plus chords of UNDIRECTED edges (both orientations become state
//      arcs). Between groups: a random tree of undirected edges (if `connect`), extra undirected
//      edges and `one_way` one-way arcs.
//   2. The fine hierarchy is the threshold-SCC reference hierarchy (c0 = 1) of all arcs. The
//      physical graph is the arc multigraph under the projection (identity, or folded onto
//      fewer physical vertices).
//   3. Atoms are disjoint internal nodes X of the fine hierarchy with no one-way arc inside X.
//      X is an SCC of the arcs of length <= h(X), and every arc inside X is an orientation of an
//      undirected edge, so the undirected edges of length <= h(X) inside X connect it: their
//      Kruskal tree has every edge <= h(X) < 2c for every cutoff c > h(X) / 2 (Def 5.11). The
//      cutoff is drawn in (h(X)/2, h(parent)] (Def 5.11: the parent's height is at least c), or in
//      (h(X)/2, 2 h(X)] for a root. Optionally, single states become one-state atoms with a cutoff
//      in (0, h(parent)].
//   4. The partition hierarchy keeps every fine node outside the atoms' subtrees; each X becomes
//      an atom leaf. The undirected edges inside an atom form its internal graph H; every other
//      arc stays in D (renumbered, and link backing arcs renumbered with it). The fine hierarchy
//      satisfies Def 4.1, and collapsing a subtree into one part changes nothing above it, so the
//      partition hierarchy satisfies Def 4.1 on its parts (Def 5.12).
//
// Weights: std::uint64_t, or quarter-integers in double (every sum stays below 2^53 and is
// exact).

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "paper/local/atom_replacement.hpp"
#include "sssp/graph.hpp"
#include "support/prng.hpp"

namespace sssp::test::atom {

using sssp::paper::ArcId;
using sssp::paper::NodeId;
using sssp::paper::OrderedHierarchy;
using sssp::paper::PartitionHierarchy;
using sssp::paper::PhysicalSource;
using sssp::paper::ReplaceableAtom;
using sssp::paper::StateArc;
using sssp::paper::StateGraph;
using sssp::paper::StateId;
using sssp::paper::StateOffer;

struct PartitionParams {
  std::size_t num_states = 30;
  std::size_t min_group = 1;
  std::size_t max_group = 8;
  std::uint64_t internal_min = 1;  ///< internal edge lengths in [internal_min, internal_max]
  std::uint64_t internal_max = 10;
  bool internal_exponential = false;  ///< internal lengths 2^U[0, 40] instead
  std::size_t chords_per_group = 1;
  std::uint64_t outer_min = 1;  ///< lengths between groups in [outer_min, outer_max]
  std::uint64_t outer_max = 20;
  bool connect = true;                   ///< join the groups by a random tree of undirected edges
  std::size_t cross_edges = 6;           ///< extra undirected edges between groups
  std::size_t one_way = 0;               ///< one-way arcs between groups
  bool fold_projection = false;          ///< project onto about 2n/3 physical vertices
  std::uint64_t atom_percent = 70;       ///< chance that a candidate node becomes an atom
  std::uint64_t singleton_percent = 10;  ///< chance that a remaining leaf becomes an atom
  bool largest_first = false;            ///< try candidates by decreasing size
};

template <class W>
struct PartitionInstance {
  std::size_t num_physical = 0;
  CsrGraph<W> physical;
  StateGraph<W> full_graph;  ///< D plus both orientations of every internal edge
  std::vector<Vertex> projection;
  OrderedHierarchy<W> fine;
  PartitionHierarchy<W> partition;
  bool all_undirected = true;  ///< no one-way arc anywhere
};

template <class W>
W length_of(std::uint64_t x, SplitMix64& rng) {
  if constexpr (std::is_floating_point_v<W>) {
    return static_cast<W>(x) + static_cast<W>(rng.between(0, 3)) * W{0.25};
  } else {
    return static_cast<W>(x);
  }
}

/// A cutoff c with lo_excl < c <= hi (both positive, lo_excl < hi). Integers: c in
/// [floor(lo_excl) + 1, hi]; doubles: a multiple of 1/8 (exact).
template <class W>
W draw_cutoff(W lo_excl, W hi, SplitMix64& rng) {
  if constexpr (std::is_floating_point_v<W>) {
    const auto lo_u = static_cast<std::uint64_t>(lo_excl * 8) + 1;
    const auto hi_u = static_cast<std::uint64_t>(hi * 8);
    return static_cast<W>(rng.between(lo_u, std::max(lo_u, hi_u))) / 8;
  } else {
    return static_cast<W>(rng.between(lo_excl + 1, std::max<W>(lo_excl + 1, hi)));
  }
}

template <class W>
W half_of(W h) {
  if constexpr (std::is_floating_point_v<W>) {
    return h / 2;
  } else {
    return h / 2;  // floor; c > floor(h/2) and c integer give 2c > h
  }
}

/// Minimal union-find for the test (independent of src/paper/primitives/dsu.hpp).
class TestDsu {
 public:
  explicit TestDsu(std::size_t n) : p_(n) { std::iota(p_.begin(), p_.end(), std::size_t{0}); }
  std::size_t find(std::size_t x) {
    while (p_[x] != x) x = p_[x] = p_[p_[x]];
    return x;
  }
  bool unite(std::size_t a, std::size_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;
    p_[b] = a;
    return true;
  }

 private:
  std::vector<std::size_t> p_;
};

/// Kruskal spanning tree (by length, then id) of a connected graph on k vertices.
template <class W>
std::vector<EdgeId> kruskal_tree(std::size_t k, const std::vector<Edge<W>>& edges) {
  std::vector<EdgeId> ids(edges.size());
  std::iota(ids.begin(), ids.end(), EdgeId{0});
  std::stable_sort(ids.begin(), ids.end(),
                   [&](EdgeId a, EdgeId b) { return edges[a].weight < edges[b].weight; });
  TestDsu dsu(k);
  std::vector<EdgeId> tree;
  for (const EdgeId id : ids) {
    if (dsu.unite(edges[id].u, edges[id].v)) tree.push_back(id);
  }
  if (k > 0 && tree.size() + 1 != k) throw std::logic_error("kruskal_tree: disconnected");
  return tree;
}

/// States below every node of `h` (DFS per node; fine for test sizes).
template <class W>
std::vector<std::vector<StateId>> states_below(const OrderedHierarchy<W>& h) {
  std::vector<std::vector<StateId>> out(h.nodes.size());
  for (NodeId v = 0; v < h.nodes.size(); ++v) {
    std::vector<NodeId> stack{v};
    while (!stack.empty()) {
      const NodeId u = stack.back();
      stack.pop_back();
      if (h.nodes[u].children.empty()) {
        out[v].push_back(h.nodes[u].state);
      } else {
        for (const NodeId c : h.nodes[u].children) stack.push_back(c);
      }
    }
    std::sort(out[v].begin(), out[v].end());
  }
  return out;
}

template <class W>
PartitionInstance<W> make_partition_instance(const PartitionParams& p, std::uint64_t seed) {
  SplitMix64 rng(seed);
  const std::size_t n = p.num_states;
  PartitionInstance<W> inst;

  // 1. Groups, undirected edges and one-way arcs.
  std::vector<std::size_t> start{0};
  while (start.back() < n) {
    start.push_back(
        std::min(n, start.back() + p.min_group +
                        static_cast<std::size_t>(rng.below(p.max_group - p.min_group + 1))));
  }
  const std::size_t groups = start.size() - 1;
  const auto internal_len = [&] {
    if (p.internal_exponential) return length_of<W>(std::uint64_t{1} << rng.between(0, 40), rng);
    return length_of<W>(rng.between(p.internal_min, p.internal_max), rng);
  };
  const auto outer_len = [&] { return length_of<W>(rng.between(p.outer_min, p.outer_max), rng); };
  const auto pick = [&](std::size_t g) {
    return static_cast<Vertex>(start[g] + rng.below(start[g + 1] - start[g]));
  };
  std::vector<Edge<W>> undirected;
  std::vector<StateArc<W>> one_way;
  for (std::size_t g = 0; g < groups; ++g) {
    const std::size_t size = start[g + 1] - start[g];
    for (std::size_t v = start[g] + 1; v < start[g + 1]; ++v) {
      undirected.push_back({static_cast<Vertex>(start[g] + rng.below(v - start[g])),
                            static_cast<Vertex>(v), internal_len()});
    }
    if (size >= 2) {
      for (std::size_t i = 0; i < p.chords_per_group; ++i) {
        const Vertex a = pick(g);
        Vertex b = pick(g);
        if (a == b) b = static_cast<Vertex>(a == start[g] ? a + 1 : a - 1);
        undirected.push_back({a, b, internal_len()});
      }
    }
  }
  if (groups >= 2) {
    if (p.connect) {
      for (std::size_t g = 1; g < groups; ++g) {
        undirected.push_back({pick(rng.below(g)), pick(g), outer_len()});
      }
    }
    for (std::size_t i = 0; i < p.cross_edges; ++i) {
      const std::size_t a = rng.below(groups);
      const std::size_t b = (a + 1 + rng.below(groups - 1)) % groups;
      undirected.push_back({pick(a), pick(b), outer_len()});
    }
    for (std::size_t i = 0; i < p.one_way; ++i) {
      const std::size_t a = rng.below(groups);
      const std::size_t b = (a + 1 + rng.below(groups - 1)) % groups;
      one_way.push_back({pick(a), pick(b), outer_len()});
    }
  }
  inst.all_undirected = one_way.empty();

  // Arcs: both orientations of every undirected edge, then the one-way arcs.
  inst.full_graph.num_states = n;
  std::vector<std::size_t> arc_edge;  // arc -> undirected edge id, or edges.size() if one-way
  for (std::size_t e = 0; e < undirected.size(); ++e) {
    inst.full_graph.arcs.push_back({undirected[e].u, undirected[e].v, undirected[e].weight});
    inst.full_graph.arcs.push_back({undirected[e].v, undirected[e].u, undirected[e].weight});
    arc_edge.push_back(e);
    arc_edge.push_back(e);
  }
  for (const auto& a : one_way) {
    inst.full_graph.arcs.push_back(a);
    arc_edge.push_back(undirected.size());
  }

  // 2. Projection, physical graph, fine hierarchy.
  inst.num_physical = p.fold_projection ? std::max<std::size_t>(1, (2 * n) / 3) : n;
  inst.projection.resize(n);
  for (std::size_t x = 0; x < n; ++x) {
    inst.projection[x] =
        static_cast<Vertex>(x < inst.num_physical ? x : rng.below(inst.num_physical));
  }
  inst.physical =
      sssp::paper::physical_graph_of_arcs(inst.full_graph, inst.projection, inst.num_physical);
  inst.fine = sssp::paper::build_threshold_scc_hierarchy(inst.full_graph, inst.projection);
  const auto& fine = inst.fine;
  const auto below = states_below(fine);

  // 3. Candidate atoms: internal nodes with no one-way arc inside.
  std::vector<std::uint8_t> has_one_way(fine.nodes.size(), 0);
  {
    std::vector<std::uint8_t> in(n, 0);
    for (NodeId v = 0; v < fine.nodes.size(); ++v) {
      if (fine.is_leaf(v)) continue;
      for (const StateId x : below[v]) in[x] = 1;
      for (std::size_t i = 0; i < inst.full_graph.arcs.size(); ++i) {
        const auto& a = inst.full_graph.arcs[i];
        if (arc_edge[i] == undirected.size() && in[a.tail] != 0 && in[a.head] != 0) {
          has_one_way[v] = 1;
        }
      }
      for (const StateId x : below[v]) in[x] = 0;
    }
  }
  std::vector<NodeId> candidates;
  for (NodeId v = 0; v < fine.nodes.size(); ++v) {
    if (!fine.is_leaf(v) && has_one_way[v] == 0) candidates.push_back(v);
  }
  for (std::size_t i = candidates.size(); i > 1; --i) {
    std::swap(candidates[i - 1], candidates[rng.below(i)]);
  }
  if (p.largest_first) {
    std::stable_sort(candidates.begin(), candidates.end(),
                     [&](NodeId a, NodeId b) { return below[a].size() > below[b].size(); });
  }
  std::vector<std::uint8_t> used(n, 0);  // state already inside a chosen atom
  std::vector<NodeId> chosen;
  for (const NodeId v : candidates) {
    if (rng.below(100) >= p.atom_percent) continue;
    bool free = true;
    for (const StateId x : below[v]) free = free && used[x] == 0;
    if (!free) continue;
    // An ancestor of a chosen node contains its states, so disjoint state sets mean neither is
    // an ancestor of the other.
    for (const StateId x : below[v]) used[x] = 1;
    chosen.push_back(v);
  }
  for (StateId x = 0; x < n; ++x) {
    if (used[x] == 0 && rng.below(100) < p.singleton_percent) {
      used[x] = 1;
      chosen.push_back(fine.leaf_of_state[x]);
    }
  }

  // 4. The partition hierarchy.
  auto& ph = inst.partition;
  ph.projection = inst.projection;
  ph.c0 = fine.c0;
  std::vector<std::uint32_t> atom_of_state(n, sssp::paper::kNoAtom);
  std::vector<std::uint32_t> atom_of_node(fine.nodes.size(), sssp::paper::kNoAtom);
  for (std::uint32_t a = 0; a < chosen.size(); ++a) {
    atom_of_node[chosen[a]] = a;
    for (const StateId x : below[chosen[a]]) atom_of_state[x] = a;
  }
  // Removed nodes: strict descendants of chosen nodes.
  std::vector<std::uint8_t> removed(fine.nodes.size(), 0);
  for (const NodeId v : chosen) {
    std::vector<NodeId> stack(fine.nodes[v].children.begin(), fine.nodes[v].children.end());
    while (!stack.empty()) {
      const NodeId u = stack.back();
      stack.pop_back();
      removed[u] = 1;
      for (const NodeId c : fine.nodes[u].children) stack.push_back(c);
    }
  }
  std::vector<NodeId> new_id(fine.nodes.size(), sssp::paper::kNoNode);
  for (NodeId v = 0; v < fine.nodes.size(); ++v) {
    if (removed[v] == 0) {
      new_id[v] = static_cast<NodeId>(ph.nodes.size());
      ph.nodes.emplace_back();
      ph.leaf_atom.push_back(atom_of_node[v]);
    }
  }
  // D: arcs not inside one atom.
  std::vector<ArcId> new_arc(inst.full_graph.arcs.size(), sssp::paper::kNoArc);
  ph.graph.num_states = n;
  for (std::size_t i = 0; i < inst.full_graph.arcs.size(); ++i) {
    const auto& a = inst.full_graph.arcs[i];
    if (atom_of_state[a.tail] != sssp::paper::kNoAtom &&
        atom_of_state[a.tail] == atom_of_state[a.head]) {
      continue;
    }
    new_arc[i] = static_cast<ArcId>(ph.graph.arcs.size());
    ph.graph.arcs.push_back(a);
  }
  for (NodeId v = 0; v < fine.nodes.size(); ++v) {
    if (removed[v] != 0) continue;
    auto& out = ph.nodes[new_id[v]];
    const auto& in = fine.nodes[v];
    out.parent = in.parent == sssp::paper::kNoNode ? sssp::paper::kNoNode : new_id[in.parent];
    if (atom_of_node[v] != sssp::paper::kNoAtom) continue;  // atom leaf: no state, no children
    out.state = in.state;
    out.height = in.height;
    for (const NodeId c : in.children) out.children.push_back(new_id[c]);
    for (auto link : in.links) {
      if (link.arc != sssp::paper::kNoArc) {
        if (new_arc[link.arc] == sssp::paper::kNoArc) {
          throw std::logic_error("make_partition_instance: link backed by an atom-internal arc");
        }
        link.arc = new_arc[link.arc];
      }
      out.links.push_back(link);
    }
  }
  // Atoms: states, internal graph, Kruskal tree, cutoff.
  for (std::uint32_t a = 0; a < chosen.size(); ++a) {
    const NodeId v = chosen[a];
    ReplaceableAtom<W> atom;
    atom.states = below[v];
    std::vector<Vertex> local(n, kNoVertex);
    for (std::size_t i = 0; i < atom.states.size(); ++i) {
      local[atom.states[i]] = static_cast<Vertex>(i);
    }
    for (const auto& e : undirected) {
      if (local[e.u] != kNoVertex && local[e.v] != kNoVertex) {
        atom.edges.push_back({local[e.u], local[e.v], e.weight});
      }
    }
    atom.spanning_tree = kruskal_tree<W>(atom.states.size(), atom.edges);
    const NodeId parent = fine.nodes[v].parent;
    if (fine.is_leaf(v)) {
      const W hi = parent == sssp::paper::kNoNode ? W{8} : fine.nodes[parent].height;
      atom.cutoff = draw_cutoff<W>(W{0}, hi, rng);
    } else {
      const W h = fine.nodes[v].height;
      const W hi = parent == sssp::paper::kNoNode ? h + h : fine.nodes[parent].height;
      atom.cutoff = draw_cutoff<W>(half_of(h), hi, rng);
    }
    for (const EdgeId id : atom.spanning_tree) {
      if (!(atom.edges[id].weight < atom.cutoff + atom.cutoff)) {
        throw std::logic_error("make_partition_instance: tree edge not below 2c");
      }
    }
    ph.atoms.push_back(std::move(atom));
  }
  return inst;
}

/// Sources and offers on the old states: `count` sources (v, offset) with offset 0 for the
/// first and in [0, 5] otherwise, each offered at state v (states 0..num_physical-1 have image
/// v = x, so the offer is feasible and tight).
template <class W>
std::pair<std::vector<PhysicalSource<W>>, std::vector<StateOffer<W>>> make_seeds(
    const PartitionInstance<W>& inst, std::size_t count, SplitMix64& rng) {
  std::vector<PhysicalSource<W>> sources;
  std::vector<StateOffer<W>> offers;
  for (std::size_t i = 0; i < count; ++i) {
    const auto v = static_cast<Vertex>(rng.below(inst.num_physical));
    const W off = i == 0 ? W{0} : static_cast<W>(rng.between(0, 5));
    sources.push_back({v, off});
    offers.push_back({static_cast<StateId>(v), off});
  }
  return {sources, offers};
}

}  // namespace sssp::test::atom
