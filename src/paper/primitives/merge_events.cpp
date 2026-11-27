// Kruskal merge forest with equal weights merged simultaneously (Sec 2 l.288-290, Alg 4 l.4,
// Sec 5.4 l.1198-1199, Sec 3.3). See merge_events.hpp.
#include "paper/primitives/merge_events.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include "paper/primitives/dsu.hpp"

namespace sssp::paper {

namespace {

constexpr std::size_t kNoSlot = std::numeric_limits<std::size_t>::max();

}  // namespace

template <class W>
MergeForest<W> build_merge_events(std::size_t n, std::span<const Edge<W>> edges,
                                  std::span<const EdgeId> candidates, Diagnostics* diag) {
  validate_edges<W>(n, edges);
  std::vector<EdgeId> order;
  order.reserve(candidates.size());
  for (const EdgeId e : candidates) {
    if (e >= edges.size()) {
      throw std::out_of_range("build_merge_events: candidate " + std::to_string(e) +
                              " is not an edge");
    }
    if (edges[e].u != edges[e].v) order.push_back(e);  // a loop never joins two components
  }
  // Order the merge events by (weight, id); Sec 3.3 charges this to the caller of the disjoint
  // sets, Sec 5.4 l.1199 sorts them ("sort the edges of its own MST").
  std::sort(order.begin(), order.end(), [&](EdgeId a, EdgeId b) {
    if (edges[a].weight < edges[b].weight) return true;
    if (edges[b].weight < edges[a].weight) return false;
    return a < b;
  });

  MergeForest<W> out;
  out.num_leaves = n;
  out.parent.assign(n, kNoMergeNode);
  DisjointSets dsu(n);                  // Sec 3.3, Lemma 3.6
  std::vector<MergeNodeId> node_of(n);  // indexed by a DSU root: its current forest node
  std::iota(node_of.begin(), node_of.end(), MergeNodeId{0});
  std::vector<std::size_t> slot(n, kNoSlot);  // DSU root -> pending node of this class
  std::vector<std::pair<DisjointSets::Item, DisjointSets::Item>> start_root;
  std::vector<std::size_t> tree;  // positions of this class's forest edges
  struct Pending {
    DisjointSets::Item root;
    std::vector<MergeNodeId> children;
    std::vector<EdgeId> edges;
  };
  std::vector<Pending> pending;
  std::uint64_t classes = 0;

  for (std::size_t begin = 0; begin < order.size();) {
    // One weight class [begin, end): all candidates of weight w (comparison model).
    const W w = edges[order[begin]].weight;
    std::size_t end = begin + 1;
    while (end < order.size() && !(w < edges[order[end]].weight)) ++end;
    ++classes;

    // The components before the class (all lighter forest edges merged) ...
    start_root.clear();
    for (std::size_t i = begin; i < end; ++i) {
      start_root.emplace_back(dsu.find(edges[order[i]].u), dsu.find(edges[order[i]].v));
    }
    // ... then the class's forest edges by Kruskal in (weight, id) order.
    tree.clear();
    for (std::size_t i = begin; i < end; ++i) {
      const auto [ru, rv] = start_root[i - begin];
      if (dsu.unite(ru, rv)) tree.push_back(i);
    }
    // Simultaneous merge: one node per component after the class that absorbed a forest edge;
    // its children are the nodes of the components before the class that it contains.
    pending.clear();
    for (const std::size_t i : tree) {
      const EdgeId e = order[i];
      const DisjointSets::Item r = dsu.find(edges[e].u);
      if (slot[r] == kNoSlot) {
        slot[r] = pending.size();
        pending.push_back(Pending{r, {}, {}});
      }
      Pending& p = pending[slot[r]];
      p.children.push_back(node_of[start_root[i - begin].first]);
      p.children.push_back(node_of[start_root[i - begin].second]);
      p.edges.push_back(e);
      out.forest_edges.push_back(e);
    }
    for (Pending& p : pending) {
      std::sort(p.children.begin(), p.children.end());
      p.children.erase(std::unique(p.children.begin(), p.children.end()), p.children.end());
      const auto id = static_cast<MergeNodeId>(out.num_nodes());
      for (const MergeNodeId c : p.children) out.parent[c] = id;
      out.parent.push_back(kNoMergeNode);
      out.internal.push_back(MergeNode<W>{w, std::move(p.children), std::move(p.edges)});
      node_of[p.root] = id;
      slot[p.root] = kNoSlot;
    }
    begin = end;
  }
  for (std::size_t x = 0; x < out.parent.size(); ++x) {
    if (out.parent[x] == kNoMergeNode) out.roots.push_back(static_cast<MergeNodeId>(x));
  }

  static_cast<void>(classes);  // only read when SSSP_DIAGNOSTICS is on
  if (diag != nullptr) {
    SSSP_COUNT(++diag->counters["merge_events.calls"]);
    SSSP_COUNT(diag->counters["merge_events.nodes"] += out.internal.size());
    SSSP_COUNT(diag->counters["merge_events.forest_edges"] += out.forest_edges.size());
    SSSP_COUNT(diag->counters["merge_events.weight_classes"] += classes);
    SSSP_COUNT(dsu.add_counters(*diag));
  }
  return out;
}

template <class W>
MergeForest<W> build_merge_events(std::size_t n, std::span<const Edge<W>> edges,
                                  Diagnostics* diag) {
  std::vector<EdgeId> all(edges.size());
  std::iota(all.begin(), all.end(), EdgeId{0});
  return build_merge_events<W>(n, edges, all, diag);
}

template MergeForest<std::uint64_t> build_merge_events<std::uint64_t>(
    std::size_t, std::span<const Edge<std::uint64_t>>, std::span<const EdgeId>, Diagnostics*);
template MergeForest<double> build_merge_events<double>(std::size_t, std::span<const Edge<double>>,
                                                        std::span<const EdgeId>, Diagnostics*);
template MergeForest<std::uint64_t> build_merge_events<std::uint64_t>(
    std::size_t, std::span<const Edge<std::uint64_t>>, Diagnostics*);
template MergeForest<double> build_merge_events<double>(std::size_t, std::span<const Edge<double>>,
                                                        Diagnostics*);

}  // namespace sssp::paper
