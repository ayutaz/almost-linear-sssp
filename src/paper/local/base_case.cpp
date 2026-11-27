// Base case of Algorithm 4 (l.2-5), Sec 2 l.286-290. See base_case.hpp.

#include "paper/local/base_case.hpp"

#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "paper/primitives/edge_pairs.hpp"
#include "paper/primitives/merge_events.hpp"

namespace sssp::paper {

template <class W>
LocalHierarchyResult<W> local_hierarchy_base_case(const CsrGraph<W>& g, Diagnostics* diag) {
  const std::size_t n = g.num_vertices();
  const auto edges = g.edges();
  // Alg 4 input contract: positive weights (zero components contracted by the caller).
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const W w = edges[i].weight;
    bool ok = w > W{0};
    if constexpr (std::is_integral_v<W>) ok = ok && w != WeightTraits<W>::infinity();
    if (!ok) {
      throw std::invalid_argument("local_hierarchy_base_case: edge " + std::to_string(i) +
                                  " does not have a positive finite weight (contract "
                                  "zero-weight components before Alg 4)");
    }
  }

  // Alg 4 l.3: one native entry state per vertex; both orientations of every nonloop edge.
  LocalHierarchyResult<W> out;
  StateGraph<W> d;
  d.num_states = n;
  std::vector<ArcId> arc_of_edge(edges.size(), kNoArc);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u == e.v) continue;
    const auto id = static_cast<EdgeId>(i);
    arc_of_edge[i] = static_cast<ArcId>(d.arcs.size());
    d.arcs.push_back({e.u, e.v, e.weight});
    d.arcs.push_back({e.v, e.u, e.weight});
    out.arc_edge.push_back(id);
    out.arc_edge.push_back(id);
  }
  std::vector<Vertex> projection(n);
  std::iota(projection.begin(), projection.end(), Vertex{0});
  out.native_entry.resize(n);
  std::iota(out.native_entry.begin(), out.native_entry.end(), StateId{0});
  out.hierarchy = make_leaf_forest(std::move(d), std::move(projection));
  auto& h = out.hierarchy;
  h.c0 = 1;  // links are MST edges of weight exactly h(C)

  // Alg 4 l.4, first half: retain a lightest edge per vertex pair (edge_pairs, two stable
  // counting sorts, Lemma 3.1). [Supplement: the tie-break by id is ours; any lightest edge is
  // allowed.]
  const std::vector<EdgeId> kept = lightest_edge_per_pair<W>(n, edges);

  // Alg 4 l.4, second half (Sec 2 l.288-290): Kruskal by (weight, id) with the MST edges of one
  // weight w merged simultaneously (merge_events, Sec 3.3 / Lemma 3.6): every component formed
  // at w that joins two or more previous components is one node of height w whose links are its
  // MST edges of weight w. Merge-forest node ids equal hierarchy node ids (leaves = states).
  const MergeForest<W> forest = build_merge_events<W>(n, edges, kept);
  for (const MergeNode<W>& node : forest.internal) {
    std::vector<PhysicalLink<W>> links;
    links.reserve(node.edges.size());
    for (const EdgeId e : node.edges) {
      links.push_back({edges[e].u, edges[e].v, node.height, arc_of_edge[e]});
    }
    add_internal_node(h, node.height, node.children, std::move(links));
  }
  collect_roots(h);

  if (diag != nullptr) {
    SSSP_COUNT(++diag->counters["local_base_case_calls"]);
    SSSP_COUNT(diag->counters["local_base_case_states"] += n);
    SSSP_COUNT(diag->counters["local_base_case_merge_nodes"] += forest.internal.size());
  }
  return out;
}

template LocalHierarchyResult<std::uint64_t> local_hierarchy_base_case<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, Diagnostics*);
template LocalHierarchyResult<double> local_hierarchy_base_case<double>(const CsrGraph<double>&,
                                                                        Diagnostics*);

}  // namespace sssp::paper
