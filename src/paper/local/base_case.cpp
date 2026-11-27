// Base case of Algorithm 4 (l.2-5), Sec 2 l.286-290. See base_case.hpp.

#include "paper/local/base_case.hpp"

#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace sssp::paper {

namespace {

constexpr std::size_t kNpos = std::numeric_limits<std::size_t>::max();

struct Dsu {
  std::vector<std::size_t> p;
  explicit Dsu(std::size_t n) : p(n) { std::iota(p.begin(), p.end(), std::size_t{0}); }
  std::size_t find(std::size_t x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  }
  bool unite(std::size_t a, std::size_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;
    p[a] = b;
    return true;
  }
};

// Stable counting sort of edge ids by an integer key in [0, n) (Lemma 3.1).
template <class Key>
std::vector<EdgeId> counting_sort(const std::vector<EdgeId>& ids, std::size_t n, Key key) {
  std::vector<std::size_t> start(n + 1, 0);
  for (const EdgeId e : ids) ++start[key(e) + 1];
  for (std::size_t k = 0; k < n; ++k) start[k + 1] += start[k];
  std::vector<EdgeId> out(ids.size());
  for (const EdgeId e : ids) out[start[key(e)]++] = e;
  return out;
}

}  // namespace

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
  std::vector<EdgeId> nonloop;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u == e.v) continue;
    const auto id = static_cast<EdgeId>(i);
    arc_of_edge[i] = static_cast<ArcId>(d.arcs.size());
    d.arcs.push_back({e.u, e.v, e.weight});
    d.arcs.push_back({e.v, e.u, e.weight});
    out.arc_edge.push_back(id);
    out.arc_edge.push_back(id);
    nonloop.push_back(id);
  }
  std::vector<Vertex> projection(n);
  std::iota(projection.begin(), projection.end(), Vertex{0});
  out.native_entry.resize(n);
  std::iota(out.native_entry.begin(), out.native_entry.end(), StateId{0});
  out.hierarchy = make_leaf_forest(std::move(d), std::move(projection));
  auto& h = out.hierarchy;
  h.c0 = 1;  // links are MST edges of weight exactly h(C)

  // Alg 4 l.4, first half: retain a lightest edge per vertex pair. Two stable counting sorts
  // order the edges by (min endpoint, max endpoint, id); the lightest (weight, id) per pair
  // survives. [Supplement: the tie-break by id is ours; any lightest edge is allowed.]
  const auto lo = [&](EdgeId e) -> std::size_t { return std::min(edges[e].u, edges[e].v); };
  const auto hi = [&](EdgeId e) -> std::size_t { return std::max(edges[e].u, edges[e].v); };
  const std::vector<EdgeId> by_pair = counting_sort(counting_sort(nonloop, n, hi), n, lo);
  std::vector<EdgeId> kept;
  for (std::size_t i = 0; i < by_pair.size();) {
    EdgeId best = by_pair[i];
    std::size_t j = i + 1;
    while (j < by_pair.size() && lo(by_pair[j]) == lo(best) && hi(by_pair[j]) == hi(best)) {
      if (edges[by_pair[j]].weight < edges[best].weight) best = by_pair[j];
      ++j;
    }
    kept.push_back(best);
    i = j;
  }

  // Alg 4 l.4, second half: Kruskal by (weight, id); the MST edges of one weight w are merged
  // simultaneously (Lemma 3.5 / Lemma 3.6 merge events): every resulting component that joins
  // two or more previous components becomes one node of height w whose links are its MST
  // edges of weight w.
  std::sort(kept.begin(), kept.end(), [&](EdgeId a, EdgeId b) {
    if (edges[a].weight < edges[b].weight) return true;
    if (edges[b].weight < edges[a].weight) return false;
    return a < b;
  });
  Dsu dsu(n);
  std::vector<NodeId> node_of(n);  // indexed by DSU representative
  std::iota(node_of.begin(), node_of.end(), NodeId{0});
  std::vector<std::size_t> slot(n, kNpos);
  struct Bucket {
    std::size_t root;
    std::vector<std::size_t> reps;
    std::vector<PhysicalLink<W>> links;
  };
  std::uint64_t merge_nodes = 0;
  for (std::size_t begin = 0; begin < kept.size();) {
    const W w = edges[kept[begin]].weight;
    std::size_t end = begin;
    while (end < kept.size() && edges[kept[end]].weight == w) ++end;
    // Components at the start of the weight class (all MST edges of weight < w merged).
    std::vector<std::pair<std::size_t, std::size_t>> start_rep;
    start_rep.reserve(end - begin);
    for (std::size_t i = begin; i < end; ++i) {
      start_rep.emplace_back(dsu.find(edges[kept[i]].u), dsu.find(edges[kept[i]].v));
    }
    std::vector<std::size_t> tree;  // positions in [begin, end) that are MST edges
    for (std::size_t i = begin; i < end; ++i) {
      const auto [ru, rv] = start_rep[i - begin];
      if (dsu.unite(ru, rv)) tree.push_back(i);
    }
    std::vector<Bucket> buckets;
    for (const std::size_t i : tree) {
      const EdgeId e = kept[i];
      const std::size_t r = dsu.find(edges[e].u);
      if (slot[r] == kNpos) {
        slot[r] = buckets.size();
        buckets.push_back(Bucket{r, {}, {}});
      }
      auto& bk = buckets[slot[r]];
      bk.reps.push_back(start_rep[i - begin].first);
      bk.reps.push_back(start_rep[i - begin].second);
      bk.links.push_back({edges[e].u, edges[e].v, w, arc_of_edge[e]});
    }
    for (auto& bk : buckets) {
      std::sort(bk.reps.begin(), bk.reps.end());
      bk.reps.erase(std::unique(bk.reps.begin(), bk.reps.end()), bk.reps.end());
      std::vector<NodeId> children;
      children.reserve(bk.reps.size());
      for (const std::size_t rep : bk.reps) children.push_back(node_of[rep]);
      node_of[bk.root] = add_internal_node(h, w, std::move(children), std::move(bk.links));
      slot[bk.root] = kNpos;
      ++merge_nodes;
    }
    begin = end;
  }
  collect_roots(h);

  static_cast<void>(merge_nodes);  // only read when SSSP_DIAGNOSTICS is on
  if (diag != nullptr) {
    SSSP_COUNT(++diag->counters["local_base_case_calls"]);
    SSSP_COUNT(diag->counters["local_base_case_states"] += n);
    SSSP_COUNT(diag->counters["local_base_case_merge_nodes"] += merge_nodes);
  }
  return out;
}

template LocalHierarchyResult<std::uint64_t> local_hierarchy_base_case<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, Diagnostics*);
template LocalHierarchyResult<double> local_hierarchy_base_case<double>(const CsrGraph<double>&,
                                                                        Diagnostics*);

}  // namespace sssp::paper
