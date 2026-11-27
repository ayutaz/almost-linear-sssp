// [Supplement, not in the paper] Threshold-SCC reference hierarchy; see the header and
// research notes 02 sec 2.11. Test support for Sec 4 (Def 4.1, l.646-658).

#include "paper/hierarchy/reference_hierarchy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace sssp::paper {

namespace {

constexpr std::size_t kNpos = std::numeric_limits<std::size_t>::max();

template <class W>
bool positive_finite(W w) noexcept {
  if constexpr (std::is_floating_point_v<W>) {
    return std::isfinite(w) && w > W{0};
  }
  return w > W{0} && w != WeightTraits<W>::infinity();
}

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

// Iterative Tarjan SCC on a CSR digraph with k vertices. Returns the SCC id of every vertex.
std::vector<std::size_t> tarjan_scc(std::size_t k, const std::vector<std::size_t>& off,
                                    const std::vector<std::size_t>& adj) {
  std::vector<std::size_t> index(k, kNpos);
  std::vector<std::size_t> low(k, 0);
  std::vector<std::size_t> comp(k, kNpos);
  std::vector<std::size_t> cursor(k, 0);
  std::vector<std::uint8_t> on_stack(k, 0);
  std::vector<std::size_t> stack;
  std::vector<std::size_t> call;
  std::size_t next = 0;
  std::size_t count = 0;
  auto open = [&](std::size_t v) {
    index[v] = low[v] = next++;
    cursor[v] = off[v];
    stack.push_back(v);
    on_stack[v] = 1;
    call.push_back(v);
  };
  for (std::size_t s = 0; s < k; ++s) {
    if (index[s] != kNpos) continue;
    open(s);
    while (!call.empty()) {
      const std::size_t v = call.back();
      if (cursor[v] < off[v + 1]) {
        const std::size_t w = adj[cursor[v]++];
        if (index[w] == kNpos) {
          open(w);
        } else if (on_stack[w] != 0) {
          low[v] = std::min(low[v], index[w]);
        }
        continue;
      }
      call.pop_back();
      if (low[v] == index[v]) {
        std::size_t w = kNpos;
        do {
          w = stack.back();
          stack.pop_back();
          on_stack[w] = 0;
          comp[w] = count;
        } while (w != v);
        ++count;
      }
      if (!call.empty()) low[call.back()] = std::min(low[call.back()], low[v]);
    }
  }
  return comp;
}

}  // namespace

template <class W>
OrderedHierarchy<W> build_threshold_scc_hierarchy(StateGraph<W> graph,
                                                  std::vector<Vertex> projection) {
  const std::size_t n = graph.num_states;
  if (projection.size() != n) {
    throw std::invalid_argument("build_threshold_scc_hierarchy: projection size differs");
  }
  std::vector<ArcId> order;
  for (std::size_t i = 0; i < graph.arcs.size(); ++i) {
    const auto& a = graph.arcs[i];
    if (a.tail >= n || a.head >= n) {
      throw std::out_of_range("build_threshold_scc_hierarchy: arc " + std::to_string(i) +
                              " has an endpoint that is not a state");
    }
    if (!positive_finite(a.length)) {
      throw std::invalid_argument("build_threshold_scc_hierarchy: arc " + std::to_string(i) +
                                  " length is not positive and finite (contract zero-length "
                                  "components first)");
    }
    if (a.tail != a.head) order.push_back(static_cast<ArcId>(i));
  }
  std::stable_sort(order.begin(), order.end(),
                   [&](ArcId x, ArcId y) { return graph.arcs[x].length < graph.arcs[y].length; });

  OrderedHierarchy<W> h = make_leaf_forest(std::move(graph), std::move(projection));
  h.c0 = 1;
  const auto& arcs = h.graph.arcs;
  Dsu comp(n);
  std::vector<NodeId> node_of(n);  // indexed by DSU representative
  std::iota(node_of.begin(), node_of.end(), NodeId{0});

  std::vector<std::size_t> local(n, kNpos);
  for (std::size_t begin = 0; begin < order.size();) {
    const W w = arcs[order[begin]].length;
    std::size_t end = begin;
    while (end < order.size() && arcs[order[end]].length == w) ++end;

    // Contract the current components (SCCs of arcs < w) and take arcs <= w between them.
    std::vector<std::size_t> reps;
    std::fill(local.begin(), local.end(), kNpos);
    for (std::size_t x = 0; x < n; ++x) {
      const std::size_t r = comp.find(x);
      if (local[r] == kNpos) {
        local[r] = reps.size();
        reps.push_back(r);
      }
    }
    const std::size_t k = reps.size();
    auto local_of = [&](StateId x) { return local[comp.find(x)]; };
    std::vector<std::size_t> off(k + 1, 0);
    for (std::size_t i = 0; i < end; ++i) {
      const auto& a = arcs[order[i]];
      if (local_of(a.tail) != local_of(a.head)) ++off[local_of(a.tail) + 1];
    }
    for (std::size_t v = 0; v < k; ++v) off[v + 1] += off[v];
    std::vector<std::size_t> adj(off[k]);
    std::vector<std::size_t> fill(off.begin(), off.end() - 1);
    for (std::size_t i = 0; i < end; ++i) {
      const auto& a = arcs[order[i]];
      const std::size_t cu = local_of(a.tail);
      const std::size_t cv = local_of(a.head);
      if (cu != cv) adj[fill[cu]++] = cv;
    }
    const std::vector<std::size_t> scc = tarjan_scc(k, off, adj);

    // Group the contracted vertices by SCC; every SCC with >= 2 members is a new node.
    std::size_t num_scc = 0;
    for (const std::size_t c : scc) num_scc = std::max(num_scc, c + 1);
    std::vector<std::vector<std::size_t>> members(num_scc);
    for (std::size_t v = 0; v < k; ++v) members[scc[v]].push_back(v);
    std::vector<std::vector<PhysicalLink<W>>> links(num_scc);
    Dsu tree(k);
    for (std::size_t i = 0; i < end; ++i) {
      const auto& a = arcs[order[i]];
      const std::size_t cu = local_of(a.tail);
      const std::size_t cv = local_of(a.head);
      if (cu == cv || scc[cu] != scc[cv]) continue;
      if (tree.unite(cu, cv)) links[scc[cu]].push_back({a.tail, a.head, a.length, order[i]});
    }
    for (std::size_t c = 0; c < num_scc; ++c) {
      if (members[c].size() < 2) continue;
      std::vector<NodeId> children;
      children.reserve(members[c].size());
      for (const std::size_t v : members[c]) children.push_back(node_of[reps[v]]);
      const NodeId id = add_internal_node(h, w, std::move(children), std::move(links[c]));
      for (const std::size_t v : members[c]) comp.unite(reps[v], reps[members[c][0]]);
      node_of[comp.find(reps[members[c][0]])] = id;
    }
    begin = end;
  }
  collect_roots(h);
  return h;
}

template <class W>
CsrGraph<W> physical_graph_of_arcs(const StateGraph<W>& graph,
                                   const std::vector<Vertex>& projection,
                                   std::size_t num_physical) {
  if (projection.size() != graph.num_states) {
    throw std::invalid_argument("physical_graph_of_arcs: projection size differs");
  }
  std::vector<Edge<W>> edges;
  edges.reserve(graph.arcs.size());
  for (const auto& a : graph.arcs) {
    if (a.tail >= graph.num_states || a.head >= graph.num_states) {
      throw std::out_of_range("physical_graph_of_arcs: arc endpoint is not a state");
    }
    edges.push_back({projection[a.tail], projection[a.head], a.length});
  }
  return CsrGraph<W>::from_edges(num_physical, std::move(edges));
}

template OrderedHierarchy<std::uint64_t> build_threshold_scc_hierarchy<std::uint64_t>(
    StateGraph<std::uint64_t>, std::vector<Vertex>);
template OrderedHierarchy<double> build_threshold_scc_hierarchy<double>(StateGraph<double>,
                                                                        std::vector<Vertex>);
template CsrGraph<std::uint64_t> physical_graph_of_arcs<std::uint64_t>(
    const StateGraph<std::uint64_t>&, const std::vector<Vertex>&, std::size_t);
template CsrGraph<double> physical_graph_of_arcs<double>(const StateGraph<double>&,
                                                         const std::vector<Vertex>&, std::size_t);

}  // namespace sssp::paper
