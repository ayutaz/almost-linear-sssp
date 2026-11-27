// Owner decomposition, Sec 5.1 (Def 5.1, Lemmas 5.2-5.3). See owner_decomposition.hpp.

#include "paper/local/owner_decomposition.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "paper/model/params.hpp"
#include "paper/primitives/dsu.hpp"
#include "paper/primitives/mst.hpp"

namespace sssp::paper {

namespace {

/// b = L_N^4 in Sec 5.1 l.967 (the exponent of the separator, not of the size threshold).
constexpr unsigned kLocalSeparatorExponent = 4;

template <class W>
void require_positive_weights(const CsrGraph<W>& g) {
  const auto edges = g.edges();
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const W w = edges[i].weight;
    bool ok = W{0} < w;
    if constexpr (std::is_integral_v<W>) ok = ok && w != WeightTraits<W>::infinity();
    if (!ok) {
      throw std::invalid_argument("owner_decomposition: edge " + std::to_string(i) +
                                  " does not have a positive finite weight (zero-length "
                                  "components are contracted first, Sec 5.1 l.962)");
    }
  }
}

/// Adjacency of a spanning tree given by edge ids: CSR over (neighbour, edge id).
struct TreeAdjacency {
  std::vector<std::size_t> begin;
  std::vector<std::pair<Vertex, EdgeId>> nbr;
};

template <class W>
TreeAdjacency tree_adjacency(std::size_t n, std::span<const Edge<W>> edges,
                             const std::vector<EdgeId>& tree) {
  TreeAdjacency a;
  a.begin.assign(n + 1, 0);
  for (const EdgeId e : tree) {
    ++a.begin[edges[e].u + 1];
    ++a.begin[edges[e].v + 1];
  }
  for (std::size_t v = 0; v < n; ++v) a.begin[v + 1] += a.begin[v];
  a.nbr.resize(a.begin[n]);
  std::vector<std::size_t> pos(a.begin.begin(), a.begin.end() - 1);
  for (const EdgeId e : tree) {
    a.nbr[pos[edges[e].u]++] = {edges[e].v, e};
    a.nbr[pos[edges[e].v]++] = {edges[e].u, e};
  }
  return a;
}

/// BFS of the tree from `root`: discovery order, parent and parent edge (Lemma 3.1).
struct RootedOrder {
  std::vector<Vertex> order;
  std::vector<Vertex> parent;
  std::vector<EdgeId> parent_edge;
};

RootedOrder root_at(const TreeAdjacency& a, std::size_t n, Vertex root) {
  RootedOrder r{{}, std::vector<Vertex>(n, kNoVertex), std::vector<EdgeId>(n, kNoEdge)};
  r.order.reserve(n);
  std::vector<std::uint8_t> seen(n, 0);
  r.order.push_back(root);
  seen[root] = 1;
  for (std::size_t i = 0; i < r.order.size(); ++i) {
    const Vertex x = r.order[i];
    for (std::size_t k = a.begin[x]; k < a.begin[x + 1]; ++k) {
      const auto [y, e] = a.nbr[k];
      if (seen[y] != 0) continue;
      seen[y] = 1;
      r.parent[y] = x;
      r.parent_edge[y] = e;
      r.order.push_back(y);
    }
  }
  return r;
}

}  // namespace

std::size_t local_separator_size(std::size_t n) {
  return static_cast<std::size_t>(model::saturating_pow(model::L(n), kLocalSeparatorExponent));
}

template <class W>
OwnerDecomposition<W> owner_decomposition(const CsrGraph<W>& g, std::size_t b, Diagnostics* diag) {
  const std::size_t n = g.num_vertices();
  if (n == 0) throw std::invalid_argument("owner_decomposition: the graph has no vertex");
  if (b == 0) throw std::invalid_argument("owner_decomposition: separator size b must be >= 1");
  require_positive_weights(g);
  const auto edges = g.edges();

  OwnerDecomposition<W> od;
  od.separator = b;
  od.hop_bound = b == std::numeric_limits<std::size_t>::max() ? b : b + 1;
  od.root = 0;  // [Interpretation, Q-31] the paper does not fix the root

  // Step 1 (l.962-966): MST (Lemma 3.5; throws on a disconnected graph), root distances, and the
  // reduced graph via LCA (Lemma 3.9).
  od.tree_edges = mst_three_pass(g, diag).edges;
  const TreeAdjacency adj = tree_adjacency(n, edges, od.tree_edges);
  RootedOrder rooted = root_at(adj, n, od.root);
  od.tree_parent = rooted.parent;
  od.tree_parent_edge = rooted.parent_edge;
  od.root_distance.assign(n, W{0});
  for (const Vertex v : rooted.order) {
    if (v == od.root) continue;
    od.root_distance[v] = WeightTraits<W>::checked_add(od.root_distance[od.tree_parent[v]],
                                                       edges[od.tree_parent_edge[v]].weight);
  }
  const ForestQueries tree_queries(od.tree_parent);
  std::vector<std::uint8_t> in_tree(edges.size(), 0);
  for (const EdgeId e : od.tree_edges) in_tree[e] = 1;
  od.reduced.assign(edges.size(), 0);
  std::uint64_t dropped = 0;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (in_tree[i] != 0) {
      od.reduced[i] = 1;
      continue;
    }
    bool keep = false;
    if (e.u != e.v) {  // a loop (positive) is always heavier than its empty tree path
      // w <= D_T(u, v) = rd(u) + rd(v) - 2 rd(lca), decided without subtraction.
      const std::optional<Vertex> lca = tree_queries.lca(e.u, e.v);
      if (!lca) throw std::logic_error("owner_decomposition: T is not spanning");
      const W lca_rd = od.root_distance[*lca];
      const W lhs =
          WeightTraits<W>::checked_add(WeightTraits<W>::checked_add(e.weight, lca_rd), lca_rd);
      const W rhs = WeightTraits<W>::checked_add(od.root_distance[e.u], od.root_distance[e.v]);
      keep = !(rhs < lhs);
    }
    od.reduced[i] = keep ? 1 : 0;
    if (!keep) ++dropped;
  }

  // Step 2 (l.968-971): postorder marking with residual sizes; one pivot if nothing is marked.
  od.marked.assign(n, 0);
  std::vector<std::size_t> residual(n, 1);
  std::size_t marks = 0;
  for (std::size_t i = n; i-- > 0;) {
    const Vertex v = rooted.order[i];
    if (residual[v] >= b) {
      od.marked[v] = 1;  // returns size zero to its parent
      ++marks;
    } else if (od.tree_parent[v] != kNoVertex) {
      residual[od.tree_parent[v]] += residual[v];
    }
  }
  const bool unmarked = marks == 0;
  if (unmarked) {
    od.marked[od.root] = 1;  // [Interpretation, Q-31] the single pivot is the root
    marks = 1;
  }

  // Step 3 (l.972-975): branch vertices of the minimal subtree S spanning A0. A vertex lies in
  // S iff it is marked or two of its directions (children, parent side) contain marks; its
  // degree in S is the number of such directions.
  std::vector<std::size_t> below(n, 0);         // marks in the subtree
  std::vector<std::uint32_t> child_dirs(n, 0);  // children whose subtree has a mark
  for (std::size_t i = n; i-- > 0;) {
    const Vertex v = rooted.order[i];
    below[v] += od.marked[v];
    const Vertex p = od.tree_parent[v];
    if (p != kNoVertex) {
      below[p] += below[v];
      if (below[v] > 0) ++child_dirs[p];
    }
  }
  od.pivot_index.assign(n, kNoPivot);
  std::uint64_t branch = 0;
  for (Vertex v = 0; v < n; ++v) {
    const std::uint32_t dirs = child_dirs[v] + (marks - below[v] > 0 ? 1U : 0U);
    const bool is_branch = dirs >= 3;
    if (od.marked[v] != 0 || is_branch) {
      od.pivot_index[v] = static_cast<std::uint32_t>(od.pivots.size());
      od.pivots.push_back(v);
      if (od.marked[v] == 0) ++branch;
    }
  }

  // Step 4 (l.983-985, l.989-991): corridor cuts. Rooted at the smallest pivot, every other
  // pivot walks up to its first pivot ancestor; the walk is the corridor of one two-pivot
  // component (or a single edge between adjacent pivots). Cut a maximum edge, smallest id on
  // ties [Interpretation, Q-31].
  const RootedOrder at_pivot = root_at(adj, n, od.pivots.front());
  od.cut.assign(edges.size(), 0);
  for (std::size_t k = 1; k < od.pivots.size(); ++k) {
    Vertex x = od.pivots[k];
    EdgeId best = kNoEdge;
    do {
      const EdgeId e = at_pivot.parent_edge[x];
      if (best == kNoEdge || edges[best].weight < edges[e].weight ||
          (!(edges[e].weight < edges[best].weight) && e < best)) {
        best = e;
      }
      x = at_pivot.parent[x];
    } while (od.pivot_index[x] == kNoPivot);
    od.corridors.push_back(
        {static_cast<std::uint32_t>(k), od.pivot_index[x], best, edges[best].weight});
    od.cut[best] = 1;
  }

  // Step 5 (l.985-988): owner cells by a search from each pivot over T minus the cut edges.
  od.owner.assign(n, kNoPivot);
  od.radius.assign(n, W{0});
  od.owner_distance.assign(n, W{0});
  od.owner_hops.assign(n, 0);
  std::vector<Vertex> queue;
  queue.reserve(n);
  for (std::uint32_t pi = 0; pi < od.pivots.size(); ++pi) {
    const Vertex p = od.pivots[pi];
    if (od.owner[p] != kNoPivot)
      throw std::logic_error("owner_decomposition: two pivots in a cell");
    od.owner[p] = pi;
    queue.clear();
    queue.push_back(p);
    for (std::size_t i = 0; i < queue.size(); ++i) {
      const Vertex x = queue[i];
      for (std::size_t k = adj.begin[x]; k < adj.begin[x + 1]; ++k) {
        const auto [y, e] = adj.nbr[k];
        if (od.cut[e] != 0 || od.owner[y] == pi) continue;
        if (od.owner[y] != kNoPivot) {
          throw std::logic_error("owner_decomposition: two pivots in a cell");
        }
        od.owner[y] = pi;
        od.radius[y] = std::max(od.radius[x], edges[e].weight);
        od.owner_distance[y] = WeightTraits<W>::checked_add(od.owner_distance[x], edges[e].weight);
        od.owner_hops[y] = od.owner_hops[x] + 1;
        queue.push_back(y);
      }
    }
  }
  for (Vertex v = 0; v < n; ++v) {
    if (od.owner[v] == kNoPivot) throw std::logic_error("owner_decomposition: cell without pivot");
  }

  // Step 6 (l.994-997): sorted merge hierarchy of T_A for sigma(p, q) = b_T(p, q).
  const std::size_t k = od.pivots.size();
  std::vector<std::uint32_t> order(od.corridors.size());
  for (std::uint32_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t c) {
    return MstEdgeKey<W>{od.corridors[a].weight, od.corridors[a].edge} <
           MstEdgeKey<W>{od.corridors[c].weight, od.corridors[c].edge};
  });
  od.merge_parent.assign(2 * k - 1, kNoVertex);
  od.merge_corridor = order;
  DisjointSets dsu(k);
  std::vector<Vertex> node_of(k);
  for (std::size_t i = 0; i < k; ++i) node_of[i] = static_cast<Vertex>(i);
  for (std::size_t i = 0; i < order.size(); ++i) {
    const auto& c = od.corridors[order[i]];
    const auto ra = dsu.find(c.lower);
    const auto rb = dsu.find(c.upper);
    if (ra == rb) throw std::logic_error("owner_decomposition: T_A is not a tree");
    const auto t = static_cast<Vertex>(k + i);
    od.merge_parent[node_of[ra]] = t;
    od.merge_parent[node_of[rb]] = t;
    dsu.unite(ra, rb);
    node_of[dsu.find(ra)] = t;
  }
  od.merge_queries.emplace(od.merge_parent);

  static_cast<void>(dropped);
  static_cast<void>(branch);
  if (diag != nullptr) {
    SSSP_COUNT(++diag->counters["owner_decomposition_calls"]);
    SSSP_COUNT(diag->counters["owner_decomposition_marks"] += unmarked ? 0 : marks);
    SSSP_COUNT(diag->counters["owner_decomposition_pivots"] += k);
    SSSP_COUNT(diag->counters["owner_decomposition_branch_pivots"] += branch);
    SSSP_COUNT(diag->counters["owner_decomposition_unmarked"] += unmarked ? 1 : 0);
    SSSP_COUNT(diag->counters["owner_decomposition_dropped_edges"] += dropped);
  }
  return od;
}

template OwnerDecomposition<std::uint64_t> owner_decomposition<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, std::size_t, Diagnostics*);
template OwnerDecomposition<double> owner_decomposition<double>(const CsrGraph<double>&,
                                                                std::size_t, Diagnostics*);

}  // namespace sssp::paper
