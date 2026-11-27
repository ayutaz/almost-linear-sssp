#pragma once
// Test support for P4-LOCAL-LARGE: Algorithm 4 (Sec 5.4, Lemma 5.10) with the paper's Mode A
// constants at sizes where the radial branch applies (k > 2 L_k^4, the least such k is 209953).
//
// Everything here is written from the paper's definitions and shares no code with src/paper:
//   * dijkstra(): binary-heap Dijkstra on an undirected edge list (physical distances) and on
//     a directed state graph (Lemma 5.10 path realization);
//   * check_def41(): Definition 4.1 (Sec 4, l.646-657) clause by clause:
//       structure  the rooted forest has exactly the states as leaves;
//       H1         every state arc projects to a physical walk no longer than the arc;
//       H2         internal nodes have >= 2 children and a positive height, nondecreasing
//                  towards the root;
//       H3         arcs between distinct children of C with length < h(C) form an acyclic
//                  quotient (each arc is charged to the lowest common ancestor of its endpoint
//                  leaves, found by Tarjan's offline LCA, and to the children below it);
//       H4         the quotient of the roots on all interroot arcs is acyclic;
//       H5         the children of C are connected by a tree of children(C) - 1 links whose
//                  endpoints lie below distinct children, each certifying a physical walk of
//                  length <= link length <= c0 h(C);
//       H6         a link that names a backing state arc names an arc joining its endpoints
//                  that is no longer than the link (certificate-only links name none).
//     Physical walk bounds are decided by "same vertex", then "a direct edge at most the bound"
//     and only otherwise by a Dijkstra cut off at the bound (counted in Def41Report).
//
// Weights: std::uint64_t, or quarter-integer doubles whose sums stay below 2^53 (every sum is
// exact, so double results are compared by exact equality as well).

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/graph.hpp"

namespace sssp::test::local_large {

using paper::ArcId;
using paper::kNoArc;
using paper::kNoNode;
using paper::NodeId;
using paper::OrderedHierarchy;
using paper::StateId;

/// Sec 2 l.284-286: L_N = ceil(log2(N + 2)); 2 * 18^4 = 209952 is the last base-case size
/// (L = 18 for 2^17 - 1 <= N <= 2^18 - 2), 2 * 19^4 = 260642 (L = 19 from N = 2^18 - 1).
inline constexpr std::size_t kLastBaseCase = 209952;
inline constexpr std::size_t kFirstRadial = 209953;
inline constexpr std::size_t kLastL18 = (std::size_t{1} << 18) - 2;   // 262142
inline constexpr std::size_t kFirstL19 = (std::size_t{1} << 18) - 1;  // 262143

template <class W>
constexpr W infinity() {
  return std::numeric_limits<W>::has_infinity ? std::numeric_limits<W>::infinity()
                                              : std::numeric_limits<W>::max();
}

/// Directed adjacency in CSR form (to, length).
template <class W>
struct Digraph {
  std::vector<std::size_t> offset;
  std::vector<std::pair<std::uint32_t, W>> out;
};

/// Undirected physical graph: both orientations of every edge (one for a loop).
template <class W>
Digraph<W> undirected(std::size_t n, const std::vector<Edge<W>>& edges) {
  Digraph<W> d;
  d.offset.assign(n + 1, 0);
  for (const auto& e : edges) {
    ++d.offset[e.u + 1];
    if (e.u != e.v) ++d.offset[e.v + 1];
  }
  std::partial_sum(d.offset.begin(), d.offset.end(), d.offset.begin());
  d.out.resize(d.offset[n]);
  std::vector<std::size_t> cur(d.offset.begin(), d.offset.end() - 1);
  for (const auto& e : edges) {
    d.out[cur[e.u]++] = {e.v, e.weight};
    if (e.u != e.v) d.out[cur[e.v]++] = {e.u, e.weight};
  }
  // Sorted by (neighbour, length) so the lightest direct edge is found by binary search.
  for (std::size_t v = 0; v < n; ++v) {
    std::sort(d.out.begin() + static_cast<std::ptrdiff_t>(d.offset[v]),
              d.out.begin() + static_cast<std::ptrdiff_t>(d.offset[v + 1]));
  }
  return d;
}

/// The state graph D of a hierarchy as a digraph.
template <class W>
Digraph<W> state_digraph(const OrderedHierarchy<W>& h) {
  const std::size_t s = h.graph.num_states;
  Digraph<W> d;
  d.offset.assign(s + 1, 0);
  for (const auto& a : h.graph.arcs) ++d.offset[a.tail + 1];
  std::partial_sum(d.offset.begin(), d.offset.end(), d.offset.begin());
  d.out.resize(d.offset[s]);
  std::vector<std::size_t> cur(d.offset.begin(), d.offset.end() - 1);
  for (const auto& a : h.graph.arcs) d.out[cur[a.tail]++] = {a.head, a.length};
  return d;
}

/// Binary-heap Dijkstra from (vertex, offset) sources; stops once the least key exceeds
/// `bound` (vertices beyond it keep infinity).
template <class W>
std::vector<W> dijkstra(const Digraph<W>& d, const std::vector<std::pair<std::uint32_t, W>>& src,
                        W bound = infinity<W>()) {
  const std::size_t n = d.offset.size() - 1;
  std::vector<W> dist(n, infinity<W>());
  using Item = std::pair<W, std::uint32_t>;
  std::priority_queue<Item, std::vector<Item>, std::greater<>> pq;
  for (const auto& [v, c] : src) {
    if (c < dist[v]) {
      dist[v] = c;
      pq.emplace(c, v);
    }
  }
  while (!pq.empty()) {
    const auto [du, u] = pq.top();
    pq.pop();
    if (du != dist[u]) continue;
    if (bound < du) break;
    for (std::size_t i = d.offset[u]; i < d.offset[u + 1]; ++i) {
      const auto [v, w] = d.out[i];
      const W nd = du + w;  // sums stay far below the type's range in these families
      if (nd < dist[v]) {
        dist[v] = nd;
        pq.emplace(nd, v);
      }
    }
  }
  return dist;
}

struct Def41Report {
  std::vector<std::string> violations;  ///< at most kMaxMessages messages are kept
  std::size_t violation_count = 0;
  std::size_t internal_nodes = 0;
  std::size_t max_depth = 0;
  std::size_t walk_searches = 0;  ///< bounds that needed a cut-off Dijkstra (not same/direct)
  std::size_t links = 0;
  std::size_t certificate_only_links = 0;
  std::size_t light_quotient_arcs = 0;  ///< arcs charged to H3 (length < h(LCA))
  std::size_t interroot_arcs = 0;

  static constexpr std::size_t kMaxMessages = 20;
  [[nodiscard]] bool ok() const { return violation_count == 0; }
  void add(const std::string& m) {
    ++violation_count;
    if (violations.size() < kMaxMessages) violations.push_back(m);
  }
  [[nodiscard]] std::string summary() const {
    std::ostringstream os;
    os << violation_count << " violation(s)";
    for (const auto& v : violations) os << "\n  " << v;
    return os.str();
  }
};

namespace detail {

/// Union-find with path halving.
struct Dsu {
  std::vector<std::uint32_t> p;
  explicit Dsu(std::size_t n) : p(n) { std::iota(p.begin(), p.end(), 0U); }
  std::uint32_t find(std::uint32_t x) {
    while (p[x] != x) {
      p[x] = p[p[x]];
      x = p[x];
    }
    return x;
  }
};

/// Kahn's algorithm on a quotient with `k` vertices; true iff acyclic.
inline bool acyclic(std::size_t k, const std::vector<std::pair<std::uint32_t, std::uint32_t>>& e) {
  std::vector<std::uint32_t> indeg(k, 0);
  std::vector<std::vector<std::uint32_t>> out(k);
  for (const auto& [a, b] : e) {
    out[a].push_back(b);
    ++indeg[b];
  }
  std::vector<std::uint32_t> q;
  for (std::uint32_t v = 0; v < k; ++v) {
    if (indeg[v] == 0) q.push_back(v);
  }
  std::size_t seen = 0;
  while (!q.empty()) {
    const auto v = q.back();
    q.pop_back();
    ++seen;
    for (const auto w : out[v]) {
      if (--indeg[w] == 0) q.push_back(w);
    }
  }
  return seen == k;
}

template <class W>
bool positive_finite(W x) {
  return W{0} < x && x < infinity<W>();
}

template <class W>
bool nonnegative_finite(W x) {
  return !(x < W{0}) && x < infinity<W>() && x == x;
}

}  // namespace detail

/// Definition 4.1 on `h` against the undirected physical graph (n, edges). `c0` is the
/// connection constant to check against (Lemma 5.10: at most two).
template <class W>
Def41Report check_def41(const OrderedHierarchy<W>& h, std::size_t n,
                        const std::vector<Edge<W>>& edges, std::uint32_t c0) {
  Def41Report rep;
  const std::size_t s = h.graph.num_states;
  const std::size_t nn = h.nodes.size();
  const auto phys = undirected(n, edges);

  // ---------------------------------------------------------------- structure
  if (h.projection.size() != s || h.leaf_of_state.size() != s) {
    rep.add("projection / leaf_of_state size differs from the number of states");
    return rep;
  }
  for (StateId x = 0; x < s; ++x) {
    if (h.projection[x] >= n) rep.add("state " + std::to_string(x) + " projects out of range");
    const NodeId l = h.leaf_of_state[x];
    if (l >= nn || !h.nodes[l].children.empty() || h.nodes[l].state != x) {
      rep.add("leaf_of_state[" + std::to_string(x) + "] is not the leaf of that state");
    }
  }
  for (const auto& a : h.graph.arcs) {
    if (a.tail >= s || a.head >= s) rep.add("arc endpoint out of range");
  }
  std::size_t leaves = 0;
  std::size_t child_entries = 0;
  std::size_t with_parent = 0;
  for (NodeId c = 0; c < nn; ++c) {
    const auto& node = h.nodes[c];
    if (node.children.empty()) {
      ++leaves;
      if (node.state >= s || h.leaf_of_state[node.state] != c) {
        rep.add("childless node " + std::to_string(c) + " is not the leaf of its state");
      }
    }
    for (const NodeId ch : node.children) {
      ++child_entries;
      if (ch >= nn || h.nodes[ch].parent != c) {
        rep.add("child " + std::to_string(ch) + " of node " + std::to_string(c) +
                " does not point back");
      }
    }
    if (node.parent != kNoNode) {
      ++with_parent;
      if (node.parent >= nn) rep.add("parent out of range at node " + std::to_string(c));
    }
  }
  if (leaves != s) rep.add("number of leaves differs from the number of states");
  if (child_entries != with_parent) rep.add("child lists and parent pointers disagree");
  {
    std::vector<NodeId> want;
    for (NodeId c = 0; c < nn; ++c) {
      if (h.nodes[c].parent == kNoNode) want.push_back(c);
    }
    auto got = h.roots;
    std::sort(got.begin(), got.end());
    if (got != want) rep.add("roots is not exactly the set of parentless nodes");
  }
  if (!rep.ok()) return rep;

  // Pre-order DFS over the forest (children in list order): tin / tout, tree id, depth. Every
  // node must be reached exactly once (no parent cycles).
  std::vector<std::uint32_t> tin(nn, 0);
  std::vector<std::uint32_t> tout(nn, 0);
  std::vector<std::uint32_t> tree(nn, std::numeric_limits<std::uint32_t>::max());
  std::vector<std::uint32_t> depth(nn, 0);
  std::vector<NodeId> order;  // pre-order
  order.reserve(nn);
  {
    std::uint32_t clock = 0;
    std::vector<std::pair<NodeId, std::size_t>> st;
    for (std::uint32_t r = 0; r < h.roots.size(); ++r) {
      st.emplace_back(h.roots[r], 0);
      tree[h.roots[r]] = r;
      tin[h.roots[r]] = clock++;
      order.push_back(h.roots[r]);
      while (!st.empty()) {
        auto& [c, i] = st.back();
        if (i < h.nodes[c].children.size()) {
          const NodeId ch = h.nodes[c].children[i++];
          if (tree[ch] != std::numeric_limits<std::uint32_t>::max()) {
            rep.add("node " + std::to_string(ch) + " reached twice");
            return rep;
          }
          tree[ch] = r;
          depth[ch] = depth[c] + 1;
          rep.max_depth = std::max<std::size_t>(rep.max_depth, depth[ch]);
          tin[ch] = clock++;
          order.push_back(ch);
          st.emplace_back(ch, 0);
        } else {
          tout[c] = clock - 1;
          st.pop_back();
        }
      }
    }
    if (order.size() != nn) {
      rep.add("parent pointers do not reach a root from every node");
      return rep;
    }
  }
  const auto below = [&](NodeId c, NodeId x) { return tin[c] <= tin[x] && tin[x] <= tout[c]; };
  // Index of the child of internal node c whose subtree contains node x (x strictly below c).
  const auto child_index = [&](NodeId c, NodeId x) {
    const auto& ch = h.nodes[c].children;
    const auto it = std::upper_bound(ch.begin(), ch.end(), tin[x],
                                     [&](std::uint32_t t, NodeId y) { return t < tin[y]; });
    return static_cast<std::uint32_t>((it - ch.begin()) - 1);
  };
  const auto walk_within = [&](Vertex u, Vertex v, W bound) {
    if (u == v) return true;
    const auto first = phys.out.begin() + static_cast<std::ptrdiff_t>(phys.offset[u]);
    const auto last = phys.out.begin() + static_cast<std::ptrdiff_t>(phys.offset[u + 1]);
    const auto it = std::lower_bound(first, last, std::pair<std::uint32_t, W>{v, W{0}});
    if (it != last && it->first == v) return !(bound < it->second);
    ++rep.walk_searches;
    const auto d = dijkstra(phys, {{u, W{0}}}, bound);
    return !(bound < d[v]);
  };

  // ---------------------------------------------------------------- H1
  for (ArcId a = 0; a < h.graph.arcs.size(); ++a) {
    const auto& arc = h.graph.arcs[a];
    if (!detail::nonnegative_finite(arc.length)) {
      rep.add("H1: arc " + std::to_string(a) + " has an invalid length");
    } else if (!walk_within(h.projection[arc.tail], h.projection[arc.head], arc.length)) {
      rep.add("H1: arc " + std::to_string(a) + " is shorter than the physical distance");
    }
  }

  // ---------------------------------------------------------------- H2
  for (NodeId c = 0; c < nn; ++c) {
    const auto& node = h.nodes[c];
    if (node.children.empty()) continue;
    ++rep.internal_nodes;
    if (node.children.size() < 2) rep.add("H2: node " + std::to_string(c) + " has one child");
    if (!detail::positive_finite(node.height)) {
      rep.add("H2: node " + std::to_string(c) + " has a non-positive or infinite height");
    }
    for (const NodeId ch : node.children) {
      if (!h.nodes[ch].children.empty() && node.height < h.nodes[ch].height) {
        rep.add("H2: child " + std::to_string(ch) + " is higher than its parent " +
                std::to_string(c));
      }
    }
  }

  // ---------------------------------------------------------------- H3 / H4
  // Tarjan's offline LCA of the endpoint leaves of every intratree arc.
  std::vector<NodeId> lca(h.graph.arcs.size(), kNoNode);
  {
    std::vector<std::vector<ArcId>> at(s);
    for (ArcId a = 0; a < h.graph.arcs.size(); ++a) {
      const auto& arc = h.graph.arcs[a];
      if (arc.tail == arc.head) continue;  // a loop lies inside one leaf
      if (tree[h.leaf_of_state[arc.tail]] != tree[h.leaf_of_state[arc.head]]) continue;
      at[arc.tail].push_back(a);
      at[arc.head].push_back(a);
    }
    detail::Dsu dsu(nn);
    std::vector<NodeId> anc(nn);
    std::iota(anc.begin(), anc.end(), 0U);
    std::vector<char> done(nn, 0);
    // Post-order: children are finished before their parent; process in reverse pre-order is
    // not enough for Tarjan, so run an explicit DFS again.
    std::vector<std::pair<NodeId, std::size_t>> st;
    for (const NodeId r : h.roots) {
      st.emplace_back(r, 0);
      while (!st.empty()) {
        auto& [c, i] = st.back();
        if (i < h.nodes[c].children.size()) {
          st.emplace_back(h.nodes[c].children[i++], 0);
          continue;
        }
        const NodeId fin = c;
        st.pop_back();
        done[fin] = 1;
        if (h.nodes[fin].children.empty()) {
          const StateId x = h.nodes[fin].state;
          for (const ArcId a : at[x]) {
            const auto& arc = h.graph.arcs[a];
            const StateId other = arc.tail == x ? arc.head : arc.tail;
            const NodeId lo = h.leaf_of_state[other];
            if (done[lo] != 0) lca[a] = anc[dsu.find(lo)];
          }
        }
        if (!st.empty()) {
          const NodeId par = st.back().first;
          dsu.p[dsu.find(fin)] = dsu.find(par);
          anc[dsu.find(par)] = par;
        }
      }
    }
  }
  {
    // Charge every light intratree arc to its LCA; group by LCA node.
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> quotient(nn);
    std::vector<std::pair<std::uint32_t, std::uint32_t>> roots_q;
    for (ArcId a = 0; a < h.graph.arcs.size(); ++a) {
      const auto& arc = h.graph.arcs[a];
      const NodeId lx = h.leaf_of_state[arc.tail];
      const NodeId ly = h.leaf_of_state[arc.head];
      if (tree[lx] != tree[ly]) {
        ++rep.interroot_arcs;
        roots_q.emplace_back(tree[lx], tree[ly]);
        continue;
      }
      if (arc.tail == arc.head) continue;
      const NodeId c = lca[a];
      if (c == kNoNode || !below(c, lx) || !below(c, ly) || h.nodes[c].children.empty()) {
        rep.add("internal: LCA of arc " + std::to_string(a) + " not found");
        continue;
      }
      if (!(arc.length < h.nodes[c].height)) continue;
      ++rep.light_quotient_arcs;
      const auto ix = child_index(c, lx);
      const auto iy = child_index(c, ly);
      if (ix == iy) {
        rep.add("internal: arc " + std::to_string(a) + " charged below one child");
        continue;
      }
      quotient[c].emplace_back(ix, iy);
    }
    for (NodeId c = 0; c < nn; ++c) {
      if (quotient[c].empty()) continue;
      if (!detail::acyclic(h.nodes[c].children.size(), quotient[c])) {
        rep.add("H3: light-arc quotient of node " + std::to_string(c) + " has a cycle");
      }
    }
    if (!detail::acyclic(h.roots.size(), roots_q)) rep.add("H4: interroot quotient has a cycle");
  }

  // ---------------------------------------------------------------- H5 / H6
  detail::Dsu link_dsu(nn);
  for (NodeId c = 0; c < nn; ++c) {
    const auto& node = h.nodes[c];
    if (node.children.empty()) {
      if (!node.links.empty()) rep.add("H5: leaf " + std::to_string(c) + " carries links");
      continue;
    }
    const std::string at = " at node " + std::to_string(c);
    if (node.links.size() + 1 != node.children.size()) rep.add("H5: wrong link count" + at);
    // c0 h(C), saturated for integers.
    W cap = node.height;
    if constexpr (std::is_integral_v<W>) {
      cap = node.height > infinity<W>() / c0 ? infinity<W>() : node.height * c0;
    } else {
      cap = node.height * static_cast<W>(c0);
    }
    for (std::size_t li = 0; li < node.links.size(); ++li) {
      const auto& link = node.links[li];
      ++rep.links;
      const std::string where = " (link " + std::to_string(li) + at + ")";
      if (link.x >= s || link.y >= s) {
        rep.add("H5: link endpoint is not a state" + where);
        continue;
      }
      const NodeId lx = h.leaf_of_state[link.x];
      const NodeId ly = h.leaf_of_state[link.y];
      if (!below(c, lx) || !below(c, ly) || lx == c || ly == c) {
        rep.add("H5: link endpoint is not below the node" + where);
        continue;
      }
      const auto ix = child_index(c, lx);
      const auto iy = child_index(c, ly);
      if (ix == iy) {
        rep.add("H5: both link endpoints below one child" + where);
        continue;
      }
      const auto rx = link_dsu.find(node.children[ix]);
      const auto ry = link_dsu.find(node.children[iy]);
      if (rx == ry) {
        rep.add("H5: links form a cycle" + where);
      } else {
        link_dsu.p[rx] = ry;
      }
      if (!detail::nonnegative_finite(link.length)) {
        rep.add("H5: invalid link length" + where);
        continue;
      }
      if (cap < link.length) rep.add("H5: link longer than c0 h(C)" + where);
      if (!walk_within(h.projection[link.x], h.projection[link.y], link.length)) {
        rep.add("H5: link shorter than the physical distance" + where);
      }
      if (link.arc == kNoArc) {
        ++rep.certificate_only_links;
      } else if (link.arc >= h.graph.arcs.size()) {
        rep.add("H6: backing arc out of range" + where);
      } else {
        const auto& arc = h.graph.arcs[link.arc];
        const bool joins = (arc.tail == link.x && arc.head == link.y) ||
                           (arc.tail == link.y && arc.head == link.x);
        if (!joins || link.length < arc.length) rep.add("H6: backing arc mismatch" + where);
      }
    }
    // Spanning: every child in one link component.
    const auto r0 = link_dsu.find(node.children[0]);
    for (const NodeId ch : node.children) {
      if (link_dsu.find(ch) != r0) {
        rep.add("H5: links do not connect the children" + at);
        break;
      }
    }
  }
  return rep;
}

}  // namespace sssp::test::local_large
