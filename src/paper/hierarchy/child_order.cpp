// Ordering the children (Sec 4.2, l.725-779, Lemma 4.4, Lemma 4.5). See child_order.hpp.

#include "paper/hierarchy/child_order.hpp"

#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "paper/model/length.hpp"
#include "paper/primitives/forest_queries.hpp"
#include "paper/primitives/search.hpp"

namespace sssp::paper {

namespace {

static_assert(kNoState == kNoVertex, "F parent arrays reuse kNoState as the root marker");

constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

template <class W>
class Arith {
 public:
  explicit Arith(Diagnostics* diag) : ops_(diag != nullptr ? &diag->ops : nullptr) {}
  [[nodiscard]] W add(W a, W b) const { return ops_.add(L::finite(a), L::finite(b)).value(); }

 private:
  using L = model::Len<W>;
  model::LenOps<W> ops_;
};

void bump([[maybe_unused]] Diagnostics* diag, [[maybe_unused]] const char* name,
          [[maybe_unused]] std::uint64_t k) {
  SSSP_COUNT(if (diag != nullptr) diag->counters[name] += k);
}

/// The running charge of the traversal as a signed pair (Sec 2, l.294-296): value = pos - neg.
template <class W>
struct Charge {
  W pos{0};
  W neg{0};
};

/// [Deviation candidate D-15, see child_order.hpp] One real from a signed pair. Both components
/// are monotone roundings of termwise-ordered sums (pos term >= neg term), so pos >= neg.
template <class W>
W materialise(const Charge<W>& c) {
  if (c.pos < c.neg) throw std::logic_error("order_children: negative charge");
  return c.pos - c.neg;
}

}  // namespace

template <class W>
LinkForest<W> build_link_forest(const NormalizedHierarchy<W>& nh, Diagnostics* diag) {
  const Arith<W> ar(diag);
  const auto& f = nh.forest;
  const std::size_t n = f.num_states();

  // Undirected adjacency of the links, in link order (one counting pass, Lemma 3.1).
  struct Half {
    StateId to;
    W len;
  };
  std::vector<std::size_t> begin(n + 1, 0);
  std::size_t num_links = 0;
  for (const auto& nd : f.nodes) {
    for (const auto& lk : nd.links) {
      if (lk.x >= n || lk.y >= n) throw std::invalid_argument("build_link_forest: bad endpoint");
      ++begin[lk.x + 1];
      ++begin[lk.y + 1];
      ++num_links;
    }
  }
  for (std::size_t s = 0; s < n; ++s) begin[s + 1] += begin[s];
  std::vector<Half> adj(begin[n]);
  {
    std::vector<std::size_t> pos(begin.begin(), begin.end() - 1);
    for (const auto& nd : f.nodes) {
      for (const auto& lk : nd.links) {
        adj[pos[lk.x]++] = {lk.y, lk.length};
        adj[pos[lk.y]++] = {lk.x, lk.length};
      }
    }
  }

  LinkForest<W> lf;
  lf.parent.assign(n, kNoState);
  lf.parent_length.assign(n, W{0});
  lf.root_distance.assign(n, W{0});
  lf.preorder.assign(n, std::numeric_limits<std::uint32_t>::max());
  lf.order.reserve(n);
  std::size_t trees = 0;
  std::vector<std::pair<StateId, std::size_t>> stack;
  for (StateId s = 0; s < n; ++s) {
    if (lf.preorder[s] != std::numeric_limits<std::uint32_t>::max()) continue;
    ++trees;  // rooted at its smallest state
    lf.preorder[s] = static_cast<std::uint32_t>(lf.order.size());
    lf.order.push_back(s);
    stack.emplace_back(s, begin[s]);
    while (!stack.empty()) {
      auto& [u, i] = stack.back();
      if (i == begin[u + 1]) {
        stack.pop_back();
        continue;
      }
      const Half e = adj[i++];
      if (lf.preorder[e.to] != std::numeric_limits<std::uint32_t>::max()) continue;
      const StateId from = u;
      lf.parent[e.to] = from;
      lf.parent_length[e.to] = e.len;
      lf.root_distance[e.to] = ar.add(lf.root_distance[from], e.len);
      lf.preorder[e.to] = static_cast<std::uint32_t>(lf.order.size());
      lf.order.push_back(e.to);
      stack.emplace_back(e.to, begin[e.to]);
    }
  }
  // A spanning forest of n states with `trees` trees has n - trees edges; more links than that
  // close a cycle (Def 4.1 H5 says each node's links form a tree on its children).
  if (num_links + trees != n) {
    throw std::invalid_argument("build_link_forest: the links do not form a forest");
  }
  return lf;
}

template <class W>
void order_children(NormalizedHierarchy<W>& nh, Diagnostics* diag) {
  if (!nh.queries.has_value()) {
    throw std::invalid_argument("order_children: hierarchy was not produced by normalize()");
  }
  const Arith<W> ar(diag);
  auto& f = nh.forest;
  const std::size_t n = f.num_states();
  const std::size_t nf = f.nodes.size();
  const auto& arcs = f.graph.arcs;

  const LinkForest<W> lf = build_link_forest(nh, diag);
  const ForestQueries fq(lf.parent);  // Lemma 3.9 on F
  const auto& pre = lf.preorder;
  std::uint64_t fdist_queries = 0;

  // Charge of an F move a -> b: rd(a) + rd(b) - 2 rd(lca_F(a, b)) as a signed pair.
  const auto move_along_f = [&](Charge<W>& c, StateId a, StateId b) {
    if (a == b) return;
    const auto l = fq.lca(a, b);
    if (!l) throw std::invalid_argument("order_children: links do not connect a node");
    ++fdist_queries;
    c.pos = ar.add(c.pos, ar.add(lf.root_distance[a], lf.root_distance[b]));
    c.neg = ar.add(c.neg, ar.add(lf.root_distance[*l], lf.root_distance[*l]));
  };
  const auto cross_arc = [&](Charge<W>& c, W len) { c.pos = ar.add(c.pos, len); };

  // Top-down order of the forest, then provisional representatives bottom-up: the state with the
  // least F preorder below each node (Q-25).
  std::vector<NodeId> order;
  order.reserve(nf);
  {
    std::vector<NodeId> stack(f.roots.rbegin(), f.roots.rend());
    while (!stack.empty()) {
      const NodeId v = stack.back();
      stack.pop_back();
      order.push_back(v);
      for (const NodeId c : f.nodes[v].children) stack.push_back(c);
    }
  }
  std::vector<StateId> prov(nf, kNoState);
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const auto& nd = f.nodes[*it];
    if (nd.children.empty()) {
      prov[*it] = nd.state;
      continue;
    }
    StateId best = kNoState;
    for (const NodeId c : nd.children) {
      if (best == kNoState || pre[prov[c]] < pre[best]) best = prov[c];
    }
    prov[*it] = best;
  }

  // DFS roots of every node: children sorted by (parent, F preorder of the provisional
  // representative), Lemma 4.5 l.769-772.
  std::vector<NodeId> non_roots;
  for (NodeId v = 0; v < nf; ++v) {
    if (f.nodes[v].parent != kNoNode) non_roots.push_back(v);
  }
  {
    const std::array<std::size_t, 2> ranges{nf, n};
    non_roots = radix_sort_stable(non_roots, ranges, [&](NodeId v, std::size_t j) -> std::size_t {
      return j == 0 ? f.nodes[v].parent : pre[prov[v]];
    });
  }
  std::vector<std::size_t> roots_begin(nf + 1, 0);
  for (const NodeId v : non_roots) ++roots_begin[f.nodes[v].parent + 1];
  for (std::size_t x = 0; x < nf; ++x) roots_begin[x + 1] += roots_begin[x];

  // Cyclic-start table, Lemma 4.5 l.773-777. Records: tails of light arcs (kind 1) and entry
  // positions (kind 0): provisional representatives (id = child) and heads of light arcs
  // (id = nf + light-arc index), sorted by (child, F preorder, kind).
  struct Rec {
    NodeId child;
    std::uint32_t pre;
    std::uint8_t kind;  // 0 = entry, 1 = tail
    std::size_t id;     // entry id, or light-arc index for a tail
  };
  const std::size_t num_light = nh.light_arcs.size();
  std::vector<Rec> recs;
  recs.reserve(non_roots.size() + 2 * num_light);
  for (const NodeId v : non_roots) recs.push_back({v, pre[prov[v]], 0, v});
  for (std::size_t li = 0; li < num_light; ++li) {
    const LightArc& la = nh.light_arcs[li];
    recs.push_back({la.tail_child, pre[arcs[la.arc].tail], 1, li});
    recs.push_back({la.head_child, pre[arcs[la.arc].head], 0, nf + li});
  }
  {
    const std::array<std::size_t, 3> ranges{nf, n, 2};
    recs = radix_sort_stable(recs, ranges, [](const Rec& r, std::size_t j) -> std::size_t {
      return j == 0 ? r.child : (j == 1 ? r.pre : r.kind);
    });
  }
  // Outgoing light arcs of each child in F preorder of their tails: outs[out_begin[i] ..].
  std::vector<std::size_t> outs;
  outs.reserve(num_light);
  std::vector<std::size_t> out_begin(nf + 1, 0);
  std::vector<std::size_t> rec_out_pos(recs.size(), kNone);
  for (std::size_t r = 0; r < recs.size(); ++r) {
    if (recs[r].kind == 1) {
      rec_out_pos[r] = outs.size();
      outs.push_back(recs[r].id);
      ++out_begin[recs[r].child + 1];
    }
  }
  for (std::size_t x = 0; x < nf; ++x) out_begin[x + 1] += out_begin[x];
  // start[entry] = position in `outs` of the first tail at or after the entry in F preorder,
  // wrapping to the child's first tail; kNone (end marker) for a child without outgoing arcs.
  std::vector<std::size_t> start(nf + num_light, kNone);
  {
    std::size_t next_tail = kNone;
    for (std::size_t r = recs.size(); r-- > 0;) {
      const Rec& rc = recs[r];
      if (r + 1 == recs.size() || recs[r + 1].child != rc.child) next_tail = kNone;
      if (rc.kind == 1) {
        next_tail = rec_out_pos[r];
      } else if (next_tail != kNone) {
        start[rc.id] = next_tail;
      } else if (out_begin[rc.child] < out_begin[rc.child + 1]) {
        start[rc.id] = out_begin[rc.child];  // wrap
      }
    }
  }

  // The DFS of Sec 4.2 at every internal node, with the charging of Lemma 4.4.
  nh.rank.assign(nf, 0);
  nh.rep.assign(nf, kNoState);
  nh.gap_next.assign(nf, W{0});
  std::vector<std::uint8_t> status(nf, 0);  // 0 unvisited, 1 on the DFS stack, 2 finished
  struct Frame {
    NodeId child;
    StateId entry;
    std::size_t start;  // position in outs, or kNone
    std::size_t seen;   // outgoing arcs inspected so far
    std::size_t back;   // light-arc index being returned over, or kNone
  };
  std::vector<Frame> stack;
  std::vector<NodeId> finish;
  std::uint64_t internal_nodes = 0;
  std::uint64_t dfs_roots = 0;
  std::uint64_t tree_arcs = 0;
  std::uint64_t inspected = 0;

  for (NodeId x = 0; x < nf; ++x) {
    auto& children = f.nodes[x].children;
    if (children.empty()) continue;
    ++internal_nodes;
    finish.clear();
    std::vector<W> gap_of_finish;  // gap_of_finish[m] = g(f_{m-1}, f_m) for m >= 1
    gap_of_finish.reserve(children.size());
    Charge<W> seg;  // charge since the last mark
    StateId cur = kNoState;

    const auto enter = [&](NodeId i, StateId entry, std::size_t entry_id) {
      status[i] = 1;
      nh.rep[i] = entry;
      stack.push_back({i, entry, start[entry_id], 0, kNone});
    };
    const auto out_count = [&](NodeId i) { return out_begin[i + 1] - out_begin[i]; };

    for (std::size_t k = roots_begin[x]; k < roots_begin[x + 1]; ++k) {
      const NodeId r = non_roots[k];
      if (status[r] != 0) continue;
      ++dfs_roots;
      // Successive DFS roots: move along F between provisional representatives (l.747).
      if (cur != kNoState) move_along_f(seg, cur, prov[r]);
      cur = prov[r];
      enter(r, prov[r], r);
      while (!stack.empty()) {
        Frame& fr = stack.back();
        if (fr.back != kNone) {
          // Return over the tree arc after its head's search finished (l.748).
          const auto& arc = arcs[nh.light_arcs[fr.back].arc];
          cross_arc(seg, arc.length);
          cur = arc.tail;
          fr.back = kNone;
        }
        const std::size_t cnt = out_count(fr.child);
        if (fr.start != kNone && fr.seen < cnt) {
          // Cyclic F preorder of the tails, from the entry (l.729-730).
          std::size_t pos = fr.start + fr.seen;
          if (pos >= out_begin[fr.child + 1]) pos -= cnt;
          ++fr.seen;
          ++inspected;
          const std::size_t li = outs[pos];
          const NodeId j = nh.light_arcs[li].head_child;
          if (status[j] == 2) continue;  // not a tree arc
          if (status[j] == 1) {
            throw std::invalid_argument("order_children: light-arc quotient of node " +
                                        std::to_string(x) + " has a cycle (Def 4.1 H3)");
          }
          const auto& arc = arcs[nh.light_arcs[li].arc];
          move_along_f(seg, cur, arc.tail);  // move inside F[i] to the exit tail
          cross_arc(seg, arc.length);
          cur = arc.head;
          fr.back = li;
          ++tree_arcs;
          enter(j, arc.head, nf + li);  // invalidates fr
          continue;
        }
        // Finish: return to the entry along F, then mark (l.746).
        move_along_f(seg, cur, fr.entry);
        cur = fr.entry;
        const NodeId i = fr.child;
        status[i] = 2;
        if (!finish.empty()) gap_of_finish.push_back(materialise(seg));
        seg = Charge<W>{};
        finish.push_back(i);
        stack.pop_back();
      }
    }
    if (finish.size() != children.size()) {
      throw std::logic_error("order_children: DFS did not finish every child");
    }
    // Ranks: reverse finishing order (Lemma 4.5). g(f_m, f_{m-1}) belongs to f_m, which is
    // ranked immediately before f_{m-1}.
    const std::size_t d = finish.size();
    for (std::size_t k = 0; k < d; ++k) {
      const NodeId c = finish[d - 1 - k];
      children[k] = c;
      nh.rank[c] = static_cast<std::uint32_t>(k);
    }
    nh.gap_next[finish[0]] = W{0};
    for (std::size_t m = 1; m < d; ++m) nh.gap_next[finish[m]] = gap_of_finish[m - 1];
  }
  for (std::size_t i = 0; i < nh.root_order.size(); ++i) {
    nh.rank[nh.root_order[i]] = static_cast<std::uint32_t>(i);
  }
  nh.children_ordered = true;

  bump(diag, "child_order_calls", 1);
  bump(diag, "child_order_nodes", internal_nodes);
  bump(diag, "child_order_light_arcs", num_light);
  bump(diag, "child_order_dfs_roots", dfs_roots);
  bump(diag, "child_order_tree_arcs", tree_arcs);
  bump(diag, "child_order_arcs_inspected", inspected);
  bump(diag, "child_order_forest_distances", fdist_queries);
}

template LinkForest<std::uint64_t> build_link_forest<std::uint64_t>(
    const NormalizedHierarchy<std::uint64_t>&, Diagnostics*);
template LinkForest<double> build_link_forest<double>(const NormalizedHierarchy<double>&,
                                                      Diagnostics*);
template void order_children<std::uint64_t>(NormalizedHierarchy<std::uint64_t>&, Diagnostics*);
template void order_children<double>(NormalizedHierarchy<double>&, Diagnostics*);

}  // namespace sssp::paper
