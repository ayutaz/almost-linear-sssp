// Normalisation of an ordered hierarchy (Sec 4.1, l.674-723, Lemma 4.3). See normalize.hpp.

#include "paper/hierarchy/normalize.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "paper/local/saturating.hpp"
#include "paper/model/length.hpp"
#include "paper/primitives/search.hpp"

namespace sssp::paper {

namespace {

template <class W>
class Arith {
 public:
  explicit Arith(Diagnostics* diag) : ops_(diag != nullptr ? &diag->ops : nullptr) {}

  [[nodiscard]] W add(W a, W b) const { return ops_.add(L::finite(a), L::finite(b)).value(); }
  /// [Deviation D-23] a + b, or +infinity if an operand is +infinity or the sum does not fit
  /// (sat::add, the D-13 helper). One counted addition if both operands are finite.
  [[nodiscard]] W sadd(W a, W b) const {
    if (a == kInf || b == kInf) return kInf;
    SSSP_COUNT(if (ops_.counters() != nullptr)++ ops_.counters()->additions);
    return sat::add(a, b);
  }
  [[nodiscard]] bool less(W a, W b) const { return ops_.less(L::finite(a), L::finite(b)); }
  [[nodiscard]] bool equal(W a, W b) const { return ops_.equal(L::finite(a), L::finite(b)); }

  /// 2a by one doubling (Sec 2, l.297), or nullopt if 2a is not representable (saturation, see
  /// the [Supplement] note in normalize.hpp).
  [[nodiscard]] std::optional<W> twice(W a) const {
    try {
      return ops_.mul_pow2(L::finite(a), 1).value();
    } catch (const std::overflow_error&) {
      return std::nullopt;
    }
  }

  /// 4a >= b, with an unrepresentable 4a (including a saturated a) counting as larger than b.
  [[nodiscard]] bool four_times_at_least(W a, W b) const {
    if (a == kInf) return true;
    const auto a2 = twice(a);
    if (!a2) return true;
    const auto a4 = twice(*a2);
    if (!a4) return true;
    return !less(*a4, b);
  }

 private:
  using L = model::Len<W>;
  static constexpr W kInf = WeightTraits<W>::infinity();
  model::LenOps<W> ops_;
};

void bump([[maybe_unused]] Diagnostics* diag, [[maybe_unused]] const char* name,
          [[maybe_unused]] std::uint64_t k) {
  SSSP_COUNT(if (diag != nullptr) diag->counters[name] += k);
}

}  // namespace

template <class W>
NormalizedHierarchy<W> normalize(const OrderedHierarchy<W>& h, Diagnostics* diag) {
  const Arith<W> ar(diag);
  if (h.c0 == 0) throw std::invalid_argument("normalize: c0 must be positive (Def 4.1)");
  if (h.c0 > std::numeric_limits<std::uint32_t>::max() / 2) {
    throw std::invalid_argument("normalize: c0 too large to double");
  }
  const std::size_t n = h.nodes.size();

  // Top-down order of the input forest (every parent before its children).
  std::vector<NodeId> order;
  order.reserve(n);
  {
    std::vector<NodeId> stack(h.roots.rbegin(), h.roots.rend());
    while (!stack.empty()) {
      const NodeId v = stack.back();
      stack.pop_back();
      if (v >= n || order.size() == n) {
        throw std::invalid_argument("normalize: malformed forest");
      }
      order.push_back(v);
      const auto& ch = h.nodes[v].children;
      for (auto it = ch.rbegin(); it != ch.rend(); ++it) stack.push_back(*it);
    }
    if (order.size() != n) throw std::invalid_argument("normalize: node unreachable from roots");
  }
  const auto internal = [&](NodeId v) { return !h.nodes[v].children.empty(); };

  // Step 1 (l.686-688): masses bottom-up. S(C) = h(C) + its links + S of internal children.
  // [Deviation D-23] A mass that does not fit W saturates at +infinity: S(C) is only an upper
  // bound (on the physical diameter of C, Lemma 4.3), and +infinity is still one.
  std::vector<W> mass(n, W{0});
  std::uint64_t saturated = 0;
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const NodeId v = *it;
    if (!internal(v)) continue;
    const auto& nd = h.nodes[v];
    if (!ar.less(W{0}, nd.height)) {
      throw std::invalid_argument("normalize: node " + std::to_string(v) +
                                  " has a non-positive height (Def 4.1 H2)");
    }
    W s = nd.height;
    for (const auto& lk : nd.links) s = ar.sadd(s, lk.length);
    for (const NodeId c : nd.children) {
      if (internal(c)) s = ar.sadd(s, mass[c]);
    }
    if (s == WeightTraits<W>::infinity()) ++saturated;
    mass[v] = s;
  }

  // Step 2 (l.688-689): normalisation links 4 S(C) >= h(P), and the link set (its top) of
  // every internal node.
  std::vector<std::uint8_t> linked(n, 0);
  std::vector<NodeId> top(n, kNoNode);
  std::uint64_t num_links = 0;
  std::uint64_t num_sets = 0;
  for (const NodeId v : order) {
    if (!internal(v)) continue;
    const NodeId p = h.nodes[v].parent;
    if (p != kNoNode && ar.four_times_at_least(mass[v], h.nodes[p].height)) {
      linked[v] = 1;
      top[v] = top[p];
      ++num_links;
    } else {
      top[v] = v;
      ++num_sets;
    }
  }

  // Step 3 (l.689-692, proof of Lemma 4.3 l.704-706): least height of each set, its doubling
  // table up to the top height, and nu by walking the table downwards from the parent's index.
  std::vector<W> hmin(n, W{0});
  for (const NodeId v : order) {
    if (!internal(v)) continue;
    const NodeId t = top[v];
    if (t == v || ar.less(h.nodes[v].height, hmin[t])) hmin[t] = h.nodes[v].height;
  }
  std::vector<W> table;
  std::vector<std::size_t> idx(n, 0);  // absolute position in `table`
  std::vector<std::size_t> table_begin(n, 0);
  std::vector<W> unit(n, W{0});
  std::uint64_t table_steps = 0;
  for (const NodeId v : order) {
    if (!internal(v)) continue;
    if (top[v] == v) {
      // Doubling from the minimum generates the whole table of the set (l.705).
      table_begin[v] = table.size();
      W t = hmin[v];
      table.push_back(t);
      for (;;) {
        const auto next = ar.twice(t);
        if (!next || ar.less(h.nodes[v].height, *next)) break;
        t = *next;
        table.push_back(t);
        ++table_steps;
      }
      idx[v] = table.size() - 1;
    } else {
      // Move downward from the parent's index, stopping at the child's height (l.705-706).
      const NodeId p = h.nodes[v].parent;
      std::size_t i = idx[p];
      while (ar.less(h.nodes[v].height, table[i])) {
        --i;  // never passes table_begin[top]: table[begin] = hmin <= h(v)
        ++table_steps;
      }
      idx[v] = i;
    }
    unit[v] = table[idx[v]];
  }

  // Step 4 (l.707-713): contract normalisation links with equal units. The survivor of a group
  // is its topmost member; it keeps its own (topmost original) mass.
  std::vector<std::uint8_t> absorbed(n, 0);
  for (const NodeId v : order) {
    if (linked[v] != 0 && ar.equal(unit[v], unit[h.nodes[v].parent])) absorbed[v] = 1;
  }
  NormalizedHierarchy<W> out;
  std::vector<NodeId> new_id(n, kNoNode);
  for (NodeId v = 0; v < n; ++v) {
    if (absorbed[v] == 0) {
      new_id[v] = static_cast<NodeId>(out.original_node.size());
      out.original_node.push_back(v);
    }
  }
  const std::size_t nf = out.original_node.size();
  out.node_of_input.assign(n, kNoNode);
  for (const NodeId v : order) {
    out.node_of_input[v] = absorbed[v] != 0 ? out.node_of_input[h.nodes[v].parent] : new_id[v];
  }

  auto& f = out.forest;
  f.graph = h.graph;
  f.projection = h.projection;
  f.c0 = 2 * h.c0;
  f.nodes.resize(nf);
  out.unit.assign(nf, W{0});
  out.mass.assign(nf, W{0});
  out.original_height.assign(nf, W{0});
  for (NodeId x = 0; x < nf; ++x) {
    const NodeId d = out.original_node[x];
    const auto& src = h.nodes[d];
    auto& dst = f.nodes[x];
    if (!internal(d)) {
      dst.state = src.state;
      continue;
    }
    out.unit[x] = unit[d];
    out.mass[x] = mass[d];
    out.original_height[x] = src.height;
    dst.height = unit[d];
    // Expand absorbed children in place: their children and links move up (l.710).
    dst.links = src.links;
    std::vector<std::pair<NodeId, std::size_t>> stack{{d, 0}};
    while (!stack.empty()) {
      auto& [m, i] = stack.back();
      if (i == h.nodes[m].children.size()) {
        stack.pop_back();
        continue;
      }
      const NodeId c = h.nodes[m].children[i++];
      if (absorbed[c] != 0) {
        dst.links.insert(dst.links.end(), h.nodes[c].links.begin(), h.nodes[c].links.end());
        stack.emplace_back(c, 0);
      } else {
        dst.children.push_back(new_id[c]);
      }
    }
    for (const NodeId c : dst.children) f.nodes[c].parent = x;
  }
  f.leaf_of_state.resize(h.leaf_of_state.size());
  for (std::size_t s = 0; s < h.leaf_of_state.size(); ++s) {
    f.leaf_of_state[s] = new_id[h.leaf_of_state[s]];
  }
  collect_roots(f);

  // Step 5 (l.719-722): discard self-arcs, root quotient, light arcs at the new LCA.
  std::vector<Vertex> parent(nf);
  for (NodeId x = 0; x < nf; ++x) parent[x] = f.nodes[x].parent;
  out.queries.emplace(parent);
  const ForestQueries& q = *out.queries;
  struct Rec {
    NodeId lca;
    LightArc arc;
  };
  std::vector<Rec> recs;
  std::uint64_t self_arcs = 0;
  for (std::size_t a = 0; a < f.graph.arcs.size(); ++a) {
    const auto& arc = f.graph.arcs[a];
    if (arc.tail == arc.head) {
      ++self_arcs;  // l.719: removing a traversal back to the same state never hurts
      continue;
    }
    const NodeId lx = f.leaf_of_state[arc.tail];
    const NodeId ly = f.leaf_of_state[arc.head];
    const auto l = q.lca(lx, ly);
    if (!l) {
      out.root_arcs.push_back(static_cast<ArcId>(a));  // l.720: straight to the root quotient
      continue;
    }
    if (ar.less(arc.length, out.unit[*l])) {
      const std::uint32_t below = q.depth(*l) + 1;
      recs.push_back({*l, LightArc{static_cast<ArcId>(a), q.level_ancestor(lx, below),
                                   q.level_ancestor(ly, below)}});
    }
  }
  // Group the light arcs by their LCA with one stable counting sort (Lemma 3.1).
  const auto sorted = counting_sort_stable(recs, nf, [](const Rec& r) { return r.lca; });
  out.light_begin.assign(nf + 1, 0);
  out.light_arcs.reserve(sorted.size());
  for (const Rec& r : sorted) {
    ++out.light_begin[r.lca + 1];
    out.light_arcs.push_back(r.arc);
  }
  for (std::size_t x = 0; x < nf; ++x) out.light_begin[x + 1] += out.light_begin[x];

  // l.722: "the same sorting argument applies to the roots" - topological order of the root
  // quotient (Lemma 3.1).
  std::vector<std::uint32_t> root_index(nf, 0);
  for (std::size_t i = 0; i < f.roots.size(); ++i) {
    root_index[f.roots[i]] = static_cast<std::uint32_t>(i);
  }
  std::vector<DiArc> rq;
  rq.reserve(out.root_arcs.size());
  for (const ArcId a : out.root_arcs) {
    const auto& arc = f.graph.arcs[a];
    rq.push_back({root_index[q.component(f.leaf_of_state[arc.tail])],
                  root_index[q.component(f.leaf_of_state[arc.head])]});
  }
  std::vector<Vertex> topo;
  try {
    topo = topological_sort(Digraph::from_arcs(f.roots.size(), rq));
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument("normalize: the root quotient has a cycle (Def 4.1 H4)");
  }
  out.root_order.reserve(topo.size());
  for (const Vertex r : topo) out.root_order.push_back(f.roots[r]);

  bump(diag, "normalize_calls", 1);
  bump(diag, "normalize_input_nodes", n);
  bump(diag, "normalize_nodes", nf);
  bump(diag, "normalize_normalization_links", num_links);
  bump(diag, "normalize_link_sets", num_sets);
  bump(diag, "normalize_table_steps", table_steps);
  bump(diag, "normalize_contracted_nodes", n - nf);
  bump(diag, "normalize_self_arcs", self_arcs);
  bump(diag, "normalize_light_arcs", out.light_arcs.size());
  bump(diag, "normalize_root_arcs", out.root_arcs.size());
  bump(diag, "normalize_saturated_masses", saturated);
  return out;
}

template NormalizedHierarchy<std::uint64_t> normalize<std::uint64_t>(
    const OrderedHierarchy<std::uint64_t>&, Diagnostics*);
template NormalizedHierarchy<double> normalize<double>(const OrderedHierarchy<double>&,
                                                       Diagnostics*);

}  // namespace sssp::paper
