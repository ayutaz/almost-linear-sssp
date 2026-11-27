// Threshold components and joins, Sec 5.3 (Lemmas 5.8-5.9). See threshold_joins.hpp.

#include "paper/local/threshold_joins.hpp"

#include <algorithm>
#include <span>
#include <stdexcept>
#include <utility>

#include "paper/primitives/mst.hpp"
#include "paper/primitives/search.hpp"

namespace sssp::paper {

namespace {

template <class W>
constexpr W kInf = WeightTraits<W>::infinity();

/// A finite (weight, attaining edge) or +infinity with kNoEdge; ordered by (weight, edge id).
template <class W>
struct Attained {
  W w = kInf<W>;
  EdgeId e = kNoEdge;

  void relax(const Attained& o) {
    if (MstEdgeKey<W>{o.w, o.e} < MstEdgeKey<W>{w, e}) *this = o;
  }
};

template <class W>
bool is_eligible(RadialArcKind k) {
  return k != RadialArcKind::kMacro;  // radial arcs and switches (l.1110-1111)
}

}  // namespace

template <class W>
ThresholdJoins<W> compute_joins(const OwnerDecomposition<W>& od, const RadialStateGraph<W>& rs,
                                Diagnostics* diag) {
  if (rs.num_vertices() != od.num_vertices() || rs.projection.size() != rs.num_states()) {
    throw std::invalid_argument("compute_joins: state graph of a different decomposition");
  }
  const std::size_t s = rs.num_states();
  const std::size_t nb = rs.num_blocks;
  const auto& arcs = rs.graph.arcs;
  ThresholdJoins<W> out;
  out.num_states = s;
  out.num_pivots = od.pivots.size();

  // Members of each block (stable counting sort by block id, Lemma 3.1).
  std::vector<std::size_t> mbegin(nb + 1, 0);
  for (StateId x = 0; x < s; ++x) ++mbegin[rs.block[x] + 1];
  for (std::size_t b = 0; b < nb; ++b) mbegin[b + 1] += mbegin[b];
  std::vector<StateId> members(s);
  {
    std::vector<std::size_t> pos(mbegin.begin(), mbegin.end() - 1);
    for (StateId x = 0; x < s; ++x) members[pos[rs.block[x]]++] = x;
  }

  // Block MSTs (Lemma 3.5) on the within-window arcs, one orientation per physical edge, and a
  // rooted order of each tree for the two sweeps (l.1117-1121).
  std::vector<std::size_t> wbegin(nb + 1, 0);
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    if (rs.arc_kind[a] == RadialArcKind::kWithinWindow && arcs[a].tail < arcs[a].head) {
      ++wbegin[rs.block[arcs[a].tail] + 1];
    }
  }
  for (std::size_t b = 0; b < nb; ++b) wbegin[b + 1] += wbegin[b];
  std::vector<ArcId> within(wbegin[nb]);
  {
    std::vector<std::size_t> pos(wbegin.begin(), wbegin.end() - 1);
    for (std::size_t a = 0; a < arcs.size(); ++a) {
      if (rs.arc_kind[a] == RadialArcKind::kWithinWindow && arcs[a].tail < arcs[a].head) {
        within[pos[rs.block[arcs[a].tail]]++] = static_cast<ArcId>(a);
      }
    }
  }
  std::vector<StateId> tree_order(s);  // per block, a BFS order of its MST (block-major)
  std::vector<StateId> tree_parent(s, kNoState);
  std::vector<W> tree_w(s, W{0});
  std::vector<std::uint32_t> local(s, 0);
  std::uint64_t block_msts = 0;
  for (std::size_t b = 0; b < nb; ++b) {
    const std::size_t k = mbegin[b + 1] - mbegin[b];
    if (k == 1) {
      tree_order[mbegin[b]] = members[mbegin[b]];
      continue;
    }
    for (std::size_t i = 0; i < k; ++i)
      local[members[mbegin[b] + i]] = static_cast<std::uint32_t>(i);
    std::vector<Edge<W>> local_edges;
    local_edges.reserve(wbegin[b + 1] - wbegin[b]);
    for (std::size_t t = wbegin[b]; t < wbegin[b + 1]; ++t) {
      const auto& a = arcs[within[t]];
      local_edges.push_back({local[a.tail], local[a.head], a.length});
    }
    const auto mst = mst_three_pass<W>(k, std::span<const Edge<W>>(local_edges));
    ++block_msts;
    std::vector<std::vector<std::pair<std::uint32_t, std::size_t>>> adj(k);
    for (const EdgeId le : mst.edges) {
      const ArcId a = within[wbegin[b] + le];
      out.joins.push_back({arcs[a].tail, arcs[a].head, arcs[a].length, JoinKind::kBlockMst, a});
      adj[local_edges[le].u].emplace_back(local_edges[le].v, le);
      adj[local_edges[le].v].emplace_back(local_edges[le].u, le);
    }
    std::vector<std::uint8_t> seen(k, 0);
    std::size_t tail = mbegin[b];
    tree_order[tail++] = members[mbegin[b]];
    seen[0] = 1;
    for (std::size_t i = mbegin[b]; i < tail; ++i) {
      const StateId x = tree_order[i];
      for (const auto& [ly, le] : adj[local[x]]) {
        if (seen[ly] != 0) continue;
        seen[ly] = 1;
        const StateId y = members[mbegin[b] + ly];
        tree_parent[y] = x;
        tree_w[y] = local_edges[le].weight;
        tree_order[tail++] = y;
      }
    }
    if (tail != mbegin[b + 1]) throw std::logic_error("compute_joins: block MST not spanning");
  }

  // The macro-free block quotient and its topological order (Lemma 5.7, Lemma 3.1).
  std::vector<DiArc> qarcs;
  std::vector<ArcId> qarc_id;
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    if (!is_eligible<W>(rs.arc_kind[a]) || rs.arc_kind[a] == RadialArcKind::kWithinWindow) continue;
    qarcs.push_back({rs.block[arcs[a].tail], rs.block[arcs[a].head]});
    qarc_id.push_back(static_cast<ArcId>(a));
  }
  const Digraph quotient = Digraph::from_arcs(nb, qarcs);
  try {
    out.block_order = topological_sort(quotient);
  } catch (const std::invalid_argument&) {
    throw std::logic_error("compute_joins: the macro-free block quotient has a cycle");
  }
  const auto out_arcs = [&](std::uint32_t b) { return quotient.arc_ids(b); };

  // I0, O0 (l.1112-1115): seeds at macro head resets / tail copies; min propagation forward and
  // backward over the blocks, ignoring radial weights.
  std::vector<Attained<W>> bi(nb), bo(nb);
  std::vector<W> post_in(s, kInf<W>), post_out(s, kInf<W>);
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    if (rs.arc_kind[a] != RadialArcKind::kMacro) continue;
    const Attained<W> seed{arcs[a].length, rs.arc_edge[a]};
    bi[rs.block[arcs[a].head]].relax(seed);
    bo[rs.block[arcs[a].tail]].relax(seed);
    post_in[arcs[a].head] = std::min(post_in[arcs[a].head], arcs[a].length);
    post_out[arcs[a].tail] = std::min(post_out[arcs[a].tail], arcs[a].length);
  }
  for (const std::uint32_t b : out.block_order) {
    for (const EdgeId q : out_arcs(b)) bi[qarcs[q].head].relax(bi[b]);
  }
  for (std::size_t i = nb; i-- > 0;) {
    const std::uint32_t b = out.block_order[i];
    for (const EdgeId q : out_arcs(b)) bo[b].relax(bo[qarcs[q].head]);
  }

  // Multisource bottleneck on a block MST (l.1117-1121): postorder takes the minimum of the own
  // seed and the child offers raised to the edge weight; preorder sends the parent's value back.
  std::vector<W> up(s, kInf<W>);
  const auto sweep = [&](std::uint32_t b, const std::vector<W>& cost, std::vector<W>& res) {
    const std::size_t lo = mbegin[b];
    const std::size_t hi = mbegin[b + 1];
    for (std::size_t i = lo; i < hi; ++i) up[tree_order[i]] = cost[tree_order[i]];
    for (std::size_t i = hi; i-- > lo + 1;) {
      const StateId v = tree_order[i];
      const StateId p = tree_parent[v];
      up[p] = std::min(up[p], std::max(up[v], tree_w[v]));
    }
    res[tree_order[lo]] = up[tree_order[lo]];
    for (std::size_t i = lo + 1; i < hi; ++i) {
      const StateId v = tree_order[i];
      res[v] = std::min(up[v], std::max(res[tree_parent[v]], tree_w[v]));
    }
  };

  // I* in block order, relaxing each outgoing transfer by max(value, length) (l.1116-1122).
  out.i_star.assign(s, kInf<W>);
  for (const std::uint32_t b : out.block_order) {
    sweep(b, post_in, out.i_star);
    for (const EdgeId q : out_arcs(b)) {
      const auto& a = arcs[qarc_id[q]];
      post_in[a.head] = std::min(post_in[a.head], std::max(out.i_star[a.tail], a.length));
    }
  }
  // O*: all directions reversed (l.1122).
  out.o_star.assign(s, kInf<W>);
  for (std::size_t i = nb; i-- > 0;) {
    const std::uint32_t b = out.block_order[i];
    for (const EdgeId q : out_arcs(b)) {
      const auto& a = arcs[qarc_id[q]];
      post_out[a.tail] = std::min(post_out[a.tail], std::max(out.o_star[a.head], a.length));
    }
    sweep(b, post_out, out.o_star);
  }

  // Per-state values and the attaining macro edges.
  out.i0.resize(s);
  out.o0.resize(s);
  out.b0.resize(s);
  out.theta.resize(s);
  out.i0_edge.resize(s);
  out.o0_edge.resize(s);
  out.b0_edge.resize(s);
  out.b0_key.assign(s, kNoCorridor);
  std::uint64_t finite_theta = 0;
  for (StateId z = 0; z < s; ++z) {
    const auto& ai = bi[rs.block[z]];
    const auto& ao = bo[rs.block[z]];
    out.i0[z] = ai.w;
    out.o0[z] = ao.w;
    out.i0_edge[z] = ai.e;
    out.o0_edge[z] = ao.e;
    const bool from_in = !(ai.w < ao.w);  // [Interpretation] I0's edge on ties
    out.b0[z] = from_in ? ai.w : ao.w;
    out.b0_edge[z] = from_in ? ai.e : ao.e;
    if (out.b0_edge[z] != kNoEdge) out.b0_key[z] = rs.classification.macro_key[out.b0_edge[z]];
    out.theta[z] = std::max(out.i_star[z], out.o_star[z]);
    if (out.theta[z] != kInf<W>) ++finite_theta;
  }

  // Join records (l.1124-1128): block MST edges (above), owner joins, corridor joins.
  for (StateId z = 0; z < s; ++z) {
    if (out.theta[z] == kInf<W>) continue;
    const std::uint32_t p = od.owner[rs.projection[z]];
    out.joins.push_back({z, out.aux_vertex(p), out.theta[z], JoinKind::kOwner, p});
  }
  for (std::uint32_t c = 0; c < od.corridors.size(); ++c) {
    const auto& cr = od.corridors[c];
    out.joins.push_back(
        {out.aux_vertex(cr.lower), out.aux_vertex(cr.upper), cr.weight, JoinKind::kCorridor, c});
  }

  static_cast<void>(block_msts);
  static_cast<void>(finite_theta);
  if (diag != nullptr) {
    SSSP_COUNT(++diag->counters["threshold_joins_calls"]);
    SSSP_COUNT(diag->counters["threshold_joins_records"] += out.joins.size());
    SSSP_COUNT(diag->counters["threshold_joins_finite_theta"] += finite_theta);
    SSSP_COUNT(diag->counters["threshold_joins_block_msts"] += block_msts);
  }
  return out;
}

template ThresholdJoins<std::uint64_t> compute_joins<std::uint64_t>(
    const OwnerDecomposition<std::uint64_t>&, const RadialStateGraph<std::uint64_t>&, Diagnostics*);
template ThresholdJoins<double> compute_joins<double>(const OwnerDecomposition<double>&,
                                                      const RadialStateGraph<double>&,
                                                      Diagnostics*);

}  // namespace sssp::paper
