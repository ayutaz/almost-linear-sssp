// Property tests for the outer state graph of arXiv 2611.34911v1, Sec 8.1 (txt l.1966-2001) and
// Lemmas 8.1-8.3 (l.2002, l.2027, l.2053), on seeded graphs built through Alg 8 l.8-11 (owner
// decomposition, route names, G', OwnerDistances). Every reference is written in this file or
// in outer_states_test_support.hpp, in exact integer images (unsigned __int128).
//
// Checked statements:
//   * l.1967: exactly the radial edges with w > J_x + J_y + d_e are deleted, with J from an
//     independent Dijkstra on G'; owner tree edges survive (l.1969-1970); the classification
//     handed to Def 5.4 is the outer one after deletion, macro keys are the route names.
//   * Lemma 8.1 (first claim): every surviving radial edge has w <= 2 (J_x + J_y).
//   * l.1984-1997: every nonpivot has v_in and v_out; the inward tree states move to their
//     owner-tree parents along physical owner edges and end at the pivot's ordinary inward copy;
//     outward tree states are entered only from their parent tree state or the pivot's ordinary
//     outward copy; entries from every ordinary inward copy and the exit to the native outward
//     copy have length 0; macro arcs have ordinary tails and heads; special states are
//     singleton blocks (l.1998). Every arc projects to a physical walk of its length (H1).
//   * Lemma 8.1 (last claim): the macro-free block quotient is acyclic, with and without the
//     special states (Kahn's algorithm), and the reported block order is topological.
//   * Lemma 8.2: I0, O0, I*, O* (and B0, theta) equal brute-force values from their definitions
//     (Sec 5.3 l.1110-1117); B0 and theta are finite together; r <= theta and
//     B0 <= theta <= 4 K^2 h B0 with K = 16.
//   * Lemma 8.3: for every distinct arc or join weight W, W + one unit, and random W, the
//     components of the joins of weight <= W (< W) restricted to states are the SCCs of the
//     arcs of length <= W (< W) (Tarjan).
//   * Work, only under SSSP_DIAGNOSTICS: states <= 6 N, arcs <= 8 m + 7 N, special states =
//     2 (N - |A|), join records <= 2 S + |A|, one scan of the m edges (O(m + N), l.1979-1980).
//
// Inputs: the generic families of outer_states_test_support.hpp (N <= 20, so h <= 6 whenever
// there are two pivots), long owner paths (h up to 65), and the hand-built three-pivot corridor
// (three_pivot_corridor) on which the SCCs of the ordinary states change when v_in / v_out are
// removed. The generic and long-path inputs never produce such a threshold (J rises by at most
// a factor h along an owner path, and blocking an ordinary walk needs a rise across two
// shells); the coverage counter `special_needed` makes sure the Lemma 8.3 check is
// discriminating on the hand-built input (without the owner-tree arcs it fails there).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "paper/outer/outer_states.hpp"
#include "paper/owner/aux_metric.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/prng.hpp"
#include "unit/paper/outer/outer_states_test_support.hpp"

using sssp::Diagnostics;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::EdgeClass;
using sssp::paper::kNoState;
using sssp::paper::Phase;
using sssp::paper::RadialArcKind;
using sssp::paper::StateId;
using sssp::paper::outer::build_outer_states;
using sssp::paper::outer::OuterArcKind;
using sssp::paper::outer::OuterStates;
using sssp::paper::owner::OuterEdgeClass;
namespace to = sssp::test::outer;
using to::Exact;
using to::exact;
using to::exact_or_inf;
using to::kExactInf;
using to::kScale;

namespace {

struct Coverage {
  std::uint64_t cases = 0;
  std::uint64_t deleted = 0;
  std::uint64_t deleted_interowner = 0;
  std::uint64_t macro_cases = 0;
  std::uint64_t finite_theta = 0;
  std::uint64_t special_states = 0;
  std::uint64_t lower_inward = 0;
  std::uint64_t far_arcs = 0;
  std::uint64_t thresholds = 0;
  std::uint64_t scc_with_special = 0;  ///< SCCs of >= 2 states containing a special state
  std::uint64_t scc_multi_block = 0;   ///< SCCs spanning >= 2 blocks
  std::uint64_t theta_ratio_max = 0;   ///< max ceil(theta / B0) seen
  std::uint64_t special_needed = 0;    ///< thresholds where the specials change ordinary SCCs
  std::uint64_t max_hop_bound = 0;     ///< largest h with >= 2 pivots
};

/// l.1967, Lemma 8.1 (first claim), l.1969-1970, and the classification handed to Def 5.4.
template <class W>
void check_deletion(const to::OuterCase<W>& c, const OuterStates<W>& os,
                    const std::vector<Exact>& J, Coverage& cov) {
  const auto& st = os.states;
  for (Vertex v = 0; v < c.n; ++v) REQUIRE(exact(st.coordinate[v]) == kScale * J[v]);
  REQUIRE(os.deleted.size() == c.g.num_edges());
  std::uint64_t deleted = 0;
  for (EdgeId e = 0; e < c.g.num_edges(); ++e) {
    const auto& ed = c.g.edge(e);
    const OuterEdgeClass oc = c.cls.edge_class[e];
    const Exact w = exact(ed.weight);
    const bool radial = oc == OuterEdgeClass::kRadial;
    const bool expect = radial && w > J[ed.u] + J[ed.v] + exact(c.cls.detour[e]);
    INFO("edge " << e);
    CHECK((os.deleted[e] != 0) == expect);
    deleted += expect ? 1 : 0;
    if (expect && c.od.owner[ed.u] != c.od.owner[ed.v]) ++cov.deleted_interowner;
    const EdgeClass want = oc == OuterEdgeClass::kMacro               ? EdgeClass::kMacro
                           : oc == OuterEdgeClass::kRadial && !expect ? EdgeClass::kRadial
                                                                      : EdgeClass::kOmitted;
    CHECK(st.classification.edge_class[e] == want);
    CHECK(st.classification.macro_key[e] ==
          (want == EdgeClass::kMacro ? c.cls.name[e] : sssp::paper::kNoCorridor));
    if (radial && !expect) CHECK(w <= 2 * (J[ed.u] + J[ed.v]));  // Lemma 8.1
  }
  CHECK(os.num_deleted == deleted);
  cov.deleted += deleted;
  for (const EdgeId e : c.od.tree_edges) {
    if (c.od.cut[e] == 0) CHECK(os.deleted[e] == 0);  // J_v >= r_v >= w_e (l.1969-1970)
  }
}

/// l.1984-1998 and H1 for every arc.
template <class W>
void check_special_states(const to::OuterCase<W>& c, const OuterStates<W>& os, Coverage& cov) {
  const auto& st = os.states;
  const std::size_t s = st.num_states();
  const auto& arcs = st.graph.arcs;
  REQUIRE(st.projection.size() == s);
  REQUIRE(os.arc_kind.size() == arcs.size());
  REQUIRE(st.arc_kind.size() == arcs.size());
  REQUIRE(st.arc_edge.size() == arcs.size());
  // H1: an arc follows its physical edge with its weight, or stays at a vertex with length 0.
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    const Vertex x = st.projection[arcs[a].tail];
    const Vertex y = st.projection[arcs[a].head];
    if (st.arc_edge[a] == kNoEdge) {
      CHECK(x == y);
      CHECK(arcs[a].length == W{0});
    } else {
      const auto& ed = c.g.edge(st.arc_edge[a]);
      CHECK(((ed.u == x && ed.v == y) || (ed.u == y && ed.v == x)));
      CHECK(arcs[a].length == ed.weight);
    }
    if (os.arc_kind[a] == OuterArcKind::kOrdinary) {
      CHECK(a < os.num_ordinary_arcs);
      CHECK_FALSE(os.is_special(arcs[a].tail));  // macros and radial arcs: ordinary copies
      CHECK_FALSE(os.is_special(arcs[a].head));
      CHECK(st.arc_kind[a] != RadialArcKind::kOwnerTree);
      cov.far_arcs += st.arc_kind[a] == RadialArcKind::kFarTransfer ? 1 : 0;
    } else {
      CHECK(a >= os.num_ordinary_arcs);
      CHECK(st.arc_kind[a] == RadialArcKind::kOwnerTree);
    }
  }
  // Per vertex: the special states and the physical owner tree.
  std::size_t specials = 0;
  for (Vertex v = 0; v < c.n; ++v) {
    INFO("vertex " << v);
    if (st.native_inward[v] != kNoState && st.lower_inward[v] != kNoState) ++cov.lower_inward;
    if (c.od.is_pivot(v)) {
      CHECK(os.special_inward[v] == kNoState);
      CHECK(os.special_outward[v] == kNoState);
      CHECK(os.owner_parent[v] == kNoVertex);
      continue;
    }
    specials += 2;
    REQUIRE(os.special_inward[v] != kNoState);
    REQUIRE(os.special_outward[v] != kNoState);
    CHECK(os.is_special(os.special_inward[v]));
    CHECK(os.is_special(os.special_outward[v]));
    CHECK(st.projection[os.special_inward[v]] == v);
    CHECK(st.projection[os.special_outward[v]] == v);
    CHECK(st.phase[os.special_inward[v]] == Phase::kInward);
    CHECK(st.phase[os.special_outward[v]] == Phase::kOutward);
    // The parent is the next vertex towards p(v) on the owner tree T[v, p(v)].
    const Vertex p = os.owner_parent[v];
    const EdgeId pe = os.owner_parent_edge[v];
    REQUIRE(p != kNoVertex);
    REQUIRE(pe != kNoEdge);
    const auto& ed = c.g.edge(pe);
    CHECK(((ed.u == v && ed.v == p) || (ed.u == p && ed.v == v)));
    CHECK(c.od.cut[pe] == 0);
    CHECK(std::ranges::find(c.od.tree_edges, pe) != c.od.tree_edges.end());
    CHECK(c.od.owner[p] == c.od.owner[v]);
    Vertex x = v;
    std::uint32_t hops = 0;
    while (!c.od.is_pivot(x) && hops <= c.n) {
      x = os.owner_parent[x];
      ++hops;
    }
    CHECK(x == c.od.pivots[c.od.owner[v]]);
    CHECK(hops == c.od.owner_hops[v]);
  }
  CHECK(s == os.num_ordinary_states + specials);
  cov.special_states += specials;
  // The special arcs, by kind.
  const auto in_of = [&](Vertex v) {
    return c.od.is_pivot(v) ? st.native_inward[v] : os.special_inward[v];
  };
  const auto out_of = [&](Vertex v) {
    return c.od.is_pivot(v) ? st.native_outward[v] : os.special_outward[v];
  };
  std::vector<std::uint32_t> entries(c.n, 0), ups(c.n, 0), downs(c.n, 0), exits(c.n, 0);
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    const auto& arc = arcs[a];
    const Vertex v = st.projection[arc.tail];
    switch (os.arc_kind[a]) {
      case OuterArcKind::kOrdinary:
        break;
      case OuterArcKind::kTreeEntry:
        CHECK_FALSE(os.is_special(arc.tail));
        CHECK(st.phase[arc.tail] == Phase::kInward);
        CHECK(arc.head == os.special_inward[v]);
        ++entries[v];
        break;
      case OuterArcKind::kInwardTree:
        CHECK(arc.tail == os.special_inward[v]);
        CHECK(arc.head == in_of(os.owner_parent[v]));
        CHECK(st.arc_edge[a] == os.owner_parent_edge[v]);
        ++ups[v];
        break;
      case OuterArcKind::kOutwardTree: {
        const Vertex h = st.projection[arc.head];
        CHECK(arc.head == os.special_outward[h]);
        CHECK(arc.tail == out_of(os.owner_parent[h]));
        CHECK(st.arc_edge[a] == os.owner_parent_edge[h]);
        ++downs[h];
        break;
      }
      case OuterArcKind::kTreeExit:
        CHECK(arc.tail == os.special_outward[v]);
        CHECK(arc.head == st.native_outward[v]);
        ++exits[v];
        break;
    }
  }
  for (Vertex v = 0; v < c.n; ++v) {
    const bool piv = c.od.is_pivot(v);
    const std::uint32_t inward_copies = st.lower_inward[v] != kNoState ? 2 : 1;
    CHECK(entries[v] == (piv ? 0 : inward_copies));
    CHECK(ups[v] == (piv ? 0 : 1));
    CHECK(downs[v] == (piv ? 0 : 1));
    CHECK(exits[v] == (piv ? 0 : 1));
  }
  // Singleton blocks (l.1998).
  std::vector<std::uint32_t> block_size(st.num_blocks, 0);
  for (StateId x = 0; x < s; ++x) ++block_size[st.block[x]];
  for (std::size_t x = os.num_ordinary_states; x < s; ++x) CHECK(block_size[st.block[x]] == 1);
  CHECK(os.num_ordinary_blocks + specials == st.num_blocks);
}

/// Lemma 8.1 (last claim): acyclic macro-free block quotients.
template <class W>
void check_acyclic(const OuterStates<W>& os) {
  const auto& st = os.states;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> all, ordinary;
  for (std::size_t a = 0; a < st.graph.arcs.size(); ++a) {
    const auto k = st.arc_kind[a];
    if (k == RadialArcKind::kMacro || k == RadialArcKind::kWithinWindow) continue;
    const auto& arc = st.graph.arcs[a];
    const std::pair q{st.block[arc.tail], st.block[arc.head]};
    CHECK(q.first != q.second);  // only within-window arcs stay inside a block
    all.push_back(q);
    if (os.arc_kind[a] == OuterArcKind::kOrdinary) ordinary.push_back(q);
  }
  CHECK(to::is_acyclic(st.num_blocks, all));
  CHECK(to::is_acyclic(st.num_blocks, ordinary));
  const auto& order = os.joins.block_order;
  REQUIRE(order.size() == st.num_blocks);
  std::vector<std::uint32_t> rank(st.num_blocks, 0);
  for (std::uint32_t i = 0; i < order.size(); ++i) rank[order[i]] = i;
  for (const auto& [t, h] : all) CHECK(rank[t] < rank[h]);
}

/// I0 / O0 by definition: least macro weight whose head reaches z (z reaches a tail copy) over
/// admissible arcs. Macro arcs are taken in increasing weight; a state is labelled once.
std::vector<Exact> least_reaching(std::size_t s, const std::vector<to::XArc>& admissible,
                                  std::vector<std::pair<Exact, std::uint32_t>> seeds) {
  std::vector<std::vector<std::uint32_t>> adj(s);
  for (const auto& a : admissible) adj[a.tail].push_back(a.head);
  std::ranges::sort(seeds);
  std::vector<Exact> label(s, kExactInf);
  std::vector<std::uint32_t> queue;
  for (const auto& [w, x] : seeds) {
    if (label[x] != kExactInf) continue;
    label[x] = w;
    queue.assign(1, x);
    for (std::size_t i = 0; i < queue.size(); ++i) {
      for (const auto y : adj[queue[i]]) {
        if (label[y] != kExactInf) continue;
        label[y] = w;
        queue.push_back(y);
      }
    }
  }
  return label;
}

/// I* / O* by definition: min over seeds of the bottleneck (largest arc, seed weight included)
/// of an admissible walk from the seed to z (naive minimax Dijkstra).
std::vector<Exact> least_bottleneck(std::size_t s, const std::vector<to::XArc>& admissible,
                                    const std::vector<std::pair<Exact, std::uint32_t>>& seeds) {
  std::vector<std::vector<std::pair<std::uint32_t, Exact>>> adj(s);
  for (const auto& a : admissible) adj[a.tail].emplace_back(a.head, a.len);
  std::vector<Exact> d(s, kExactInf);
  for (const auto& [w, x] : seeds) d[x] = std::min(d[x], w);
  std::vector<std::uint8_t> done(s, 0);
  for (;;) {
    std::size_t u = s;
    for (std::size_t v = 0; v < s; ++v) {
      if (done[v] == 0 && d[v] != kExactInf && (u == s || d[v] < d[u])) u = v;
    }
    if (u == s) break;
    done[u] = 1;
    for (const auto& [y, len] : adj[u]) d[y] = std::min(d[y], std::max(d[u], len));
  }
  return d;
}

/// Lemma 8.2 with brute-force I0, O0, I*, O*.
template <class W>
void check_lemma_8_2(const to::OuterCase<W>& c, const OuterStates<W>& os, Coverage& cov) {
  const auto& st = os.states;
  const auto& jn = os.joins;
  const std::size_t s = st.num_states();
  REQUIRE(jn.num_states == s);
  std::vector<to::XArc> fwd, rev;
  std::vector<std::pair<Exact, std::uint32_t>> heads, tails;
  for (std::size_t a = 0; a < st.graph.arcs.size(); ++a) {
    const auto& arc = st.graph.arcs[a];
    const Exact len = exact(arc.length);
    if (st.arc_kind[a] == RadialArcKind::kMacro) {
      heads.emplace_back(len, arc.head);
      tails.emplace_back(len, arc.tail);
      continue;
    }
    fwd.push_back({arc.tail, arc.head, len});
    rev.push_back({arc.head, arc.tail, len});
  }
  if (!heads.empty()) ++cov.macro_cases;
  const auto i0 = least_reaching(s, fwd, heads);
  const auto o0 = least_reaching(s, rev, tails);
  const auto is = least_bottleneck(s, fwd, heads);
  const auto os_ = least_bottleneck(s, rev, tails);
  const Exact h = c.od.hop_bound;
  for (StateId z = 0; z < s; ++z) {
    INFO("state " << z << " (vertex " << st.projection[z] << ")");
    CHECK(exact_or_inf(jn.i0[z]) == i0[z]);
    CHECK(exact_or_inf(jn.o0[z]) == o0[z]);
    CHECK(exact_or_inf(jn.i_star[z]) == is[z]);
    CHECK(exact_or_inf(jn.o_star[z]) == os_[z]);
    const Exact b0 = std::max(i0[z], o0[z]);
    const Exact theta = std::max(is[z], os_[z]);
    CHECK(exact_or_inf(jn.b0[z]) == b0);
    CHECK(exact_or_inf(jn.theta[z]) == theta);
    CHECK((b0 == kExactInf) == (theta == kExactInf));  // finite together
    if (theta == kExactInf || b0 == kExactInf) continue;
    ++cov.finite_theta;
    CHECK(exact(c.od.radius[st.projection[z]]) <= theta);  // r <= theta
    CHECK(b0 <= theta);
    // theta <= 4 K^2 h B0 = 1024 h B0 (h saturates in Mode A: the bound then holds trivially).
    if (h < (Exact{1} << 64)) CHECK(theta <= Exact{1024} * h * b0);
    if (b0 > 0)
      cov.theta_ratio_max = std::max<std::uint64_t>(
          cov.theta_ratio_max, static_cast<std::uint64_t>((theta + b0 - 1) / b0));
  }
}

/// Lemma 8.3: join components = SCCs at <= W and < W.
template <class W>
void check_lemma_8_3(const OuterStates<W>& os, std::uint64_t seed, Coverage& cov) {
  const auto& st = os.states;
  const auto& jn = os.joins;
  const std::size_t s = st.num_states();
  std::vector<to::XArc> arcs;
  std::set<Exact> lengths;
  for (const auto& arc : st.graph.arcs) {
    arcs.push_back({arc.tail, arc.head, exact(arc.length)});
    lengths.insert(exact(arc.length));
  }
  for (const auto& j : jn.joins) {
    REQUIRE(j.a < jn.num_join_vertices());
    REQUIRE(j.b < jn.num_join_vertices());
    lengths.insert(exact(j.weight));
  }
  std::set<Exact> thresholds;
  for (const Exact x : lengths) {
    if (x > 0) thresholds.insert(x);
    thresholds.insert(x + 1);
  }
  sssp::test::SplitMix64 rng(seed);
  const Exact top = lengths.empty() ? 1 : *lengths.rbegin() + 2;
  for (int i = 0; i < 8; ++i) {
    thresholds.insert(1 + static_cast<Exact>(rng.below(static_cast<std::uint64_t>(top))));
  }
  for (const Exact t : thresholds) {
    for (const bool strict : {false, true}) {
      INFO("W = " << static_cast<std::uint64_t>(t) << (strict ? " (strict)" : " (inclusive)"));
      const auto keep = [&](Exact x) { return strict ? x < t : x <= t; };
      std::vector<std::pair<std::uint32_t, std::uint32_t>> sel;
      for (const auto& a : arcs) {
        if (keep(a.len)) sel.emplace_back(a.tail, a.head);
      }
      const auto scc = to::scc_ids(s, sel);
      to::UnionFind uf(jn.num_join_vertices());
      for (const auto& j : jn.joins) {
        if (keep(exact(j.weight))) uf.unite(j.a, j.b);
      }
      std::vector<std::uint32_t> comp(s);
      for (StateId z = 0; z < s; ++z) comp[z] = uf.find(z);
      CHECK(to::same_partition(scc, comp));
      ++cov.thresholds;
      // Coverage: thresholds at which the special states change the SCCs of the ordinary
      // states (the owner-tree states are needed because J is not monotone on the owner tree,
      // l.1984); without such thresholds the test would not notice their absence.
      std::vector<std::pair<std::uint32_t, std::uint32_t>> sel_ordinary;
      for (std::size_t a = 0; a < os.num_ordinary_arcs; ++a) {
        if (keep(arcs[a].len)) sel_ordinary.emplace_back(arcs[a].tail, arcs[a].head);
      }
      const auto scc_ordinary = to::scc_ids(os.num_ordinary_states, sel_ordinary);
      const std::vector<std::uint32_t> scc_restricted(
          scc.begin(), scc.begin() + static_cast<std::ptrdiff_t>(os.num_ordinary_states));
      if (!to::same_partition(scc_ordinary, scc_restricted)) ++cov.special_needed;
      // Coverage: nontrivial SCCs with special states or several blocks.
      std::vector<std::uint32_t> size(s, 0), first_block(s, 0xffffffffu);
      std::vector<std::uint8_t> special(s, 0), multi(s, 0);
      for (StateId z = 0; z < s; ++z) {
        const auto k = scc[z];
        ++size[k];
        if (z >= os.num_ordinary_states) special[k] = 1;
        if (first_block[k] == 0xffffffffu) first_block[k] = st.block[z];
        if (first_block[k] != st.block[z]) multi[k] = 1;
      }
      for (std::size_t k = 0; k < s; ++k) {
        if (size[k] >= 2 && special[k] != 0) ++cov.scc_with_special;
        if (multi[k] != 0) ++cov.scc_multi_block;
      }
    }
  }
}

template <class W>
std::uint64_t counter([[maybe_unused]] const Diagnostics& d, [[maybe_unused]] const char* name) {
  const auto it = d.counters.find(name);
  return it == d.counters.end() ? 0 : it->second;
}

/// Work bounds through the counters (SSSP_DIAGNOSTICS only).
template <class W>
void check_work([[maybe_unused]] const to::OuterCase<W>& c,
                [[maybe_unused]] const OuterStates<W>& os, [[maybe_unused]] const Diagnostics& d) {
#if SSSP_DIAGNOSTICS
  const std::uint64_t n = c.n;
  const std::uint64_t m = c.g.num_edges();
  const std::uint64_t pivots = c.od.pivots.size();
  CHECK(counter<W>(d, "outer_states.calls") == 1);
  CHECK(counter<W>(d, "outer_states.edges_scanned") == m);
  CHECK(counter<W>(d, "outer_states.deleted_radial_edges") == os.num_deleted);
  CHECK(counter<W>(d, "outer_states.special_states") == 2 * (n - pivots));
  CHECK(counter<W>(d, "outer_states.special_arcs") ==
        os.states.graph.arcs.size() - os.num_ordinary_arcs);
  const std::uint64_t states = counter<W>(d, "outer_states.states");
  const std::uint64_t arcs = counter<W>(d, "outer_states.arcs");
  CHECK(states == os.states.num_states());
  CHECK(arcs == os.states.graph.arcs.size());
  CHECK(states <= 6 * n);
  CHECK(arcs <= 8 * m + 7 * n);
  CHECK(counter<W>(d, "threshold_joins_records") <= 2 * states + pivots);
  CHECK(counter<W>(d, "radial_states_calls") == 1);
  CHECK(counter<W>(d, "threshold_joins_calls") == 1);
#endif
}

/// Every check of this file on one case.
template <class W>
void check_all(const to::OuterCase<W>& c, Coverage& cov) {
  Diagnostics diag;
  const auto os = build_outer_states(c.g, c.od, c.cls, c.aux, c.jr.scaled_distance, &diag);
  CHECK(diag.fallbacks.empty());  // Sec 8.1 never substitutes a search
  const auto J = to::reference_j(c);
  check_deletion(c, os, J, cov);
  check_special_states(c, os, cov);
  check_acyclic(os);
  check_lemma_8_2(c, os, cov);
  check_lemma_8_3(os, c.n * 131 + c.g.num_edges(), cov);
  check_work(c, os, diag);
  if (c.od.pivots.size() >= 2) {
    cov.max_hop_bound = std::max<std::uint64_t>(cov.max_hop_bound, c.od.hop_bound);
  }
  ++cov.cases;
}

/// Long owner paths: 2..4 chains of 60..140 light edges (weight 1..2) hanging from vertex 0 or
/// from an earlier chain, plus 4 chords per chain of weight 3 * 2^U[0, 12]. With separators
/// b = 24..100 the owner cells hold long tree paths, so h > 16 and the bound
/// theta <= 4 K^2 h B0 of Lemma 8.2 is exercised with larger h than in the generic families.
template <class W>
std::vector<sssp::Edge<W>> long_owner_paths(std::uint64_t seed, std::size_t& n) {
  namespace lt = sssp::test::local;
  sssp::test::SplitMix64 rng(seed);
  std::vector<sssp::Edge<W>> e;
  std::vector<std::vector<Vertex>> chains(2 + rng.below(3));
  Vertex next = 1;
  for (std::size_t c = 0; c < chains.size(); ++c) {
    const std::size_t len = rng.between(60, 140);
    for (std::size_t i = 0; i < len; ++i) chains[c].push_back(next++);
    Vertex top = 0;
    if (c > 0 && rng.below(2) == 0) {
      const auto& up = chains[rng.below(c)];
      top = up[rng.below(up.size())];
    }
    e.push_back({top, chains[c][0], lt::to_weight<W>(rng.between(1, 2), rng)});
    for (std::size_t i = 1; i < len; ++i) {
      e.push_back({chains[c][i - 1], chains[c][i], lt::to_weight<W>(rng.between(1, 2), rng)});
    }
  }
  n = next;
  const auto any = [&]() { return static_cast<Vertex>(rng.below(n)); };
  for (std::size_t i = 0; i < chains.size() * 4; ++i) {
    e.push_back({any(), any(), lt::to_weight<W>(std::uint64_t{3} << rng.between(0, 12), rng)});
  }
  return e;
}

}  // namespace

TEMPLATE_TEST_CASE("Sec 8.1 / Lemmas 8.1-8.3: outer states on seeded graphs",
                   "[property][outer_states]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  to::for_each_case<W>({1, 2, 3, 5, 8, 13, 20}, 2, [&](to::OuterCase<W>& c) { check_all(c, cov); });
  WARN("cases " << cov.cases << ", deleted " << cov.deleted << " (interowner "
                << cov.deleted_interowner << "), macro cases " << cov.macro_cases
                << ", finite theta " << cov.finite_theta << ", lower inward copies "
                << cov.lower_inward << ", far arcs " << cov.far_arcs << ", thresholds "
                << cov.thresholds << ", SCCs with special states " << cov.scc_with_special
                << ", multi-block SCCs " << cov.scc_multi_block << ", max theta/B0 "
                << cov.theta_ratio_max << ", thresholds where the specials matter "
                << cov.special_needed << ", max h with >= 2 pivots " << cov.max_hop_bound);
  CHECK(cov.cases > 0);
  CHECK(cov.deleted > 0);
  CHECK(cov.deleted_interowner > 0);
  CHECK(cov.macro_cases > 0);
  CHECK(cov.finite_theta > 0);
  CHECK(cov.special_states > 0);
  CHECK(cov.lower_inward > 0);
  CHECK(cov.far_arcs > 0);
  CHECK(cov.scc_with_special > 0);
  CHECK(cov.scc_multi_block > 0);
}

TEMPLATE_TEST_CASE("Lemmas 8.1-8.3: long owner paths (h > 16) and the three-pivot corridor",
                   "[property][outer_states]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  for (std::uint64_t seed = 1; seed <= 2; ++seed) {
    std::size_t n = 0;
    const auto edges = long_owner_paths<W>(seed * 15485863, n);
    for (const std::size_t b : {std::size_t{24}, std::size_t{64}}) {
      for (const auto names : {to::Names::kExactP1, to::Names::kNameRoutesPaper}) {
        INFO("long-owner-paths seed="
             << seed << " n=" << n << " b=" << b
             << (names == to::Names::kExactP1 ? " exact-P1" : " name-routes"));
        auto c = to::make_outer_case<W>(n, edges, b, names, seed + b);
        check_all(c, cov);
      }
    }
  }
  const std::uint64_t needed_before = cov.special_needed;
  // The hand-built input of outer_states_test_support.hpp on which v_in / v_out decide the SCCs.
  for (const std::size_t k : {std::size_t{34}, std::size_t{40}, std::size_t{57}}) {
    for (const auto names : {to::Names::kExactP1, to::Names::kNameRoutesPaper}) {
      INFO("three-pivot-corridor k="
           << k << (names == to::Names::kExactP1 ? " exact-P1" : " name-routes"));
      to::ThreePivotCorridor r;
      const auto edges = to::three_pivot_corridor<W>(k, r);
      auto c = to::make_outer_case<W>(r.n, edges, r.b, names, k);
      check_all(c, cov);
    }
  }
  WARN("cases " << cov.cases << ", deleted " << cov.deleted << ", finite theta " << cov.finite_theta
                << ", thresholds " << cov.thresholds << ", thresholds where the specials matter "
                << cov.special_needed << " (random " << needed_before << "), max theta/B0 "
                << cov.theta_ratio_max << ", max h with >= 2 pivots " << cov.max_hop_bound);
  CHECK(cov.cases == 14);
  CHECK(cov.max_hop_bound > 16);
  CHECK(cov.finite_theta > 0);
  CHECK(cov.special_needed > needed_before);
}
