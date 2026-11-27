// Property tests for the numerical scales of the outer construction: arXiv 2611.34911v1, Sec
// 8.2 (txt l.2076-2104) and Lemma 8.4 (l.2105-2142), on seeded graphs (SC-11, SC-12 of
// docs/research/04_outer_scales_pipeline.md Sec 6.1).
//
// Every claim is checked twice: by check_outer_scales (src) and by the assertions written here
// from the paper's statements. The merge plan is checked against Tarjan's SCCs of the state
// arcs (Lemma 8.3, strict form): the states merged by a class [g, g') are one SCC of the arcs
// shorter than g', and the parts they had before are the SCCs of the arcs shorter than g.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/atom_replacement.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "paper/local/threshold_joins.hpp"
#include "paper/outer/outer_states.hpp"
#include "paper/outer/scales.hpp"
#include "paper/outer/scales_numeric.hpp"
#include "paper/owner/aux_metric.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "support/prng.hpp"
#include "unit/paper/outer/outer_scales_test_support.hpp"
#include "unit/paper/outer/outer_states_test_support.hpp"

using sssp::Edge;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::JoinKind;
using sssp::paper::kNoAtom;
using sssp::paper::StateId;
using sssp::paper::outer::build_outer_merge_plan;
using sssp::paper::outer::check_outer_scales;
using sssp::paper::outer::OuterBlockKind;
using sssp::paper::outer::OuterEventGrid;
using sssp::paper::outer::OuterMergePlan;
using sssp::test::SplitMix64;
namespace to = sssp::test::outer;

namespace {

template <class W>
constexpr W kInf = WeightTraits<W>::infinity();

struct Coverage {
  std::uint64_t cases = 0;
  std::uint64_t aligned = 0;
  std::uint64_t unaligned = 0;
  std::uint64_t unattached = 0;
  std::uint64_t atoms = 0;
  std::uint64_t unaligned_atoms = 0;
  std::uint64_t aligned_atoms = 0;
  std::uint64_t max_atom_states = 0;
  std::uint64_t separate_events = 0;
  std::uint64_t global_events = 0;
  std::uint64_t lowered_supports = 0;
  std::uint64_t classes = 0;
  std::uint64_t scc_checks = 0;
  std::uint64_t merges = 0;  ///< class groups with at least two parts
  std::uint64_t atom_base_cases = 0;
  std::uint64_t atom_radial = 0;
  std::uint64_t distance_checks = 0;
};

/// x + x without overflow concerns for the comparisons below (saturating).
template <class W>
W twice(W x) {
  if constexpr (std::is_integral_v<W>) {
    return x > kInf<W> / 2 ? kInf<W> : static_cast<W>(x + x);
  } else {
    return x + x;
  }
}

/// The statements of Sec 8.2 / Lemma 8.4 on the scales, written from the paper.
template <class W>
void check_scales(const to::ScalesCase<W>& r, Coverage& cov) {
  const auto& st = r.os.states;
  const auto& sc = r.sc;
  const auto& od = r.c.od;
  const auto& jr = r.os.joins.joins;
  const auto& theta = r.os.joins.theta;
  {
    const auto bad = check_outer_scales(r.c.g, od, r.c.cls, r.c.aux, r.c.jr.calendars, r.os, sc);
    for (const auto& msg : bad) UNSCOPED_INFO(msg);
    REQUIRE(bad.empty());
  }
  REQUIRE(r.diag.fallbacks.empty());
  REQUIRE(sc.blocks.size() == st.num_blocks);
  // s_U = least original radius (Q-35) and s_U / 2 < t_U <= s_U (Q-36).
  std::vector<W> least(st.num_blocks, kInf<W>);
  std::vector<std::size_t> size(st.num_blocks, 0);
  for (StateId x = 0; x < st.num_states(); ++x) {
    ++size[st.block[x]];
    least[st.block[x]] = std::min(least[st.block[x]], od.radius[st.projection[x]]);
  }
  for (std::size_t u = 0; u < st.num_blocks; ++u) {
    const auto& b = sc.blocks[u];
    REQUIRE((b.kind == OuterBlockKind::kSingleton) == (size[u] == 1));
    if (size[u] == 1) continue;
    REQUIRE(b.s == least[u]);
    REQUIRE(b.t <= b.s);
    REQUIRE(b.s < twice(b.t));
    cov.aligned += b.kind == OuterBlockKind::kAligned ? 1 : 0;
    cov.unaligned += b.kind == OuterBlockKind::kUnaligned ? 1 : 0;
    cov.unattached += b.kind == OuterBlockKind::kUnattached ? 1 : 0;
    if (b.kind == OuterBlockKind::kAligned && b.s < sc.key_weight[b.key]) ++cov.lowered_supports;
  }
  for (std::size_t j = 0; j < jr.size(); ++j) {
    const auto& ev = sc.events[j];
    cov.separate_events += ev.grid == OuterEventGrid::kSeparate ? 1 : 0;
    cov.global_events += ev.grid == OuterEventGrid::kGlobal ? 1 : 0;
    if (jr[j].kind != JoinKind::kBlockMst) {
      // Finite joins are classified on the global grid, in a class that brackets them.
      REQUIRE(ev.grid == OuterEventGrid::kGlobal);
      REQUIRE(ev.lower <= jr[j].weight);
      REQUIRE(jr[j].weight < ev.upper);
      continue;
    }
    const auto& b = sc.blocks[st.block[jr[j].a]];
    const W w = jr[j].weight;
    // Subthreshold edges are local; significant ones are classified at or above t_U.
    REQUIRE((w < b.t) == (ev.grid == OuterEventGrid::kAtom));
    if (ev.grid != OuterEventGrid::kAtom) {
      REQUIRE(b.t <= ev.lower);
      REQUIRE(ev.lower <= w);
      REQUIRE(w < ev.upper);
    }
    if (b.kind == OuterBlockKind::kUnaligned) {
      // l.2086-2087: below a_U / 2, hence below half the earliest join of the block.
      REQUIRE(twice(w) < sc.name_weight[b.name]);
      REQUIRE(twice(w) < theta[b.z]);
    }
    if (b.kind == OuterBlockKind::kAligned && !(w < b.t)) {
      // Lemma 8.4 item 1: significant aligned weights are covered by the name's index range.
      const auto& sp = sc.grid.supports[b.key];
      REQUIRE(sc.grid.points[sp.lower_index] <= w);
      REQUIRE(w < sc.grid.points[sp.upper_index]);
    }
  }
  // Lemma 8.4: each subthreshold atom has fewer than b vertices, one owner and one radius.
  for (std::size_t y = 0; y < sc.atoms.size(); ++y) {
    const auto& atom = sc.atoms[y];
    REQUIRE(atom.states.size() >= 2);
    REQUIRE(atom.states.size() < od.separator);
    const Vertex v0 = st.projection[atom.states[0]];
    for (const StateId x : atom.states) {
      REQUIRE(od.owner[st.projection[x]] == od.owner[v0]);
      REQUIRE(od.radius[st.projection[x]] == od.radius[v0]);
      REQUIRE_FALSE(od.is_pivot(st.projection[x]));
    }
    REQUIRE(atom.cutoff <= od.radius[v0]);
    for (const auto e : atom.spanning_tree) REQUIRE(atom.edges[e].weight < atom.cutoff);
    ++cov.atoms;
    cov.max_atom_states = std::max<std::uint64_t>(cov.max_atom_states, atom.states.size());
    const auto kind = sc.blocks[sc.atom_block[y]].kind;
    cov.unaligned_atoms += kind == OuterBlockKind::kUnaligned ? 1 : 0;
    cov.aligned_atoms += kind == OuterBlockKind::kAligned ? 1 : 0;
  }
  ++cov.cases;
}

/// The merge plan: parts, outside arcs, and the classes in processing order. With `sccs`, the
/// SCC statements of note 04 Sec 2.9 against Tarjan (SC-12).
template <class W>
void check_plan(const to::ScalesCase<W>& r, const OuterMergePlan<W>& plan, bool sccs,
                Coverage& cov) {
  const auto& st = r.os.states;
  const auto& sc = r.sc;
  const auto& jr = r.os.joins.joins;
  const std::size_t ns = st.num_states();
  const std::size_t np = r.os.joins.num_pivots;

  // Parts: one per state outside the atoms, one per atom (Def 5.12, Q-74).
  REQUIRE(plan.part.size() == ns);
  REQUIRE(plan.part_atom.size() == plan.num_parts);
  REQUIRE(plan.part_state.size() == plan.num_parts);
  std::vector<std::size_t> part_size(plan.num_parts, 0);
  for (StateId x = 0; x < ns; ++x) {
    const std::uint32_t p = plan.part[x];
    REQUIRE(p < plan.num_parts);
    ++part_size[p];
    REQUIRE(plan.part_atom[p] == sc.state_atom[x]);
    if (sc.state_atom[x] == kNoAtom) REQUIRE(plan.part_state[p] == x);
  }
  std::size_t atom_parts = 0;
  for (std::uint32_t p = 0; p < plan.num_parts; ++p) {
    if (plan.part_atom[p] == kNoAtom) {
      REQUIRE(part_size[p] == 1);
    } else {
      ++atom_parts;
      REQUIRE(part_size[p] == sc.atoms[plan.part_atom[p]].states.size());
      REQUIRE(plan.part_state[p] == sssp::paper::kNoState);
    }
  }
  REQUIRE(atom_parts == sc.atoms.size());

  // D: exactly the arcs that do not lie inside an atom, in their order.
  REQUIRE(plan.graph.num_states == ns);
  REQUIRE(plan.arc_of.size() == plan.graph.arcs.size());
  std::size_t next = 0;
  for (std::uint32_t a = 0; a < st.graph.arcs.size(); ++a) {
    const auto& arc = st.graph.arcs[a];
    const bool inside =
        sc.state_atom[arc.tail] != kNoAtom && sc.state_atom[arc.tail] == sc.state_atom[arc.head];
    if (inside) continue;
    REQUIRE(next < plan.arc_of.size());
    REQUIRE(plan.arc_of[next] == a);
    REQUIRE(plan.graph.arcs[next].tail == arc.tail);
    REQUIRE(plan.graph.arcs[next].head == arc.head);
    REQUIRE(plan.graph.arcs[next].length == arc.length);
    ++next;
  }
  REQUIRE(next == plan.arc_of.size());

  // Classes: every join outside the atoms exactly once, with the class of its event; separate
  // grids first, then the global classes in increasing index (Q-68).
  std::vector<std::uint8_t> seen(jr.size(), 0);
  bool global_started = false;
  std::size_t last_global = 0;
  to::UnionFind uf(ns + np);
  std::vector<W> height(ns + np, W{0});
  for (std::size_t j = 0; j < jr.size(); ++j) {
    if (sc.events[j].grid == OuterEventGrid::kAtom) uf.unite(jr[j].a, jr[j].b);
  }
  for (const auto& atom : sc.atoms) height[uf.find(atom.states[0])] = atom.cutoff;

  std::vector<std::pair<std::uint32_t, std::uint32_t>> arcs_lower, arcs_upper;
  for (const auto& cls : plan.classes) {
    ++cov.classes;
    REQUIRE_FALSE(cls.joins.empty());
    REQUIRE(W{0} < cls.height);
    REQUIRE(cls.height < cls.upper);
    if (cls.grid == OuterEventGrid::kGlobal) {
      REQUIRE((!global_started || last_global < cls.index));
      global_started = true;
      last_global = cls.index;
    } else {
      REQUIRE(cls.grid == OuterEventGrid::kSeparate);
      REQUIRE_FALSE(global_started);
    }
    std::vector<std::uint32_t> before;
    if (sccs) {
      before.resize(ns);
      for (StateId x = 0; x < ns; ++x) before[x] = uf.find(x);
    }
    for (const std::uint32_t j : cls.joins) {
      REQUIRE(j < jr.size());
      REQUIRE(seen[j] == 0);
      seen[j] = 1;
      REQUIRE(sc.events[j].grid == cls.grid);
      REQUIRE(sc.events[j].lower == cls.height);
      REQUIRE(sc.events[j].upper == cls.upper);
      if (cls.grid == OuterEventGrid::kSeparate) REQUIRE(st.block[jr[j].a] == cls.block);
      // H2 and Def 5.12: heights do not decrease towards the root, and an atom's cutoff is at
      // most its parent's height.
      REQUIRE(height[uf.find(jr[j].a)] <= cls.height);
      REQUIRE(height[uf.find(jr[j].b)] <= cls.height);
    }
    for (const std::uint32_t j : cls.joins) uf.unite(jr[j].a, jr[j].b);
    std::vector<std::uint32_t> merged;
    for (const std::uint32_t j : cls.joins) {
      const auto root = uf.find(jr[j].a);
      height[root] = cls.height;
      merged.push_back(root);
    }
    std::ranges::sort(merged);
    merged.erase(std::unique(merged.begin(), merged.end()), merged.end());
    if (!sccs) continue;
    // Lemma 8.3 (strict): arcs shorter than the class upper end / shorter than the height.
    arcs_lower.clear();
    arcs_upper.clear();
    for (const auto& arc : st.graph.arcs) {
      if (arc.length < cls.upper) arcs_upper.emplace_back(arc.tail, arc.head);
      if (arc.length < cls.height) arcs_lower.emplace_back(arc.tail, arc.head);
    }
    const auto scc_upper = to::scc_ids(ns, arcs_upper);
    const auto scc_lower = to::scc_ids(ns, arcs_lower);
    std::vector<std::uint32_t> group_now, group_before, lower_ids, upper_ids;
    for (StateId x = 0; x < ns; ++x) {
      const auto root = uf.find(x);
      if (!std::ranges::binary_search(merged, root)) {
        // A state outside the merged groups is in none of their SCCs (checked through the
        // partition equality on the union below: its SCC id must not occur there).
        continue;
      }
      group_now.push_back(root);
      group_before.push_back(before[x]);
      lower_ids.push_back(scc_lower[x]);
      upper_ids.push_back(scc_upper[x]);
    }
    // The states below each new node are one SCC of the arcs < upper ...
    REQUIRE(to::same_partition(group_now, upper_ids));
    std::vector<std::uint32_t> inside_ids = upper_ids;
    std::ranges::sort(inside_ids);
    for (StateId x = 0; x < ns; ++x) {
      if (!std::ranges::binary_search(merged, uf.find(x))) {
        REQUIRE_FALSE(std::ranges::binary_search(inside_ids, scc_upper[x]));
      }
    }
    // ... and its children are the SCCs of the arcs < height among them.
    REQUIRE(to::same_partition(group_before, lower_ids));
    std::ranges::sort(group_before);
    group_before.erase(std::unique(group_before.begin(), group_before.end()), group_before.end());
    std::ranges::sort(group_now);
    group_now.erase(std::unique(group_now.begin(), group_now.end()), group_now.end());
    cov.merges += group_before.size() > group_now.size() ? 1 : 0;
    ++cov.scc_checks;
  }
  for (std::size_t j = 0; j < jr.size(); ++j) {
    REQUIRE((seen[j] != 0) == (sc.events[j].grid != OuterEventGrid::kAtom));
  }
  if (sccs) {
    // After the last class the groups are the SCCs of the whole state graph (the roots).
    std::vector<std::pair<std::uint32_t, std::uint32_t>> all;
    for (const auto& arc : st.graph.arcs) all.emplace_back(arc.tail, arc.head);
    std::vector<std::uint32_t> final_group(ns);
    for (StateId x = 0; x < ns; ++x) final_group[x] = uf.find(x);
    REQUIRE(to::same_partition(final_group, to::scc_ids(ns, all)));
  }
}

/// S8 and S9: the partition hierarchy (Def 5.12), its substitution (Lemma 5.15), Def 4.1 with
/// c0 <= 2 on the result, and the distances of Evaluate seeded as in Alg 8 l.14 against
/// Bellman-Ford (SC-14).
template <class W>
void check_partition(const to::ScalesCase<W>& r, const sssp::paper::LocalHierarchyParams& local,
                     bool distances, Coverage& cov) {
  using sssp::paper::kNoNode;
  using sssp::paper::NodeId;
  const auto& st = r.os.states;
  const auto& sc = r.sc;
  sssp::Diagnostics diag;
  const auto part = sssp::paper::outer::build_outer_partition(r.os, sc, &diag);
  const auto& ph = part.hierarchy;
  REQUIRE(diag.fallbacks.empty());  // S8 never substitutes a search
  REQUIRE(ph.c0 == 2);
  REQUIRE(ph.graph.num_states == st.num_states());
  REQUIRE(ph.projection == st.projection);
  REQUIRE(ph.atoms.size() == sc.atoms.size());
  REQUIRE(ph.leaf_atom.size() == ph.nodes.size());
  REQUIRE(part.node_upper.size() == ph.nodes.size());
  REQUIRE(part.leaf_of_state.size() == st.num_states());
  std::size_t roots = 0;
  for (NodeId v = 0; v < ph.nodes.size(); ++v) {
    const auto& node = ph.nodes[v];
    roots += node.parent == kNoNode ? 1 : 0;
    if (node.children.empty()) {
      // A leaf is a part: a state outside the atoms, or an atom.
      REQUIRE((node.state == sssp::paper::kNoState) == (ph.leaf_atom[v] != kNoAtom));
      if (node.parent != kNoNode && ph.leaf_atom[v] != kNoAtom) {
        REQUIRE(ph.atoms[ph.leaf_atom[v]].cutoff <= ph.nodes[node.parent].height);
      }
      continue;
    }
    REQUIRE(ph.leaf_atom[v] == kNoAtom);
    REQUIRE(node.children.size() >= 2);
    REQUIRE(node.links.size() + 1 == node.children.size());
    REQUIRE(W{0} < node.height);
    REQUIRE(node.height < part.node_upper[v]);
    if (node.parent != kNoNode) REQUIRE(node.height <= ph.nodes[node.parent].height);
    for (const auto& link : node.links) {
      // Connectivity by arcs strictly below the class upper end 2a (l.2098-2099): c0 = 2.
      REQUIRE(link.length < part.node_upper[v]);
      REQUIRE(link.arc < ph.graph.arcs.size());
      REQUIRE(ph.graph.arcs[link.arc].length == link.length);
    }
  }
  REQUIRE(part.roots.size() == roots);
  for (StateId x = 0; x < st.num_states(); ++x) {
    const NodeId leaf = part.leaf_of_state[x];
    REQUIRE(leaf < ph.nodes.size());
    REQUIRE(ph.leaf_atom[leaf] == sc.state_atom[x]);
    if (sc.state_atom[x] == kNoAtom) REQUIRE(ph.nodes[leaf].state == x);
  }

  // S9 (Lemma 5.15): one Alg 4 call per atom; the result satisfies Def 4.1 with c0 <= 2.
  sssp::Diagnostics sub_diag;
  const auto sub = sssp::paper::substitute_atoms(ph, &sub_diag, local);
  REQUIRE(sub.hierarchy.c0 <= 2);
  {
    const auto report = sssp::paper::validate_ordered_hierarchy(sub.hierarchy, r.c.g);
    INFO(report.summary());
    REQUIRE(report.ok());
  }
  const std::string base_case = sssp::paper::local_base_case_reason(local);
  for (const auto& fb : sub_diag.fallbacks) {
    REQUIRE(fb.reason == base_case);  // only Alg 4's recorded base case, never another search
    ++cov.atom_base_cases;
  }
  REQUIRE(sub_diag.fallbacks.size() <= sc.atoms.size());
  cov.atom_radial += sc.atoms.size() - sub_diag.fallbacks.size();

  // H4 as produced: the roots are listed in condensation order of the arcs between them.
  {
    std::vector<NodeId> root_of(ph.nodes.size(), kNoNode);
    std::vector<std::size_t> rank(ph.nodes.size(), 0);
    for (std::size_t i = 0; i < part.roots.size(); ++i) rank[part.roots[i]] = i;
    const auto root = [&](NodeId v) {
      NodeId top = v;
      while (root_of[top] == kNoNode && ph.nodes[top].parent != kNoNode) top = ph.nodes[top].parent;
      if (root_of[top] != kNoNode) top = root_of[top];
      for (NodeId w = v; w != top && root_of[w] == kNoNode; w = ph.nodes[w].parent)
        root_of[w] = top;
      return top;
    };
    for (const auto& arc : ph.graph.arcs) {
      REQUIRE(rank[root(part.leaf_of_state[arc.tail])] <= rank[root(part.leaf_of_state[arc.head])]);
    }
  }

  // Alg 8 l.14-15: one seed at the native inward copy of s (redirected into its atom), then the
  // minimum emitted value per vertex over all states projecting to it (Q-73).
  if (!distances) return;
  const std::size_t n = r.c.n;
  for (const Vertex s : {Vertex{0}, static_cast<Vertex>(n / 2), static_cast<Vertex>(n - 1)}) {
    const std::vector<sssp::paper::StateOffer<W>> seed{{st.native_inward[s], W{0}}};
    sssp::paper::SeededOrderedHierarchy<W> seeded{
        sub.hierarchy, {{s, W{0}}}, sssp::paper::redirect_offers<W>(sub, seed)};
    const auto ev = sssp::paper::evaluate(seeded);
    std::vector<W> dist(n, kInf<W>);
    for (StateId x = 0; x < sub.hierarchy.num_states(); ++x) {
      const Vertex v = st.projection[sub.old_state[x]];
      REQUIRE(sub.hierarchy.projection[x] == v);
      dist[v] = std::min(dist[v], ev.value[x]);
    }
    const auto oracle = sssp::baseline::bellman_ford(r.c.g, s);
    for (Vertex v = 0; v < n; ++v) {
      INFO("source " << s << " vertex " << v);
      REQUIRE(dist[v] == oracle.dist[v]);
    }
    ++cov.distance_checks;
  }
}

template <class W>
void check_all(to::OuterCase<W>& c, bool sccs, Coverage& cov,
               const sssp::paper::LocalHierarchyParams& local = sssp::paper::kPaperLocalParams,
               bool distances = true) {
  const auto r = to::run_scales(std::move(c));
  check_scales(r, cov);
  check_plan(r, build_outer_merge_plan(r.os, r.sc), sccs, cov);
  check_partition(r, local, distances, cov);
}

/// Exact P = 1 route names with 128-bit path sums (the all-pairs oracle of exact_names would
/// overflow on weights near 2^63): lambda_e = the least instance weight that is at least the
/// G_P distance of its endpoints.
std::pair<std::vector<std::uint64_t>, std::vector<std::uint32_t>> wide_names(
    const sssp::paper::pivots::PivotRoutingInstance<std::uint64_t>& inst) {
  using to::Exact;
  const std::size_t t = inst.t;
  std::vector<std::vector<Exact>> d(t, std::vector<Exact>(t, to::kExactInf));
  for (std::size_t i = 0; i < t; ++i) d[i][i] = 0;
  for (const auto& e : inst.edges) {
    d[e.u][e.v] = std::min<Exact>(d[e.u][e.v], e.weight);
    d[e.v][e.u] = std::min<Exact>(d[e.v][e.u], e.weight);
  }
  for (std::size_t k = 0; k < t; ++k) {
    for (std::size_t i = 0; i < t; ++i) {
      for (std::size_t j = 0; j < t; ++j) {
        if (d[i][k] != to::kExactInf && d[k][j] != to::kExactInf) {
          d[i][j] = std::min(d[i][j], d[i][k] + d[k][j]);
        }
      }
    }
  }
  std::vector<std::uint64_t> weights;
  for (const auto& e : inst.edges) weights.push_back(e.weight);
  std::ranges::sort(weights);
  std::vector<std::uint64_t> chosen;
  for (const auto& e : inst.edges) {
    std::uint64_t best = e.weight;
    for (const std::uint64_t a : weights) {
      if (a >= d[e.u][e.v]) {
        best = a;
        break;
      }
    }
    chosen.push_back(best);
  }
  std::vector<std::uint64_t> dictionary = chosen;
  std::ranges::sort(dictionary);
  dictionary.erase(std::unique(dictionary.begin(), dictionary.end()), dictionary.end());
  std::vector<std::uint32_t> name;
  for (const std::uint64_t c : chosen) {
    name.push_back(
        static_cast<std::uint32_t>(std::ranges::lower_bound(dictionary, c) - dictionary.begin()));
  }
  return {std::move(dictionary), std::move(name)};
}

/// Two light clusters (weights 1 and 2) joined by one tree edge of weight H > 2^63 and by
/// heavier nontree edges: with b = 3 and P = 1 (f = 31) a block of scale s_U <= 2 that reaches
/// a macro edge is unaligned (F^2 s_U <= 2^63 < a_U), and its edges of weight 1 form atoms.
std::vector<Edge<std::uint64_t>> two_scale(std::uint64_t seed, std::size_t& n) {
  SplitMix64 rng(seed);
  const std::size_t left = rng.between(3, 7);
  const std::size_t right = rng.between(3, 7);
  n = left + right;
  std::vector<Edge<std::uint64_t>> e;
  const auto v = [](std::size_t x) { return static_cast<Vertex>(x); };
  for (std::size_t x = 1; x < left; ++x) e.push_back({v(rng.below(x)), v(x), rng.between(1, 2)});
  for (std::size_t x = 1; x < right; ++x) {
    e.push_back({v(left + rng.below(x)), v(left + x), rng.between(1, 2)});
  }
  const std::uint64_t heavy = (std::uint64_t{1} << 63) + 1 + rng.below(4);
  e.push_back({v(rng.below(left)), v(left + rng.below(right)), heavy});
  const std::size_t extra = rng.between(1, 4);
  for (std::size_t i = 0; i < extra; ++i) {
    e.push_back({v(rng.below(left)), v(left + rng.below(right)), heavy + 1 + rng.below(6)});
  }
  for (std::size_t i = 0; i < 2; ++i) {
    const std::size_t a = rng.below(left);
    const std::size_t b = rng.below(left);
    if (a != b) e.push_back({v(a), v(b), rng.between(1, 2)});
  }
  return e;
}

void report(const Coverage& cov) {
  WARN("cases " << cov.cases << ", aligned " << cov.aligned << ", unaligned " << cov.unaligned
                << ", unattached " << cov.unattached << ", atoms " << cov.atoms << " (aligned "
                << cov.aligned_atoms << ", unaligned " << cov.unaligned_atoms << ", max states "
                << cov.max_atom_states << "), lowered supports " << cov.lowered_supports
                << ", separate events " << cov.separate_events << ", global events "
                << cov.global_events << ", classes " << cov.classes << ", SCC checks "
                << cov.scc_checks << ", merges " << cov.merges);
}

}  // namespace

TEMPLATE_TEST_CASE("Sec 8.2 / Lemma 8.4: scales, atoms and SCCs on seeded graphs (SC-11, SC-12)",
                   "[property][outer_scales]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  to::for_each_case<W>({1, 2, 3, 5, 8, 13, 20, 32}, 2,
                       [&](to::OuterCase<W>& c) { check_all(c, true, cov); });
  report(cov);
  CHECK(cov.cases > 0);
  CHECK(cov.aligned > 0);
  CHECK(cov.unattached > 0);
  CHECK(cov.atoms > 0);
  CHECK(cov.aligned_atoms > 0);
  CHECK(cov.max_atom_states >= 3);
  CHECK(cov.lowered_supports > 0);
  CHECK(cov.separate_events > 0);
  CHECK(cov.global_events > 0);
  CHECK(cov.merges > 0);
}

TEMPLATE_TEST_CASE("Sec 8.2 / Lemma 8.4: scales and atoms up to N = 200 (SC-11)",
                   "[property][outer_scales]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  // Alg 4 parameters {2, 2, 4} (Q-57; Mode B values, D-7): an atom with at least 99 states takes
  // Alg 4's radial branch instead of the base case, so Lemma 5.15 is exercised with both.
  const sssp::paper::LocalHierarchyParams local{2, 2, 4};
  to::for_each_case<W>({64, 200}, 1, [&](to::OuterCase<W>& c) { check_all(c, false, cov, local); });
  report(cov);
  CHECK(cov.aligned > 0);
  CHECK(cov.atoms > 0);
  CHECK(cov.max_atom_states >= 99);
  CHECK(cov.atom_base_cases > 0);
  CHECK(cov.atom_radial > 0);
  CHECK(cov.distance_checks > 0);
}

TEST_CASE("Sec 8.2 / Lemma 8.4: unaligned blocks and their atoms on uint64 (P = 1, b = 3)",
          "[property][outer_scales]") {
  using W = std::uint64_t;
  Coverage cov;
  for (std::uint64_t seed = 1; seed <= 300; ++seed) {
    std::size_t n = 0;
    const auto edges = two_scale(seed * 2654435761ULL, n);
    INFO("two-scale seed=" << seed << " n=" << n);
    const auto g = sssp::CsrGraph<W>::from_edges(n, edges);
    const auto od = sssp::paper::owner_decomposition(g, 3);
    const auto pe = sssp::paper::owner::projected_interowner_edges(g, od);
    auto [dictionary, name] = wide_names(pe.instance);
    auto c = to::make_outer_case_named<W>(n, edges, 3, std::move(dictionary), std::move(name));
    // Everything up to Def 4.1 on the substituted hierarchy. The distances are not compared
    // here: Evaluate's Sec 4 bookkeeping sums do not fit 64 bits for heights near 2^63 (it
    // throws std::overflow_error), which is outside Sec 8.2.
    check_all(c, true, cov, sssp::paper::kPaperLocalParams, false);
  }
  report(cov);
  CHECK(cov.distance_checks == 0);
  CHECK(cov.unaligned > 0);
  CHECK(cov.unaligned_atoms > 0);
  CHECK(cov.aligned > 0);
  CHECK(cov.merges > 0);
}
