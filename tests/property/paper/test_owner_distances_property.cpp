// Property tests for Sec 7.4 / Algorithm 7 (OwnerDistances) of arXiv 2611.34911v1 on seeded
// auxiliary graphs, against brute-force references written in this file in exact integer
// arithmetic (unsigned __int128; std::uint64_t lengths as is, quarter-integer doubles times 4).
//
// Checked statements:
//   * Alg 7 l.4-5: families are exactly the classes of (tree piece, width calendar) of the
//     representatives; `preorder` is a DFS preorder of T' (parents first, subtrees contiguous);
//     representatives increase in preorder inside a family; gaps equal the T' distances of
//     consecutive representatives (parent walking); horizons equal max (R_rep + L_Y); the family
//     horizon bound R_v <= R_rep + L_Y holds for every member (Sec 7.4 l.1786-1794); the number
//     of activation records equals the cells of the family calendars with start <= horizon.
//   * Lemma 7.8: every effective decrease (interatom exits and family proposals) arrives before
//     the starting marker of its destination cell; interatom arrivals satisfy z + w' <= R_v,
//     respect the direction rule between graph pieces, and join distinct atoms.
//   * Lemma 7.7 / 4.6: no emitted value is below the true J (an O(n^2) Dijkstra on G' written in
//     the test), J is the least emission, and every local nominal slot is scanned at most twice.
//   * Lemma 7.9 (deadline): z <= J_v <= z + L_Y for every member of a registered atom, where z
//     is its registration time; each atom is visited in at most L_Y / q + 2 cells.
//   * Lemma 7.9 (work), only under SSSP_DIAGNOSTICS (counters are compiled out otherwise):
//     atom cell visits = sum of the per-atom counts = Lemma 7.7's C (eval_cells) <= 4 N;
//     unsuccessful run checks <= activation records + sum_i (4 + (2 L_i + g(i-1, i) + g(i, i+1))
//     / q) (the charging argument of l.1806-1812); activation records <= 300 N (Lemma 7.4-style
//     bound: 2 per family + 136 per tree piece + 128 per retained tree edge + 2 |Y| per family).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/local/local_hierarchy.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/calendars.hpp"
#include "paper/owner/owner_distances.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::LocalHierarchyParams;
using sssp::paper::owner::ArrivalKind;
using sssp::paper::owner::AuxGraph;
using sssp::paper::pivots::PivotRoutingInstance;
using sssp::paper::pivots::RouteNames;
namespace lt = sssp::test::local;
namespace ow = sssp::paper::owner;

namespace {

using Exact = unsigned __int128;
constexpr Exact kScale = 128;

/// Exact integer image of a length: x (integers) or 4 x (quarter-integer doubles).
template <class W>
Exact exact(W x) {
  if constexpr (std::is_floating_point_v<W>) {
    const W y = x * W{4};
    REQUIRE(y == std::floor(y));
    REQUIRE(y < W{0x1p100});
    return static_cast<Exact>(y);
  } else {
    return static_cast<Exact>(x);
  }
}

std::vector<std::size_t> separator_sizes(std::size_t n) {
  std::size_t lg = 0;
  while ((std::size_t{1} << lg) < n) ++lg;
  std::set<std::size_t> s{1, 2, 3, 4, std::max<std::size_t>(lg, 1), n, n + 1};
  return {s.begin(), s.end()};
}

template <class W>
RouteNames<W> exact_names(const PivotRoutingInstance<W>& inst) {
  const auto dist = sssp::test::floyd_warshall<W>(inst.t, inst.edges);
  std::vector<W> weights;
  for (const auto& e : inst.edges) weights.push_back(e.weight);
  std::ranges::sort(weights);
  std::vector<W> chosen;
  for (const auto& e : inst.edges) {
    W best = e.weight;
    for (const W a : weights) {
      if (a > e.weight) break;
      if (a >= dist[e.u][e.v]) {
        best = a;
        break;
      }
    }
    chosen.push_back(best);
  }
  RouteNames<W> names;
  names.p_log2 = 0;
  names.dictionary = chosen;
  std::ranges::sort(names.dictionary);
  names.dictionary.erase(std::unique(names.dictionary.begin(), names.dictionary.end()),
                         names.dictionary.end());
  for (const W c : chosen) {
    const auto it = std::ranges::lower_bound(names.dictionary, c);
    names.name.push_back(static_cast<std::uint32_t>(it - names.dictionary.begin()));
  }
  return names;
}

template <class W>
AuxGraph<W> make_aux_graph(const std::vector<Edge<W>>& edges, std::size_t n, std::size_t b,
                           bool exact_p1, std::uint64_t seed) {
  const auto g = CsrGraph<W>::from_edges(n, edges);
  const auto od = sssp::paper::owner_decomposition(g, b);
  const auto pe = ow::projected_interowner_edges(g, od);
  ow::OuterEdgeClassification<W> cls;
  if (exact_p1) {
    cls = ow::classify_outer_edges(g, od, pe, exact_names(pe.instance));
  } else {
    sssp::paper::model::WordRng rng(seed);
    cls = ow::classify_interowner_edges(g, od, rng);
  }
  return ow::build_aux_graph(g, od, cls);
}

/// Light clusters under heavy spokes (as in test_calendars_property.cpp).
template <class W>
std::vector<Edge<W>> heavy_spokes(std::size_t n, std::uint64_t seed) {
  sssp::test::SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  if (n <= 1) return e;
  const std::size_t clusters = 1 + rng.below(3);
  std::vector<std::vector<Vertex>> members(clusters);
  for (Vertex v = 1; v < n; ++v) members[rng.below(clusters)].push_back(v);
  for (const auto& mem : members) {
    if (mem.empty()) continue;
    e.push_back({0, mem[0], lt::to_weight<W>(1000ULL << rng.between(0, 10), rng)});
    for (std::size_t i = 1; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(i)], mem[i], lt::to_weight<W>(rng.between(1, 3), rng)});
    }
    for (std::size_t i = 0; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(mem.size())], mem[rng.below(mem.size())],
                   lt::to_weight<W>(rng.between(1, 3), rng)});
    }
  }
  return e;
}

/// Chains of clusters of 37..44 vertices: cluster c > 0 hangs below cluster c - 1 by a medium
/// edge of weight 200..400 >= 2q (one tree piece and one width: one family of several atoms),
/// plus cross edges. Every atom is a union of clusters, so LocalHierarchyParams{1, 2, b} takes
/// Alg 4's radial branch.
template <class W>
std::vector<Edge<W>> cluster_chain(std::size_t clusters, std::uint64_t seed, std::size_t& n) {
  sssp::test::SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  std::vector<std::vector<Vertex>> members(clusters);
  Vertex next = 1;
  for (std::size_t c = 0; c < clusters; ++c) {
    const std::size_t size = rng.between(37, 44);
    for (std::size_t i = 0; i < size; ++i) members[c].push_back(next++);
  }
  n = next;
  for (std::size_t c = 0; c < clusters; ++c) {
    const auto& mem = members[c];
    if (c == 0) {
      e.push_back({0, mem[0], lt::to_weight<W>(2000, rng)});
    } else {
      const auto& up = members[c - 1];
      e.push_back({up[rng.below(up.size())], mem[0], lt::to_weight<W>(rng.between(200, 400), rng)});
    }
    for (std::size_t i = 1; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(i)], mem[i], lt::to_weight<W>(rng.between(1, 3), rng)});
    }
    for (std::size_t i = 0; i < mem.size() / 4; ++i) {
      e.push_back({mem[rng.below(mem.size())], mem[rng.below(mem.size())],
                   lt::to_weight<W>(rng.between(1, 60), rng)});
    }
  }
  for (std::size_t i = 0; i + 1 < clusters; ++i) {
    const auto& a = members[i];
    const auto& b = members[i + 1 + rng.below(clusters - i - 1)];
    e.push_back({a[rng.below(a.size())], b[rng.below(b.size())],
                 lt::to_weight<W>(rng.between(100, 600), rng)});
    const auto& m = members[rng.below(clusters)];
    e.push_back({0, m[rng.below(m.size())], lt::to_weight<W>(rng.between(2000, 4000), rng)});
  }
  return e;
}

/// O(n^2) Dijkstra on G' in exact images (independent of src/).
template <class W>
std::vector<Exact> exact_distances(const AuxGraph<W>& aux) {
  const std::size_t n = aux.num_vertices();
  constexpr Exact kInf = ~Exact{0};
  std::vector<Exact> d(n, kInf);
  std::vector<std::uint8_t> done(n, 0);
  d[ow::kAuxRoot] = 0;
  for (std::size_t it = 0; it < n; ++it) {
    std::size_t u = n;
    for (std::size_t v = 0; v < n; ++v) {
      if (done[v] == 0 && d[v] != kInf && (u == n || d[v] < d[u])) u = v;
    }
    if (u == n) break;
    done[u] = 1;
    for (const auto& arc : aux.graph.arcs(static_cast<Vertex>(u))) {
      const Exact c = d[u] + exact(arc.weight);
      if (c < d[arc.to]) d[arc.to] = c;
    }
  }
  return d;
}

/// T' distance by parent walking.
template <class W>
Exact tree_distance(const AuxGraph<W>& aux, Vertex x, Vertex y) {
  std::map<Vertex, Exact> up;
  Exact s = 0;
  for (Vertex v = x;; v = aux.tree_parent[v]) {
    up.emplace(v, s);
    if (aux.tree_parent[v] == kNoVertex) break;
    s += exact(aux.graph.edge(aux.tree_parent_edge[v]).weight);
  }
  Exact t = 0;
  for (Vertex v = y;; v = aux.tree_parent[v]) {
    const auto it = up.find(v);
    if (it != up.end()) return it->second + t;
    t += exact(aux.graph.edge(aux.tree_parent_edge[v]).weight);
  }
}

template <class W>
std::uint64_t counter([[maybe_unused]] const Diagnostics& d, [[maybe_unused]] const char* name) {
  const auto it = d.counters.find(name);
  return it == d.counters.end() ? 0 : it->second;
}

struct Coverage {
  std::uint64_t cases = 0;
  std::uint64_t arrivals = 0;
  std::uint64_t supplemental = 0;
  std::uint64_t multi_atom_families = 0;
  std::uint64_t radial_atoms = 0;
  std::uint64_t resumed_slots = 0;
  std::uint64_t max_slot_scans = 0;
  double max_activation_per_n = 0;
  double max_visits_per_n = 0;
};

template <class W>
void check_instance(const AuxGraph<W>& aux, const LocalHierarchyParams& lp, Coverage& cov) {
  Diagnostics diag;
  ow::OwnerDistancesOptions opt;
  opt.local = lp;
  opt.record_arrivals = true;
  opt.record_emissions = true;
  opt.count_slot_scans = true;
  const auto res = ow::owner_distances(aux, &diag, opt);
  const std::size_t na = aux.num_vertices();
  const auto& cal = res.calendars;
  const auto& atoms = res.atoms;
  const auto& fam = res.families;
  const std::size_t num_atoms = atoms.num_atoms();
  const auto J = exact_distances(aux);

  // J and the emissions (Lemma 4.6 through Lemma 7.7).
  REQUIRE(res.scaled_distance.size() == na);
  for (Vertex v = 0; v < na; ++v) REQUIRE(exact(res.scaled_distance[v]) == kScale * J[v]);
  for (Vertex v = 0; v < na; ++v) REQUIRE(exact(res.distance[v]) == J[v]);
  std::vector<Exact> least(na, ~Exact{0});
  for (const auto& [v, value] : res.emissions) {
    REQUIRE(v < na);
    REQUIRE(exact(value) >= kScale * J[v]);
    least[v] = std::min(least[v], exact(value));
  }
  for (Vertex v = 1; v < na; ++v) REQUIRE(least[v] == kScale * J[v]);
  REQUIRE(res.max_slot_scans <= 2);

  // Alg 7 l.4-5: preorder of T'.
  REQUIRE(fam.preorder.size() == na);
  std::vector<std::size_t> subtree(na, 1);
  for (std::size_t k = na; k-- > 1;) {
    const Vertex v = aux.tree_order[k];
    subtree[aux.tree_parent[v]] += subtree[v];
  }
  std::vector<std::uint8_t> used(na, 0);
  for (Vertex v = 0; v < na; ++v) {
    REQUIRE(fam.preorder[v] < na);
    REQUIRE(used[fam.preorder[v]] == 0);
    used[fam.preorder[v]] = 1;
    if (v == ow::kAuxRoot) continue;
    const Vertex p = aux.tree_parent[v];
    REQUIRE(fam.preorder[v] > fam.preorder[p]);
    REQUIRE(fam.preorder[v] + subtree[v] <= fam.preorder[p] + subtree[p]);
  }
  REQUIRE(fam.preorder[ow::kAuxRoot] == 0);

  // Families: the classes of (tree piece, calendar) of the representatives.
  REQUIRE(fam.family.size() == num_atoms);
  REQUIRE(fam.order.size() == num_atoms);
  REQUIRE(fam.offsets.size() == fam.num_families() + 1);
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> key_family;
  for (std::uint32_t a = 0; a < num_atoms; ++a) {
    const Vertex rep = atoms.representative[a];
    const auto key = std::make_pair(cal.tree_piece[rep], cal.calendar[rep]);
    const auto [it, fresh] = key_family.try_emplace(key, fam.family[a]);
    REQUIRE(it->second == fam.family[a]);
    REQUIRE(fam.tree_piece[fam.family[a]] == key.first);
    REQUIRE(fam.calendar[fam.family[a]] == key.second);
  }
  REQUIRE(key_family.size() == fam.num_families());
  std::uint64_t expected_records = 0;
  for (std::uint32_t f = 0; f < fam.num_families(); ++f) {
    REQUIRE(fam.offsets[f] < fam.offsets[f + 1]);
    Exact horizon = 0;
    for (std::uint32_t k = fam.offsets[f]; k < fam.offsets[f + 1]; ++k) {
      const std::uint32_t a = fam.order[k];
      REQUIRE(fam.family[a] == f);
      REQUIRE(fam.position[a] == k);
      const Vertex rep = atoms.representative[a];
      const Exact reach = exact(aux.owner_distance[rep]) + exact(atoms.length[a]);
      horizon = std::max(horizon, reach);
      for (const Vertex v : atoms.members(a)) REQUIRE(exact(aux.owner_distance[v]) <= reach);
      if (k + 1 < fam.offsets[f + 1]) {
        const Vertex next = atoms.representative[fam.order[k + 1]];
        REQUIRE(fam.preorder[rep] < fam.preorder[next]);
        REQUIRE(exact(fam.gap[k]) == kScale * tree_distance(aux, rep, next));
      } else {
        REQUIRE(fam.gap[k] == W{0});
      }
    }
    REQUIRE(exact(fam.horizon[f]) == kScale * horizon);
    if (fam.offsets[f + 1] - fam.offsets[f] > 1) ++cov.multi_atom_families;
    const auto& wc = cal.calendars[fam.calendar[f]];
    const Exact q = exact(wc.width);
    std::uint64_t cells = 0;
    for (std::size_t i = 0; i < wc.num_cells(); ++i) {
      if (q * i > kScale * horizon) break;
      ++cells;
    }
    expected_records += cells;
  }
  REQUIRE(res.activation_records == expected_records);

  // Lemma 7.8: the arrival log.
  for (const auto& r : res.arrivals) {
    INFO("arrival tail " << r.tail << " head " << r.head << " marker " << r.marker);
    REQUIRE(r.on_time);
    REQUIRE(r.head < na);
    REQUIRE(exact(r.value) <= kScale * exact(aux.owner_distance[r.head]));
    REQUIRE(exact(r.value) >= kScale * J[r.head]);
    if (r.kind == ArrivalKind::kSupplemental) {
      ++cov.supplemental;
      continue;
    }
    REQUIRE(r.tail < na);
    REQUIRE(atoms.atom[r.tail] != atoms.atom[r.head]);
    if (cal.graph_piece[r.tail] != cal.graph_piece[r.head]) {
      // Direction rule: only from the smaller range, r_head > 16 R_tail.
      REQUIRE(exact(aux.radius[r.head]) > 16 * exact(aux.owner_distance[r.tail]));
    }
  }
  cov.arrivals += res.arrivals.size();

  // Lemma 7.9: registrations and deadlines.
  REQUIRE(res.registration.size() == num_atoms);
  REQUIRE(res.atom_cells.size() == num_atoms);
  std::uint64_t visits = 0;
  [[maybe_unused]] Exact run_bound = res.activation_records;
  for (std::uint32_t a = 0; a < num_atoms; ++a) {
    REQUIRE(res.registration[a] != sssp::WeightTraits<W>::infinity());
    const Exact z = exact(res.registration[a]);
    const Exact len = kScale * exact(atoms.length[a]);
    for (const Vertex v : atoms.members(a)) {
      REQUIRE(z <= kScale * J[v]);
      REQUIRE(kScale * J[v] <= z + len);
    }
    const Exact q = exact(cal.width[atoms.representative[a]]);
    REQUIRE(static_cast<Exact>(res.atom_cells[a]) * q <= len + 2 * q);
    visits += res.atom_cells[a];
    const std::uint32_t k = fam.position[a];
    const std::uint32_t f = fam.family[a];
    Exact around = 2 * len;
    if (k > fam.offsets[f]) around += exact(fam.gap[k - 1]);
    if (k + 1 < fam.offsets[f + 1]) around += exact(fam.gap[k]);
    run_bound += 4 + around / q;
    if (res.radial_atom[a] != 0) ++cov.radial_atoms;
  }

  // Fallbacks: one Alg 4 base-case record per atom that did not take the radial branch.
  std::size_t base = 0;
  for (std::uint32_t a = 0; a < num_atoms; ++a) {
    if (res.radial_atom[a] == 0) ++base;
  }
  REQUIRE(diag.fallbacks.size() == base);

  const double n = static_cast<double>(std::max<std::size_t>(na, 1));
  cov.max_activation_per_n =
      std::max(cov.max_activation_per_n, static_cast<double>(res.activation_records) / n);
  cov.max_visits_per_n = std::max(cov.max_visits_per_n, static_cast<double>(visits) / n);
  cov.max_slot_scans = std::max<std::uint64_t>(cov.max_slot_scans, res.max_slot_scans);
  ++cov.cases;
#if SSSP_DIAGNOSTICS
  REQUIRE(counter<W>(diag, "owner_distances.atom_cell_visits") == visits);
  REQUIRE(counter<W>(diag, "eval_cells") == visits);
  REQUIRE(visits <= 4 * na);
  REQUIRE(counter<W>(diag, "owner_distances.registrations") == num_atoms);
  REQUIRE(counter<W>(diag, "owner_distances.activation_records") == res.activation_records);
  REQUIRE(res.activation_records <= 300 * na);
  REQUIRE(static_cast<Exact>(counter<W>(diag, "owner_distances.run_checks_failed")) <= run_bound);
  REQUIRE(counter<W>(diag, "owner_distances.run_checks") ==
          counter<W>(diag, "owner_distances.run_checks_failed") + num_atoms);
  cov.resumed_slots += counter<W>(diag, "eval_resumed_slot_visits");
#endif
}

}  // namespace

TEMPLATE_TEST_CASE("Alg 7 invariants (Lemmas 7.7-7.9) on seeded auxiliary graphs",
                   "[property][owner_distances]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  const lt::Family families[] = {lt::Family::kRandomTree, lt::Family::kSparse, lt::Family::kGrid,
                                 lt::Family::kComplete};
  const lt::Weights weights[] = {lt::Weights::kTies, lt::Weights::kUniform,
                                 lt::Weights::kExponential};
  const std::size_t sizes[] = {2, 7, 20, 33};
  for (const auto f : families) {
    for (const auto wd : weights) {
      for (const std::size_t n : sizes) {
        const std::uint64_t gseed =
            n * 131 + static_cast<std::uint64_t>(f) * 17 + static_cast<std::uint64_t>(wd) + 5;
        const auto edges = lt::make_graph<W>(f, wd, n, gseed);
        for (const std::size_t b : separator_sizes(n)) {
          for (const bool exact_p1 : {true, false}) {
            INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " seed=" << gseed << " b=" << b
                             << (exact_p1 ? " exact-P1" : " name-routes-paper"));
            check_instance(make_aux_graph(edges, n, b, exact_p1, gseed + b),
                           sssp::paper::kPaperLocalParams, cov);
          }
        }
      }
    }
  }
  for (const std::size_t n : {7, 20, 33, 60}) {
    for (std::uint64_t seed = 1; seed <= 3; ++seed) {
      const std::uint64_t gseed = seed * 7919 + n;
      const auto edges = heavy_spokes<W>(n, gseed);
      for (const std::size_t b : separator_sizes(n)) {
        INFO("heavy-spokes n=" << n << " seed=" << gseed << " b=" << b);
        check_instance(make_aux_graph(edges, n, b, true, gseed), sssp::paper::kPaperLocalParams,
                       cov);
      }
    }
  }
  INFO("cases " << cov.cases << ", arrivals " << cov.arrivals << ", family decreases "
                << cov.supplemental << ", multi-atom families " << cov.multi_atom_families);
  CHECK(cov.arrivals > 0);
  CHECK(cov.supplemental > 0);
  CHECK(cov.multi_atom_families > 0);
  WARN("owner distances (Mode A atoms): cases "
       << cov.cases << ", max activation records/N " << cov.max_activation_per_n
       << ", max atom cell visits/N " << cov.max_visits_per_n << ", max scans per slot "
       << cov.max_slot_scans);
}

TEMPLATE_TEST_CASE("Alg 7 invariants with radial atoms (Mode B atom parameters, Q-04)",
                   "[property][owner_distances]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  for (std::size_t clusters = 2; clusters <= 4; ++clusters) {
    for (std::uint64_t seed = 1; seed <= 2; ++seed) {
      std::size_t n = 0;
      const std::uint64_t gseed = seed * 65537 + clusters;
      const auto edges = cluster_chain<W>(clusters, gseed, n);
      for (const std::size_t local_b : {std::size_t{2}, std::size_t{5}}) {
        INFO("cluster-chain k=" << clusters << " seed=" << gseed << " local b=" << local_b);
        check_instance(make_aux_graph(edges, n, n + 1, true, gseed),
                       LocalHierarchyParams{1, 2, local_b}, cov);
      }
    }
  }
  INFO("cases " << cov.cases << ", radial atoms " << cov.radial_atoms << ", multi-atom families "
                << cov.multi_atom_families << ", resumed slot visits " << cov.resumed_slots);
  CHECK(cov.radial_atoms > 0);
  CHECK(cov.multi_atom_families > 0);
  WARN("owner distances (radial atoms): cases "
       << cov.cases << ", max activation records/N " << cov.max_activation_per_n
       << ", max atom cell visits/N " << cov.max_visits_per_n << ", max scans per slot "
       << cov.max_slot_scans << ", resumed slot visits " << cov.resumed_slots);
}
