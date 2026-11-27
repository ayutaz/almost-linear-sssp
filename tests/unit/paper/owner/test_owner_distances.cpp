// Sec 7.4 and Algorithm 7 (OwnerDistances) of arXiv 2611.34911v1 on hand-computed auxiliary
// graphs: families, representative preorder, gaps and horizons (Alg 7 l.4-5), the run splitting
// with padded family proposals (l.1801-1815), root offers (l.9, Q-47), interatom exits with the
// direction rule and the endpoint caps (l.1819-1824, Lemma 7.8), and J (l.14).
//
// Calendar values are 2^7 times the lengths they stand for (calendars.hpp); the expectations of
// scaled quantities are written as 128 * (real value). Every expected J below is the auxiliary
// distance d_{w'}(v, a) computed by hand.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/calendars.hpp"
#include "paper/owner/owner_distances.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::owner::ArrivalKind;
using sssp::paper::owner::AuxGraph;
using sssp::paper::owner::kAuxRoot;
using sssp::paper::owner::owner_distances;
using sssp::paper::owner::OwnerDistancesOptions;

namespace {

/// G' with root 0 and the given T' edge ids; T' is rooted by a BFS over the tree arcs in CSR
/// order, r / R are the root bottleneck / distance along it (Sec 7.3 l.1693). (Same builder as
/// tests/unit/paper/owner/test_calendars.cpp.)
template <class W>
AuxGraph<W> make_aux(std::size_t n, std::vector<Edge<W>> edges, std::vector<EdgeId> tree,
                     std::size_t hop_bound) {
  AuxGraph<W> aux;
  aux.hop_bound = hop_bound;
  aux.graph = CsrGraph<W>::from_edges(n, std::move(edges));
  for (EdgeId i = 0; i < aux.graph.num_edges(); ++i) aux.input_edge.push_back(i);
  aux.aux_edge = aux.input_edge;
  aux.input_vertex.assign(n, kNoVertex);
  aux.aux_vertex.assign(n, kNoVertex);
  for (Vertex v = 1; v < n; ++v) aux.input_vertex[v] = aux.aux_vertex[v] = v;
  aux.aux_vertex[0] = kAuxRoot;
  aux.tree_edges = tree;
  std::vector<std::uint8_t> in_tree(aux.graph.num_edges(), 0);
  for (const EdgeId t : tree) in_tree[t] = 1;
  aux.tree_parent.assign(n, kNoVertex);
  aux.tree_parent_edge.assign(n, kNoEdge);
  aux.radius.assign(n, W{0});
  aux.owner_distance.assign(n, W{0});
  std::vector<std::uint8_t> seen(n, 0);
  aux.tree_order.push_back(kAuxRoot);
  seen[kAuxRoot] = 1;
  for (std::size_t k = 0; k < aux.tree_order.size(); ++k) {
    const Vertex x = aux.tree_order[k];
    for (const auto& arc : aux.graph.arcs(x)) {
      if (in_tree[arc.edge] == 0 || seen[arc.to] != 0) continue;
      seen[arc.to] = 1;
      aux.tree_parent[arc.to] = x;
      aux.tree_parent_edge[arc.to] = arc.edge;
      aux.radius[arc.to] = x == kAuxRoot ? arc.weight : std::max(aux.radius[x], arc.weight);
      aux.owner_distance[arc.to] = aux.owner_distance[x] + arc.weight;
      aux.tree_order.push_back(arc.to);
    }
  }
  REQUIRE(aux.tree_order.size() == n);
  return aux;
}

template <class W>
W scaled(W x) {
  return W{128} * x;
}

template <class W>
std::uint64_t counter(const Diagnostics& d, const std::string& name) {
  const auto it = d.counters.find(name);
  return it == d.counters.end() ? 0 : it->second;
}

/// Every fallback is Alg 4's base case of one atom (Mode A parameters).
void require_base_case_records(const Diagnostics& diag, std::size_t atoms, std::size_t vertices) {
  REQUIRE(diag.fallbacks.size() == atoms);
  std::size_t covered = 0;
  for (const auto& f : diag.fallbacks) {
    CHECK(f.reason == "local: k <= 2*L_k^4 (paper base case)");
    covered += f.vertices;
  }
  CHECK(covered == vertices);
}

// Example A (singleton atoms, one family, a root shortcut and a padded proposal).
//   id: 0 (0,1) 10   1 (1,2) 10   2 (2,3) 10   3 (3,4) 10   4 (0,3) 11;  T' = {0,1,2,3}.
//   r = 10, R = 10, 20, 30, 40. Grid origin 10, q = 128 * 10 / 32 = 40 (real 0.3125): no
//   defining edge (w' < 0.625), four singleton atoms in one family (P0, q), order 1, 2, 3, 4,
//   gaps 10, 10, 10, horizon max R = 40; cells k * 0.3125, starts <= 40 for k = 0..128.
//   Run: cell 32 registers {1} (z 10), proposal to rep 2: u = 20 = R_2, no decrease. Cell 35
//   [10.9375, 11.25) registers {3} (z 11 from the root edge 0-3); proposals: to 2 u = 21 > R_2
//   (capped), to 4 u = 21 < 40 (decrease, destination cell 67). J = 0, 10, 20, 11, 21.
template <class W>
AuxGraph<W> example_a() {
  return make_aux<W>(5, {{0, 1, W{10}}, {1, 2, W{10}}, {2, 3, W{10}}, {3, 4, W{10}}, {0, 3, W{11}}},
                     {0, 1, 2, 3}, 5);
}

// Example B (one family of two atoms, an internal chord, an entry by a root edge, an interatom
// exit inside one graph piece).
//   id: 0 (0,1) 1000  1 (1,2) 4  2 (2,3) 2  3 (3,4) 100  4 (4,5) 4  5 (1,3) 7  6 (0,3) 1001;
//   T' = {0,1,2,3,4}. r = 1000, R = 1000, 1004, 1006, 1106, 1110. q = 4000 scaled (real 31.25),
//   defining edges w' < 62.5: 1, 2, 4 and the chord 5. Atoms {1,2,3} (L 6) and {4,5} (L 4), one
//   family, gap D_T'(1, 4) = 106, horizon max(1006, 1110) = 1110, 36 activation records.
//   Cell 32 [1000, 1031.25) registers {1,2,3} (z 1000): proposal to 4, u = 1112 > R_4 (capped);
//   it emits 1000, 1001 (root edge into 3), 1003 (3 -> 2); exit 3 -> 4 gives 1101 < R_4.
//   Cell 35 registers {4,5} (z 1101). J = 0, 1000, 1003, 1001, 1101, 1105.
template <class W>
AuxGraph<W> example_b() {
  return make_aux<W>(6,
                     {{0, 1, W{1000}},
                      {1, 2, W{4}},
                      {2, 3, W{2}},
                      {3, 4, W{100}},
                      {4, 5, W{4}},
                      {1, 3, W{7}},
                      {0, 3, W{1001}}},
                     {0, 1, 2, 3, 4}, 6);
}

// Example C (two graph pieces).
//   id: 0 (0,1) 1  1 (1,2) 100  2 (2,3) 1  3 (0,3) 101;  T' = {0,1,2}. r = 1, 100, 100,
//   R = 1, 101, 102. Edge 1 does not overlap (r_2 = 100 > 16 R_1): pieces {1} and {2,3}, the
//   interpiece arc has tail 1 and head 2, so the exit 2 -> 1 is dropped. Atoms {1} and {2,3}.
//   J = 0, 1, 101, 101.
template <class W>
AuxGraph<W> example_c() {
  return make_aux<W>(4, {{0, 1, W{1}}, {1, 2, W{100}}, {2, 3, W{1}}, {0, 3, W{101}}}, {0, 1, 2}, 4);
}

template <class W>
void require_arrivals_on_time(const sssp::paper::owner::OwnerDistancesResult<W>& res) {
  for (const auto& a : res.arrivals) {
    INFO("arrival at " << a.head << " value " << a.value << " marker " << a.marker);
    CHECK(a.on_time);
  }
}

}  // namespace

TEMPLATE_TEST_CASE("Alg 7 on example A: root shortcut and padded family proposals",
                   "[unit][owner_distances]", std::uint64_t, double) {
  using W = TestType;
  const auto aux = example_a<W>();
  Diagnostics diag;
  OwnerDistancesOptions opt;
  opt.record_arrivals = true;
  opt.count_slot_scans = true;
  const auto res = owner_distances(aux, &diag, opt);
  CHECK(res.distance == std::vector<W>{W{0}, W{10}, W{20}, W{11}, W{21}});
  CHECK(res.scaled_distance ==
        std::vector<W>{W{0}, scaled(W{10}), scaled(W{20}), scaled(W{11}), scaled(W{21})});

  // Alg 7 l.4-5.
  const auto& fam = res.families;
  REQUIRE(res.atoms.num_atoms() == 4);
  REQUIRE(fam.num_families() == 1);
  CHECK(fam.offsets == std::vector<std::uint32_t>{0, 4});
  CHECK(fam.order == std::vector<std::uint32_t>{0, 1, 2, 3});
  CHECK(fam.gap == std::vector<W>{scaled(W{10}), scaled(W{10}), scaled(W{10}), W{0}});
  CHECK(fam.horizon == std::vector<W>{scaled(W{40})});
  CHECK(fam.preorder == std::vector<std::uint32_t>{0, 1, 2, 3, 4});
  CHECK(res.activation_records == 129);

  // Registrations at the due minima; every singleton atom lives for one cell.
  CHECK(res.registration ==
        std::vector<W>{scaled(W{10}), scaled(W{20}), scaled(W{11}), scaled(W{21})});
  CHECK(res.atom_cells == std::vector<std::uint32_t>{1, 1, 1, 1});
  CHECK(res.radial_atom == std::vector<std::uint8_t>{0, 0, 0, 0});

  // The padded proposal to 4 is the only effective family decrease (Lemma 7.8).
  std::size_t supplemental = 0;
  for (const auto& a : res.arrivals) {
    if (a.kind != ArrivalKind::kSupplemental) continue;
    ++supplemental;
    CHECK(a.head == 4);
    CHECK(a.value == scaled(W{21}));
    CHECK(a.tail == kNoVertex);
  }
  CHECK(supplemental == 1);
  require_arrivals_on_time(res);
  CHECK(res.max_slot_scans <= 2);
  require_base_case_records(diag, 4, 4);
#if SSSP_DIAGNOSTICS
  CHECK(counter<W>(diag, "owner_distances.registrations") == 4);
  CHECK(counter<W>(diag, "owner_distances.supplemental_proposals") == 3);
  CHECK(counter<W>(diag, "owner_distances.supplemental_capped") == 1);
  CHECK(counter<W>(diag, "owner_distances.supplemental_decreases") == 1);
  CHECK(counter<W>(diag, "owner_distances.root_offers") == 2);
  CHECK(counter<W>(diag, "owner_distances.initial_offers") == 4);
  CHECK(counter<W>(diag, "owner_distances.activation_records") == 129);
  CHECK(counter<W>(diag, "owner_distances.atom_cell_visits") == 4);
#endif
}

TEMPLATE_TEST_CASE("Alg 7 on example B: a two-atom family with an internal chord",
                   "[unit][owner_distances]", std::uint64_t, double) {
  using W = TestType;
  const auto aux = example_b<W>();
  Diagnostics diag;
  OwnerDistancesOptions opt;
  opt.record_arrivals = true;
  opt.record_emissions = true;
  opt.count_slot_scans = true;
  const auto res = owner_distances(aux, &diag, opt);
  CHECK(res.distance == std::vector<W>{W{0}, W{1000}, W{1003}, W{1001}, W{1101}, W{1105}});

  REQUIRE(res.atoms.num_atoms() == 2);
  CHECK(res.atoms.length == std::vector<W>{W{6}, W{4}});
  const auto& fam = res.families;
  REQUIRE(fam.num_families() == 1);
  CHECK(fam.order == std::vector<std::uint32_t>{0, 1});
  CHECK(fam.gap == std::vector<W>{scaled(W{106}), W{0}});
  CHECK(fam.horizon == std::vector<W>{scaled(W{1110})});
  CHECK(res.activation_records == 36);
  CHECK(res.registration == std::vector<W>{scaled(W{1000}), scaled(W{1101})});
  CHECK(res.atom_cells == std::vector<std::uint32_t>{1, 1});

  // The exit 3 -> 4 is the only effective interatom decrease; the family proposal is capped.
  std::size_t interatom = 0;
  for (const auto& a : res.arrivals) {
    CHECK(a.kind == ArrivalKind::kInteratom);
    if (a.kind != ArrivalKind::kInteratom) continue;
    ++interatom;
    CHECK(a.tail == 3);
    CHECK(a.head == 4);
    CHECK(a.value == scaled(W{1101}));
  }
  CHECK(interatom == 1);
  require_arrivals_on_time(res);

  // No emission below J (Lemma 4.6 via Lemma 7.7), and every vertex emits.
  std::vector<std::uint8_t> emitted(aux.num_vertices(), 0);
  for (const auto& [v, value] : res.emissions) {
    CHECK(value >= res.scaled_distance[v]);
    emitted[v] = 1;
  }
  CHECK(emitted == std::vector<std::uint8_t>{0, 1, 1, 1, 1, 1});
  CHECK(res.max_slot_scans <= 2);
  // Atom sizes 3 and 2, internal edges 3 (chord included) and 1.
  REQUIRE(diag.fallbacks.size() == 2);
  CHECK(diag.fallbacks[0].vertices == 3);
  CHECK(diag.fallbacks[0].edges == 3);
  CHECK(diag.fallbacks[1].vertices == 2);
  CHECK(diag.fallbacks[1].edges == 1);
#if SSSP_DIAGNOSTICS
  CHECK(counter<W>(diag, "owner_distances.supplemental_proposals") == 1);
  CHECK(counter<W>(diag, "owner_distances.supplemental_capped") == 1);
  CHECK(counter<W>(diag, "owner_distances.internal_edges") == 4);
  CHECK(counter<W>(diag, "owner_distances.interatom_decreases") == 1);
#endif
}

TEMPLATE_TEST_CASE("Alg 7 on example C: two graph pieces and the direction rule",
                   "[unit][owner_distances]", std::uint64_t, double) {
  using W = TestType;
  const auto aux = example_c<W>();
  Diagnostics diag;
  OwnerDistancesOptions opt;
  opt.record_arrivals = true;
  const auto res = owner_distances(aux, &diag, opt);
  CHECK(res.distance == std::vector<W>{W{0}, W{1}, W{101}, W{101}});
  REQUIRE(res.calendars.interpiece_arcs.size() == 1);
  CHECK(res.atoms.num_atoms() == 2);
  CHECK(res.families.num_families() == 2);
  require_arrivals_on_time(res);
  require_base_case_records(diag, 2, 3);
#if SSSP_DIAGNOSTICS
  CHECK(counter<W>(diag, "owner_distances.exits_wrong_direction") == 1);
  CHECK(counter<W>(diag, "owner_distances.exits") == 1);
  CHECK(counter<W>(diag, "owner_distances.root_offers") == 2);
#endif
}

TEMPLATE_TEST_CASE("Alg 7 on degenerate auxiliary graphs", "[unit][owner_distances]", std::uint64_t,
                   double) {
  using W = TestType;
  SECTION("only the root") {
    Diagnostics diag;
    const auto res = owner_distances(make_aux<W>(1, {}, {}, 2), &diag);
    CHECK(res.distance == std::vector<W>{W{0}});
    CHECK(res.atoms.num_atoms() == 0);
    CHECK(res.families.num_families() == 0);
    CHECK(diag.fallbacks.empty());
  }
  SECTION("one child of the root") {
    Diagnostics diag;
    const auto res = owner_distances(make_aux<W>(2, {{1, 0, W{3}}}, {0}, 3), &diag);
    CHECK(res.distance == std::vector<W>{W{0}, W{3}});
    CHECK(res.atoms.num_atoms() == 1);
    CHECK(res.families.num_families() == 1);
    CHECK(res.atom_cells == std::vector<std::uint32_t>{1});
    require_base_case_records(diag, 1, 1);
  }
  SECTION("one atom, one family, one cell") {
    // A star of light edges under one spoke: r = 640 everywhere, q = 20 (real), the edges of
    // weight 1 and 2 are defining, all J lie in [640, 642] inside the cell [640, 660).
    const auto aux =
        make_aux<W>(4, {{0, 1, W{640}}, {1, 2, W{1}}, {1, 3, W{2}}, {2, 3, W{2}}}, {0, 1, 2}, 4);
    Diagnostics diag;
    OwnerDistancesOptions opt;
    opt.count_slot_scans = true;
    const auto res = owner_distances(aux, &diag, opt);
    CHECK(res.distance == std::vector<W>{W{0}, W{640}, W{641}, W{642}});
    CHECK(res.atoms.num_atoms() == 1);
    CHECK(res.families.num_families() == 1);
    CHECK(res.atom_cells == std::vector<std::uint32_t>{1});
    CHECK(res.max_slot_scans <= 2);
    require_base_case_records(diag, 1, 3);
  }
}

TEMPLATE_TEST_CASE("Alg 7 after the paper's outer separator (Mode A pipeline)",
                   "[unit][owner_distances]", std::uint64_t, double) {
  using W = TestType;
  // b = L_N^64 gives one pivot (vertex 0); G' is the input graph rooted at 0 (Q-48), J is the
  // distance from 0, and every atom takes Alg 4's base case (recorded).
  const std::vector<Edge<W>> edges{
      {0, 1, W{4}}, {1, 2, W{1}}, {2, 3, W{5000}}, {3, 4, W{2}}, {0, 4, W{9000}}};
  const auto g = CsrGraph<W>::from_edges(5, edges);
  const auto od = sssp::paper::owner_decomposition(
      g, sssp::paper::owner::outer_separator_size(g.num_vertices()));
  REQUIRE(od.pivots.size() == 1);
  sssp::paper::model::WordRng rng(7);
  const auto cls = sssp::paper::owner::classify_interowner_edges(g, od, rng);
  const auto aux = sssp::paper::owner::build_aux_graph(g, od, cls);
  Diagnostics diag;
  const auto res = owner_distances(aux, &diag);
  CHECK(res.distance == std::vector<W>{W{0}, W{4}, W{5}, W{5005}, W{5007}});
  require_base_case_records(diag, res.atoms.num_atoms(), 4);
}

TEMPLATE_TEST_CASE("Alg 7 rejects lengths whose 2^7 multiple does not fit",
                   "[unit][owner_distances]", std::uint64_t, double) {
  using W = TestType;
  const W huge =
      std::is_floating_point_v<W> ? static_cast<W>(1e307) : static_cast<W>(std::uint64_t{1} << 60);
  const auto aux = make_aux<W>(3, {{0, 1, W{10}}, {1, 2, W{10}}, {0, 2, huge}}, {0, 1}, 3);
  CHECK_THROWS_AS(owner_distances(aux), std::overflow_error);
}
