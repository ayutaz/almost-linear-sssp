// Sec 7.1-7.2 (Lemma 7.1, Lemma 7.2) of arXiv 2611.34911v1 on hand-computed examples: the
// projected interowner instance, the radial / macro classification from given route names, the
// auxiliary graph G' with the identified pivot root a, the owner-spoke tree T', and the owner
// coordinate J (reference Dijkstra on G').
//
// The separator sizes b used here are small so that the owner decomposition has several pivots on
// a small graph (Mode B style, Q-03). The paper's outer b = L_N^64 (outer_separator_size) gives a
// single pivot on every graph small enough to write by hand; that case is tested separately.
// Hand-written route names use P = 1 (p_log2 = 0) and are valid certificates in the sense of
// Definition 6.1 (checked with check_route_names against Floyd-Warshall on G_P).
//
// Counters are compiled out with SSSP_DIAGNOSTICS=0 (ci-release, release): assertions on them run
// only when they are counted. Every counted quantity is also checked on the returned structures,
// so the behaviour is tested in every build.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::owner_decomposition;
using sssp::paper::model::WordRng;
using sssp::paper::owner::build_aux_graph;
using sssp::paper::owner::check_aux_metric;
using sssp::paper::owner::classify_interowner_edges;
using sssp::paper::owner::classify_outer_edges;
using sssp::paper::owner::kAuxRoot;
using sssp::paper::owner::kNoInstanceEdge;
using sssp::paper::owner::kNoName;
using sssp::paper::owner::outer_separator_size;
using sssp::paper::owner::OuterEdgeClass;
using sssp::paper::owner::owner_coordinates_ref;
using sssp::paper::owner::projected_interowner_edges;
using sssp::paper::pivots::check_route_names;
using sssp::paper::pivots::RouteNames;

namespace {

// The 7-vertex example. Tree path 0-1-2-3-4-5-6 (edges 0..5) plus five extra edges.
//   id: 0 (0,1) 3   1 (1,2) 1   2 (2,3) 4   3 (3,4) 2   4 (4,5) 5   5 (5,6) 1
//       6 (0,4) 10  7 (2,6) 12  8 (0,6) 30  9 (4,4) 7 (loop)  10 (1,3) 5
// MST = edges 0..5 (edge 10 ties with edge 4 at weight 5 but closes the cycle 1-2-3).
// Reduced graph (w <= D_T): edge 6 (10 <= 10), 7 (12 <= 12), 10 (5 <= 5); edge 8 (30 > 16) and
// the loop are omitted.
// b = 2, postorder from root 0: 6 (1), 5 (2: mark), 4 (1), 3 (2: mark), 2 (1), 1 (2: mark), 0 (1).
// A = {1, 3, 5} (pivot indices 0, 1, 2), no branch vertices. Corridor 1-2-3 cuts edge 2 (4 > 1),
// corridor 3-4-5 cuts edge 4 (5 > 2). Cells {0,1,2} -> 1, {3,4} -> 3, {5,6} -> 5; h = 3.
// r = R = (3, 0, 1, 0, 2, 0, 1).
template <class W>
std::vector<Edge<W>> seven_vertex_edges() {
  return {{0, 1, W{3}},  {1, 2, W{1}},  {2, 3, W{4}},  {3, 4, W{2}}, {4, 5, W{5}}, {5, 6, W{1}},
          {0, 4, W{10}}, {2, 6, W{12}}, {0, 6, W{30}}, {4, 4, W{7}}, {1, 3, W{5}}};
}

// Projected interowner edges (input id order) and their lengths A_e = R_x + w + R_y:
//   edge 2: pivots (0,1), A = 1 + 4 + 0 = 5      edge 4: (1,2), A = 2 + 5 + 0 = 7
//   edge 6: (0,1), A = 3 + 10 + 2 = 15           edge 7: (0,2), A = 1 + 12 + 1 = 14
//   edge 10: (0,1), A = 0 + 5 + 0 = 5
// G_P distances: d(0,1) = 5, d(1,2) = 7, d(0,2) = 12. With P = 1 the least valid names are
// 5, 7, 5, 14, 5 (edge 7 needs lambda >= 12 among the instance weights <= 14).
template <class W>
RouteNames<W> seven_vertex_names() {
  RouteNames<W> names;
  names.dictionary = {W{5}, W{7}, W{14}};
  names.name = {0, 1, 0, 2, 0};
  names.p_log2 = 0;
  return names;
}

}  // namespace

TEST_CASE("outer separator size is L_N^64, saturating") {
  CHECK(outer_separator_size(0) == 1);  // L_0 = 1
  // L_1 = 2: 2^64 does not fit, the size saturates (every feasible N then has one pivot).
  CHECK(outer_separator_size(1) == std::numeric_limits<std::size_t>::max());
  CHECK(outer_separator_size(1000) == std::numeric_limits<std::size_t>::max());
}

TEMPLATE_TEST_CASE("projected interowner edges of the 7-vertex example", "", std::uint64_t,
                   double) {
  using W = TestType;
  const auto g = CsrGraph<W>::from_edges(7, seven_vertex_edges<W>());
  const auto od = owner_decomposition(g, 2);
  REQUIRE(od.pivots == std::vector<Vertex>{1, 3, 5});
  REQUIRE(od.cut == std::vector<std::uint8_t>{0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0});
  REQUIRE(od.reduced == std::vector<std::uint8_t>{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1});

  const auto pe = projected_interowner_edges(g, od);
  CHECK(pe.instance.t == 3);
  REQUIRE(pe.instance.edges.size() == 5);
  CHECK(pe.input_edge == std::vector<EdgeId>{2, 4, 6, 7, 10});
  const std::vector<W> lengths{W{5}, W{7}, W{15}, W{14}, W{5}};
  const std::vector<std::pair<Vertex, Vertex>> ends{{0, 1}, {1, 2}, {0, 1}, {0, 2}, {0, 1}};
  for (std::size_t i = 0; i < 5; ++i) {
    CHECK(pe.instance.edges[i].weight == lengths[i]);
    CHECK(pe.instance.edges[i].u == ends[i].first);
    CHECK(pe.instance.edges[i].v == ends[i].second);
  }
  const std::uint32_t no = kNoInstanceEdge;
  CHECK(pe.instance_edge == std::vector<std::uint32_t>{no, no, 0, no, 1, no, 2, 3, no, no, 4});

  // The hand names are valid route certificates for P = 1 (Definition 6.1).
  const auto names = seven_vertex_names<W>();
  const auto dist = sssp::test::floyd_warshall<W>(3, pe.instance.edges);
  CHECK(check_route_names(pe.instance, names, dist).empty());
}

TEMPLATE_TEST_CASE("classification of the 7-vertex example from given names", "", std::uint64_t,
                   double) {
  using W = TestType;
  const auto g = CsrGraph<W>::from_edges(7, seven_vertex_edges<W>());
  const auto od = owner_decomposition(g, 2);
  const auto pe = projected_interowner_edges(g, od);
  Diagnostics diag;
  const auto cls = classify_outer_edges(g, od, pe, seven_vertex_names<W>(), &diag);

  const auto R = OuterEdgeClass::kRadial;
  const auto M = OuterEdgeClass::kMacro;
  const auto O = OuterEdgeClass::kOmitted;
  // d_e = P lambda_e = lambda_e. Edge 6: 2 * 5 = 10 <= 10 is radial (the test is not strict).
  // Edges 2, 4, 7, 10: 10 > 4, 14 > 5, 28 > 12, 10 > 5 are macro. Intraowner edges are radial.
  CHECK(cls.edge_class == std::vector<OuterEdgeClass>{R, R, M, R, M, R, R, M, O, O, M});
  CHECK(cls.interowner == std::vector<std::uint8_t>{0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1});
  const std::uint32_t no = kNoName;
  CHECK(cls.name == std::vector<std::uint32_t>{no, no, 0, no, 1, no, 0, 2, no, no, 0});
  CHECK(cls.dictionary == std::vector<W>{W{5}, W{7}, W{14}});
  CHECK(cls.dictionary_detour == std::vector<W>{W{5}, W{7}, W{14}});
  CHECK(cls.detour ==
        std::vector<W>{W{0}, W{0}, W{5}, W{0}, W{7}, W{0}, W{5}, W{14}, W{0}, W{0}, W{5}});
  CHECK(cls.p_log2 == 0);
  CHECK(cls.hop_bound == 3);
  // 5 interowner edges (1 radial, 4 macro) and 2 omitted ones, fixed by the vectors above.
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["aux_metric.interowner_edges"] == 5);
  CHECK(diag.counters["aux_metric.radial_interowner_edges"] == 1);
  CHECK(diag.counters["aux_metric.macro_edges"] == 4);
  CHECK(diag.counters["aux_metric.omitted_edges"] == 2);
#endif
  CHECK(diag.fallbacks.empty());
}

TEMPLATE_TEST_CASE("detours are formed by doubling once per name", "", std::uint64_t, double) {
  using W = TestType;
  const auto g = CsrGraph<W>::from_edges(7, seven_vertex_edges<W>());
  const auto od = owner_decomposition(g, 2);
  const auto pe = projected_interowner_edges(g, od);
  auto names = seven_vertex_names<W>();
  names.p_log2 = 2;  // P = 4: every interowner edge becomes macro (2 * 4 * 5 > 10)
  const auto cls = classify_outer_edges(g, od, pe, names);
  CHECK(cls.dictionary_detour == std::vector<W>{W{20}, W{28}, W{56}});
  CHECK(cls.detour[6] == W{20});
  CHECK(cls.edge_class[6] == OuterEdgeClass::kMacro);
  CHECK(cls.p_log2 == 2);
}

TEST_CASE("integer detours saturate at the infinity sentinel and are then macro") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(7, seven_vertex_edges<W>());
  const auto od = owner_decomposition(g, 2);
  const auto pe = projected_interowner_edges(g, od);
  auto names = seven_vertex_names<W>();
  names.p_log2 = 63;  // 5 * 2^63 does not fit in 64 bits
  Diagnostics diag;
  const auto cls = classify_outer_edges(g, od, pe, names, &diag);
  // All 3 dictionary detours saturate (the size check keeps the loop from passing vacuously).
  REQUIRE(cls.dictionary_detour.size() == 3);
  for (const W d : cls.dictionary_detour) CHECK(d == sssp::WeightTraits<W>::infinity());
  CHECK(cls.edge_class[6] == OuterEdgeClass::kMacro);
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["aux_metric.saturated_detours"] == 3);
#endif
}

TEMPLATE_TEST_CASE("invalid route names are rejected", "", std::uint64_t, double) {
  using W = TestType;
  const auto g = CsrGraph<W>::from_edges(7, seven_vertex_edges<W>());
  const auto od = owner_decomposition(g, 2);
  const auto pe = projected_interowner_edges(g, od);
  {
    auto names = seven_vertex_names<W>();
    names.name.pop_back();  // one name per instance edge is required
    CHECK_THROWS_AS(classify_outer_edges(g, od, pe, names), std::invalid_argument);
  }
  {
    auto names = seven_vertex_names<W>();
    names.name[0] = 3;  // outside the dictionary
    CHECK_THROWS_AS(classify_outer_edges(g, od, pe, names), std::invalid_argument);
  }
  {
    auto names = seven_vertex_names<W>();
    names.name[0] = 1;  // lambda = 7 > A_e = 5 violates Definition 6.1
    CHECK_THROWS_AS(classify_outer_edges(g, od, pe, names), std::invalid_argument);
  }
}

TEMPLATE_TEST_CASE("auxiliary graph, owner-spoke tree and J of the 7-vertex example", "",
                   std::uint64_t, double) {
  using W = TestType;
  const auto g = CsrGraph<W>::from_edges(7, seven_vertex_edges<W>());
  const auto od = owner_decomposition(g, 2);
  const auto pe = projected_interowner_edges(g, od);
  const auto cls = classify_outer_edges(g, od, pe, seven_vertex_names<W>());
  Diagnostics diag;
  const auto aux = build_aux_graph(g, od, cls, &diag);

  // Pivots 1, 3, 5 become the root a = 0; nonpivots 0, 2, 4, 6 become 1, 2, 3, 4.
  CHECK(aux.num_vertices() == 5);
  CHECK(aux.hop_bound == 3);
  CHECK(aux.aux_vertex == std::vector<Vertex>{1, 0, 2, 0, 3, 0, 4});
  CHECK(aux.input_vertex == std::vector<Vertex>{kNoVertex, 0, 2, 4, 6});
  // Radial edges 0, 1, 3, 5 (intraowner, w' = w) and 6 (w' = 10 + 5).
  CHECK(aux.input_edge == std::vector<EdgeId>{0, 1, 3, 5, 6});
  CHECK(aux.aux_edge ==
        std::vector<EdgeId>{0, 1, kNoEdge, 2, kNoEdge, 3, 4, kNoEdge, kNoEdge, kNoEdge, kNoEdge});
  const auto edges = aux.graph.edges();
  REQUIRE(edges.size() == 5);
  const std::vector<std::pair<Vertex, Vertex>> ends{{1, 0}, {0, 2}, {0, 3}, {0, 4}, {1, 3}};
  const std::vector<W> wprime{W{3}, W{1}, W{2}, W{1}, W{15}};
  for (std::size_t i = 0; i < 5; ++i) {
    CHECK(edges[i].u == ends[i].first);
    CHECK(edges[i].v == ends[i].second);
    CHECK(edges[i].weight == wprime[i]);
  }
  // T' = the owner spokes (uncut tree edges), rooted at a.
  CHECK(aux.tree_edges == std::vector<EdgeId>{0, 1, 2, 3});
  CHECK(aux.tree_parent == std::vector<Vertex>{kNoVertex, 0, 0, 0, 0});
  CHECK(aux.tree_parent_edge == std::vector<EdgeId>{kNoEdge, 0, 1, 2, 3});
  CHECK(aux.tree_order.front() == kAuxRoot);
  CHECK(aux.tree_order.size() == 5);
  CHECK(aux.radius == std::vector<W>{W{0}, W{3}, W{1}, W{2}, W{1}});
  CHECK(aux.owner_distance == std::vector<W>{W{0}, W{3}, W{1}, W{2}, W{1}});
  // 5 vertices, 5 edges and no loop at a: every radial edge has an aux edge (aux_edge above).
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["aux_metric.aux_vertices"] == 5);
  CHECK(diag.counters["aux_metric.aux_edges"] == 5);
  CHECK(diag.counters["aux_metric.root_loops"] == 0);
#endif

  const auto J = owner_coordinates_ref(aux);
  CHECK(J == std::vector<W>{W{0}, W{3}, W{1}, W{2}, W{1}});
  CHECK(check_aux_metric(g, od, cls, aux, J).empty());
}

TEMPLATE_TEST_CASE("one owner: J is the distance to the single pivot and can be below R", "",
                   std::uint64_t, double) {
  using W = TestType;
  // 0 -2- 1 -2- 2 and 0 -3- 2: MST = path, edge 2 kept (3 <= D_T = 4). With the paper's outer
  // b = L_3^64 (saturated) the root 0 is the only pivot (no mark is made).
  const std::vector<Edge<W>> e{{0, 1, W{2}}, {1, 2, W{2}}, {0, 2, W{3}}};
  const auto g = CsrGraph<W>::from_edges(3, e);
  const auto od = owner_decomposition(g, outer_separator_size(3));
  REQUIRE(od.pivots == std::vector<Vertex>{0});
  WordRng rng(7);
  Diagnostics diag;
  const auto cls = classify_interowner_edges(g, od, rng, {}, &diag);
  CHECK(cls.edge_class == std::vector<OuterEdgeClass>(3, OuterEdgeClass::kRadial));
  CHECK(cls.detour == std::vector<W>(3, W{0}));
  CHECK(cls.dictionary.empty());
  CHECK(cls.interowner == std::vector<std::uint8_t>(3, 0));
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["aux_metric.interowner_edges"] == 0);
#endif
  const auto aux = build_aux_graph(g, od, cls);
  CHECK(aux.input_vertex == std::vector<Vertex>{kNoVertex, 1, 2});
  CHECK(aux.radius == std::vector<W>{W{0}, W{2}, W{2}});
  CHECK(aux.owner_distance == std::vector<W>{W{0}, W{2}, W{4}});
  const auto J = owner_coordinates_ref(aux);
  CHECK(J == std::vector<W>{W{0}, W{2}, W{3}});  // r_2 = 2 <= J_2 = 3 < R_2 = 4
  CHECK(check_aux_metric(g, od, cls, aux, J).empty());
}

TEMPLATE_TEST_CASE("a radial edge between two pivots becomes a discarded loop at a", "",
                   std::uint64_t, double) {
  using W = TestType;
  // Unit path 0-1-2-3-4 (edges 0..3), a unit detour 0-5-4 (edges 4, 5) and the chord 0 -4- 4
  // (edge 6). All unit edges tie; edge 5 closes the cycle and is the nontree edge, so
  // D_T(0, 4) = 4 and the chord is reduced (4 <= 4). b = 1 marks every vertex: all are pivots,
  // every radius is 0 and A_e = w_e. G_P distance of (0, 4) is 2 (through 5), so with P = 2 the
  // name lambda = 1 certifies the chord (2 * 1 >= 2, 1 <= 4): d = 2, 2d = 4 <= 4 is radial. Both
  // ends are pivots, so in G' the chord is a positive loop at a and is discarded (Sec 7.2).
  const std::vector<Edge<W>> e{{0, 1, W{1}}, {1, 2, W{1}}, {2, 3, W{1}}, {3, 4, W{1}},
                               {0, 5, W{1}}, {5, 4, W{1}}, {0, 4, W{4}}};
  const auto g = CsrGraph<W>::from_edges(6, e);
  const auto od = owner_decomposition(g, 1);
  REQUIRE(od.pivots == std::vector<Vertex>{0, 1, 2, 3, 4, 5});
  REQUIRE(od.tree_edges == std::vector<EdgeId>{0, 1, 2, 3, 4});
  REQUIRE(od.reduced == std::vector<std::uint8_t>{1, 1, 1, 1, 1, 1, 1});
  const auto pe = projected_interowner_edges(g, od);
  REQUIRE(pe.input_edge == std::vector<EdgeId>{0, 1, 2, 3, 4, 5, 6});
  RouteNames<W> names;
  names.dictionary = {W{1}};
  names.name = std::vector<std::uint32_t>(7, 0);
  names.p_log2 = 1;
  const auto dist = sssp::test::floyd_warshall<W>(6, pe.instance.edges);
  REQUIRE(check_route_names(pe.instance, names, dist).empty());

  const auto cls = classify_outer_edges(g, od, pe, names);
  const auto M = OuterEdgeClass::kMacro;
  // Unit edges: d = 2, 2d = 4 > 1. The chord: 2d = 4 <= 4.
  CHECK(cls.edge_class == std::vector<OuterEdgeClass>{M, M, M, M, M, M, OuterEdgeClass::kRadial});
  Diagnostics diag;
  const auto aux = build_aux_graph(g, od, cls, &diag);
  CHECK(aux.num_vertices() == 1);
  CHECK(aux.graph.num_edges() == 0);
  // The radial chord (both ends pivots) is the one discarded loop at a.
  CHECK(aux.aux_edge[6] == kNoEdge);
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["aux_metric.root_loops"] == 1);
#endif
  const auto J = owner_coordinates_ref(aux);
  CHECK(J == std::vector<W>{W{0}});
  CHECK(check_aux_metric(g, od, cls, aux, J).empty());
}

TEMPLATE_TEST_CASE("check_aux_metric reports violations of Lemmas 7.1 and 7.2", "", std::uint64_t,
                   double) {
  using W = TestType;
  const auto g = CsrGraph<W>::from_edges(7, seven_vertex_edges<W>());
  const auto od = owner_decomposition(g, 2);
  const auto pe = projected_interowner_edges(g, od);
  const auto cls = classify_outer_edges(g, od, pe, seven_vertex_names<W>());
  const auto aux = build_aux_graph(g, od, cls);
  const auto J = owner_coordinates_ref(aux);
  {
    auto bad = J;
    bad[1] = W{2};  // below r = 3
    CHECK_FALSE(check_aux_metric(g, od, cls, aux, bad).empty());
  }
  {
    auto bad = J;
    bad[3] = W{17};  // above R = 2
    CHECK_FALSE(check_aux_metric(g, od, cls, aux, bad).empty());
  }
  {
    auto bad = cls;
    bad.edge_class[2] = OuterEdgeClass::kRadial;  // a cut corridor edge must be macro
    CHECK_FALSE(check_aux_metric(g, od, bad, aux, J).empty());
  }
  {
    auto bad = cls;
    bad.detour[6] = W{6};  // radial with 2d = 12 > w = 10
    CHECK_FALSE(check_aux_metric(g, od, bad, aux, J).empty());
  }
}
