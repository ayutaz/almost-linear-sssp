// The outer state graph of arXiv 2611.34911v1, Sec 8.1 (txt l.1966-2001; Lemmas 8.1-8.3):
// deletion of radial edges with w > J_x + J_y + d_e, compared in the calendar units of Alg 7 as
// 2^7 w > J'_x + J'_y + 2^7 d_e (J' = 2^7 J, D-21), the two-shell construction on J with K = 16,
// and the directed owner-tree states v_in / v_out (singleton blocks).
//
// Expected values are derived by hand from the paper's rules, never copied from a run.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include "paper/local/radial_states.hpp"
#include "paper/outer/outer_states.hpp"
#include "paper/owner/aux_metric.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "unit/paper/outer/outer_states_test_support.hpp"

using sssp::Diagnostics;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::Vertex;
using sssp::paper::EdgeClass;
using sssp::paper::EdgeClassification;
using sssp::paper::kNoState;
using sssp::paper::RadialArcKind;
using sssp::paper::StateId;
using sssp::paper::outer::build_outer_states;
using sssp::paper::outer::kOuterShellRatioLog2;
using sssp::paper::outer::outer_radial_deleted;
using sssp::paper::outer::OuterArcKind;
namespace to = sssp::test::outer;

namespace {

/// Hand example (Mode A: b = L_5^64 saturates, one pivot, vertex 0 = the root of T, Q-31).
///   e0 0-4 w1, e1 2-3 w1, e2 0-1 w2, e3 1-2 w2, e4 0-2 w3, e5 3-4 w6, e6 1-4 w3.
/// T = {e0, e1, e2, e3} (Kruskal; e4, e5, e6 close cycles). Every nontree edge has
/// w <= D_T (D_T(0,2) = 4, D_T(3,4) = 6, D_T(1,4) = 3), so all are reduced, intraowner, radial
/// with d_e = 0, and G' = G. J = (0, 2, 3, 4, 1): J_1 via e2, J_2 via e4, J_3 = J_2 + 1,
/// J_4 via e0. Deletion w > J_x + J_y: e5 (6 > 4 + 1) only; e6 is the equality case 3 = 2 + 1.
template <class W>
std::vector<Edge<W>> hand_edges() {
  return {{0, 4, W{1}}, {2, 3, W{1}}, {0, 1, W{2}}, {1, 2, W{2}},
          {0, 2, W{3}}, {3, 4, W{6}}, {1, 4, W{3}}};
}

using ArcRow = std::tuple<StateId, StateId, double, EdgeId, OuterArcKind>;

}  // namespace

TEMPLATE_TEST_CASE("outer_radial_deleted: 2^7 w > J'_x + J'_y + 2^7 d_e, exact",
                   "[outer_states][unit]", std::uint64_t, double) {
  using W = TestType;
  // J in input units 4 and 1, scaled by 2^7.
  CHECK(outer_radial_deleted<W>(W{6}, W{0}, W{512}, W{128}));
  CHECK_FALSE(outer_radial_deleted<W>(W{5}, W{0}, W{512}, W{128}));  // equality: kept
  CHECK_FALSE(outer_radial_deleted<W>(W{3}, W{0}, W{256}, W{128}));  // equality: kept
  CHECK(outer_radial_deleted<W>(W{4}, W{1}, W{128}, W{128}));        // 512 > 256 + 128
  CHECK_FALSE(outer_radial_deleted<W>(W{4}, W{2}, W{128}, W{128}));  // 512 = 256 + 256
  if constexpr (std::is_floating_point_v<W>) {
    // Quarter-integers (Q-55 plan A): 2^7 2.25 = 288 = 192 + 64 + 2^7 0.25.
    CHECK_FALSE(outer_radial_deleted<W>(W{2.25}, W{0.25}, W{192}, W{64}));
    CHECK(outer_radial_deleted<W>(W{2.5}, W{0.25}, W{192}, W{64}));
    // 2^7 w is not finite: the comparison cannot be made exactly.
    CHECK_THROWS_AS(outer_radial_deleted<W>(std::numeric_limits<W>::max(), W{0}, W{1}, W{1}),
                    std::overflow_error);
  } else {
    // No wrap-around: 2^7 w and the sums exceed 64 bits but are compared exactly.
    const W big = W{1} << 62;
    CHECK(outer_radial_deleted<W>(big, big / 2, W{0}, W{0}));      // 2^69 > 2^68
    const W j = (W{1} << 63) + (W{1} << 62);                       // J'_x + J'_y = 2^64 + 2^63
    CHECK_FALSE(outer_radial_deleted<W>(W{1} << 57, W{0}, j, j));  // 2^64 < 2^65 + 2^64
    CHECK(outer_radial_deleted<W>(W{3} << 57, W{0}, j, j - 1));    // 3 2^64 > 3 2^64 - 1
  }
}

TEMPLATE_TEST_CASE("build_outer_states: hand example in Mode A (one pivot)", "[outer_states][unit]",
                   std::uint64_t, double) {
  using W = TestType;
  const auto edges = hand_edges<W>();
  const std::size_t n = 5;
  auto c = to::make_outer_case<W>(n, edges, sssp::paper::owner::outer_separator_size(n),
                                  to::Names::kExactP1, 1);
  REQUIRE(c.od.pivots == std::vector<Vertex>{0});
  // Sanity of the inputs (Lemma 7.2, by hand above).
  const std::vector<W> j_expected{W{0}, W{2}, W{3}, W{4}, W{1}};
  for (Vertex v = 0; v < n; ++v) {
    REQUIRE(c.jr.scaled_distance[c.aux.aux_vertex[v]] == W{128} * j_expected[v]);
  }

  const std::size_t fallbacks_before = c.diag.fallbacks.size();
  Diagnostics diag;
  const auto os = build_outer_states(c.g, c.od, c.cls, c.aux, c.jr.scaled_distance, &diag);
  CHECK(diag.fallbacks.empty());  // Sec 8.1 itself never substitutes a search
  CHECK(c.diag.fallbacks.size() == fallbacks_before);

  // Step 1: deletion.
  CHECK(os.deleted == std::vector<std::uint8_t>{0, 0, 0, 0, 0, 1, 0});
  CHECK(os.num_deleted == 1);

  // Step 2: the ordinary part is Def 5.4 on (G, J', classification after deletion, K = 16).
  CHECK(kOuterShellRatioLog2 == 4);
  EdgeClassification expected_cls{std::vector<EdgeClass>(7, EdgeClass::kRadial),
                                  std::vector<std::uint32_t>(7, sssp::paper::kNoCorridor)};
  expected_cls.edge_class[5] = EdgeClass::kOmitted;
  std::vector<W> coord;
  for (const W x : j_expected) coord.push_back(W{128} * x);
  const auto ref = sssp::paper::build_radial_states(c.g, coord, expected_cls, 4);
  CHECK(os.states.shell_ratio_log2 == 4);
  CHECK(os.states.coordinate == coord);
  CHECK(os.states.classification.edge_class == expected_cls.edge_class);
  CHECK(os.states.classification.macro_key == expected_cls.macro_key);
  REQUIRE(os.num_ordinary_states == ref.num_states());
  REQUIRE(os.num_ordinary_arcs == ref.graph.arcs.size());
  CHECK(os.num_ordinary_blocks == ref.num_blocks);
  for (std::size_t a = 0; a < ref.graph.arcs.size(); ++a) {
    CHECK(os.states.graph.arcs[a].tail == ref.graph.arcs[a].tail);
    CHECK(os.states.graph.arcs[a].head == ref.graph.arcs[a].head);
    CHECK(os.states.graph.arcs[a].length == ref.graph.arcs[a].length);
    CHECK(os.states.arc_kind[a] == ref.arc_kind[a]);
    CHECK(os.arc_kind[a] == OuterArcKind::kOrdinary);
  }
  for (StateId x = 0; x < ref.num_states(); ++x) {
    CHECK(os.states.projection[x] == ref.projection[x]);
    CHECK(os.states.phase[x] == ref.phase[x]);
    CHECK(os.states.block[x] == ref.block[x]);
  }
  CHECK(os.states.native_inward == ref.native_inward);
  CHECK(os.states.native_outward == ref.native_outward);
  // J in 1..4 with anchor J_4 = 1: one near component, every positive vertex in shell 0.
  CHECK(os.states.shell == std::vector<std::uint32_t>(5, 0));

  // Step 3: the special states and their arcs (l.1984-1997).
  REQUIRE(os.special_inward[0] == kNoState);
  REQUIRE(os.special_outward[0] == kNoState);
  CHECK(os.states.num_states() == ref.num_states() + 8);
  const std::vector<Vertex> parent{sssp::kNoVertex, 0, 1, 2, 0};
  const std::vector<EdgeId> parent_edge{kNoEdge, 2, 3, 1, 0};
  CHECK(os.owner_parent == parent);
  CHECK(os.owner_parent_edge == parent_edge);
  std::vector<ArcRow> expected;
  const auto in = [&](Vertex v) {
    return v == 0 ? os.states.native_inward[0] : os.special_inward[v];
  };
  const auto out = [&](Vertex v) {
    return v == 0 ? os.states.native_outward[0] : os.special_outward[v];
  };
  for (Vertex v = 1; v < n; ++v) {
    REQUIRE(os.special_inward[v] != kNoState);
    REQUIRE(os.special_outward[v] != kNoState);
    CHECK(os.states.projection[os.special_inward[v]] == v);
    CHECK(os.states.projection[os.special_outward[v]] == v);
    CHECK(os.states.phase[os.special_inward[v]] == sssp::paper::Phase::kInward);
    CHECK(os.states.phase[os.special_outward[v]] == sssp::paper::Phase::kOutward);
    // Shell 0: one ordinary inward copy, the native one.
    REQUIRE(os.states.lower_inward[v] == kNoState);
    const double w = static_cast<double>(edges[parent_edge[v]].weight);
    expected.emplace_back(os.states.native_inward[v], os.special_inward[v], 0.0, kNoEdge,
                          OuterArcKind::kTreeEntry);
    expected.emplace_back(os.special_inward[v], in(parent[v]), w, parent_edge[v],
                          OuterArcKind::kInwardTree);
    expected.emplace_back(out(parent[v]), os.special_outward[v], w, parent_edge[v],
                          OuterArcKind::kOutwardTree);
    expected.emplace_back(os.special_outward[v], os.states.native_outward[v], 0.0, kNoEdge,
                          OuterArcKind::kTreeExit);
  }
  std::vector<ArcRow> actual;
  for (std::size_t a = ref.graph.arcs.size(); a < os.states.graph.arcs.size(); ++a) {
    const auto& arc = os.states.graph.arcs[a];
    CHECK(os.states.arc_kind[a] == RadialArcKind::kOwnerTree);
    actual.emplace_back(arc.tail, arc.head, static_cast<double>(arc.length), os.states.arc_edge[a],
                        os.arc_kind[a]);
  }
  std::ranges::sort(expected);
  std::ranges::sort(actual);
  CHECK(actual == expected);

  // Special states are singleton blocks (l.1998).
  CHECK(os.states.num_blocks == ref.num_blocks + 8);
  for (Vertex v = 1; v < n; ++v) {
    for (const StateId x : {os.special_inward[v], os.special_outward[v]}) {
      CHECK(os.states.block[x] >= ref.num_blocks);
      for (StateId y = 0; y < os.states.num_states(); ++y) {
        if (y != x) CHECK(os.states.block[y] != os.states.block[x]);
      }
    }
  }

  // Step 4: no macro edges with one pivot, so every theta is infinite and the joins are the
  // block MST edges only (no owner joins, no corridors).
  CHECK(os.joins.num_states == os.states.num_states());
  CHECK(os.joins.num_pivots == 1);
  for (StateId z = 0; z < os.joins.num_states; ++z) {
    CHECK(os.joins.theta[z] == sssp::WeightTraits<W>::infinity());
  }
  for (const auto& jr : os.joins.joins) CHECK(jr.kind == sssp::paper::JoinKind::kBlockMst);
  CHECK(os.hop_bound == c.od.hop_bound);
}

TEMPLATE_TEST_CASE("build_outer_states: rejects inconsistent inputs", "[outer_states][unit]",
                   std::uint64_t, double) {
  using W = TestType;
  const auto edges = hand_edges<W>();
  auto c = to::make_outer_case<W>(5, edges, sssp::paper::owner::outer_separator_size(5),
                                  to::Names::kExactP1, 1);
  auto j = c.jr.scaled_distance;
  {
    auto bad = j;
    bad.pop_back();
    CHECK_THROWS_AS(build_outer_states(c.g, c.od, c.cls, c.aux, bad), std::invalid_argument);
  }
  {
    auto bad = j;
    bad[sssp::paper::owner::kAuxRoot] = W{1};  // J_a = 0 (Lemma 7.2)
    CHECK_THROWS_AS(build_outer_states(c.g, c.od, c.cls, c.aux, bad), std::invalid_argument);
  }
  {
    auto bad = j;
    bad[c.aux.aux_vertex[3]] = W{0};  // a nonpivot has J >= r > 0 (Lemma 7.2)
    CHECK_THROWS_AS(build_outer_states(c.g, c.od, c.cls, c.aux, bad), std::invalid_argument);
  }
  {
    auto bad = j;
    bad[c.aux.aux_vertex[3]] = sssp::WeightTraits<W>::infinity();
    CHECK_THROWS_AS(build_outer_states(c.g, c.od, c.cls, c.aux, bad), std::invalid_argument);
  }
  {
    auto cls = c.cls;
    cls.edge_class.pop_back();
    CHECK_THROWS_AS(build_outer_states(c.g, c.od, cls, c.aux, j), std::invalid_argument);
  }
  {
    const auto other = sssp::CsrGraph<W>::from_edges(6, edges);
    CHECK_THROWS_AS(build_outer_states(other, c.od, c.cls, c.aux, j), std::invalid_argument);
  }
}

TEMPLATE_TEST_CASE("build_outer_states: a single vertex has no radial edges and no specials",
                   "[outer_states][unit]", std::uint64_t, double) {
  using W = TestType;
  auto c = to::make_outer_case<W>(1, {}, sssp::paper::owner::outer_separator_size(1),
                                  to::Names::kExactP1, 1);
  const auto os = build_outer_states(c.g, c.od, c.cls, c.aux, c.jr.scaled_distance);
  CHECK(os.num_deleted == 0);
  CHECK(os.states.num_states() == 2);  // the pivot copies in, out
  CHECK(os.num_ordinary_states == 2);
  CHECK(os.special_inward == std::vector<StateId>{kNoState});
  CHECK(os.joins.joins.empty());
}

TEMPLATE_TEST_CASE("build_outer_states: the owner-tree states join the cores (Lemma 8.3)",
                   "[outer_states][unit]", std::uint64_t, double) {
  using W = TestType;
  for (const std::size_t k : {std::size_t{34}, std::size_t{36}, std::size_t{45}}) {
    for (const auto names : {to::Names::kExactP1, to::Names::kNameRoutesPaper}) {
      INFO("k = " << k << (names == to::Names::kExactP1 ? " exact-P1" : " name-routes"));
      to::ThreePivotCorridor r;
      const auto edges = to::three_pivot_corridor<W>(k, r);
      auto c = to::make_outer_case<W>(r.n, edges, r.b, names, k);
      // Preconditions of the construction (outer_states_test_support.hpp), by hand.
      std::vector<Vertex> pivots{r.p, r.q, r.p2};
      std::ranges::sort(pivots);
      REQUIRE(c.od.pivots == pivots);
      REQUIRE(c.od.cut[r.corridor_yp] == 1);
      REQUIRE(c.od.cut[r.corridor_qp2] == 1);
      REQUIRE(c.cls.edge_class[r.corridor_yp] == sssp::paper::owner::OuterEdgeClass::kMacro);
      REQUIRE(c.cls.edge_class[r.corridor_qp2] == sssp::paper::owner::OuterEdgeClass::kMacro);
      REQUIRE(c.jr.scaled_distance[c.aux.aux_vertex[r.y]] == W{128 * 260});
      REQUIRE(c.jr.scaled_distance[c.aux.aux_vertex[r.t]] == W{128 * 17});

      const auto os = build_outer_states(c.g, c.od, c.cls, c.aux, c.jr.scaled_distance);
      const auto& st = os.states;
      CHECK(os.deleted[r.shortcut_yq] == 0);  // 260 = J_y + J_q + 0: kept (l.1967)
      CHECK(st.shell[r.y] == 0);
      CHECK(st.shell[r.path[0]] == 0);
      CHECK(st.shell[r.path[k / 2]] == 2);

      const StateId in_p = st.native_inward[r.p];
      const StateId in_q = st.native_inward[r.q];
      const StateId in_p2 = st.native_inward[r.p2];
      const StateId in_y = st.native_inward[r.y];
      const auto scc_at = [&](W threshold, bool with_special) {
        std::vector<std::pair<std::uint32_t, std::uint32_t>> sel;
        const std::size_t arcs = with_special ? st.graph.arcs.size() : os.num_ordinary_arcs;
        for (std::size_t a = 0; a < arcs; ++a) {
          const auto& arc = st.graph.arcs[a];
          if (!(threshold < arc.length)) sel.emplace_back(arc.tail, arc.head);
        }
        return to::scc_ids(st.num_states(), sel);
      };
      const auto joins_at = [&](W threshold) {
        to::UnionFind uf(os.joins.num_join_vertices());
        for (const auto& j : os.joins.joins) {
          if (!(threshold < j.weight)) uf.unite(j.a, j.b);
        }
        return uf;
      };
      // W = 256: one SCC through the three cores with the special states, two without them.
      const auto full = scc_at(W{256}, true);
      CHECK(full[in_q] == full[in_y]);
      CHECK(full[in_q] == full[in_p]);
      CHECK(full[in_q] == full[in_p2]);
      CHECK(full[in_q] == full[os.special_inward[r.y]]);
      CHECK(full[in_q] == full[os.special_outward[r.y]]);
      const auto ordinary = scc_at(W{256}, false);
      CHECK(ordinary[in_y] == ordinary[in_p]);   // y - p macro cycle
      CHECK(ordinary[in_q] == ordinary[in_p2]);  // q - p2 macro cycle
      CHECK(ordinary[in_q] != ordinary[in_y]);   // separated without v_in / v_out
      auto uf256 = joins_at(W{256});
      CHECK(uf256.find(in_q) == uf256.find(in_y));
      CHECK(uf256.find(in_q) == uf256.find(in_p));
      // W = 255: the y - p macro (256) is out; the q - p2 cycle (200) stays.
      const auto below = scc_at(W{255}, true);
      CHECK(below[in_q] == below[in_p2]);
      CHECK(below[in_q] != below[in_y]);
      CHECK(below[in_p] != below[in_y]);
      auto uf255 = joins_at(W{255});
      CHECK(uf255.find(in_q) == uf255.find(in_p2));
      CHECK(uf255.find(in_q) != uf255.find(in_y));
    }
  }
}
