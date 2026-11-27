// Tests for the threshold components and joins of Sec 5.3 (Lemmas 5.8 and 5.9) on two
// hand-computed local state graphs (owner decomposition with b = 2, paper K and macro test).
// Expected values were derived by hand from the definitions of I0, O0, I*, O* (l.1110-1117).

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::Vertex;
using sssp::paper::compute_joins;
using sssp::paper::JoinKind;
using sssp::paper::StateId;

namespace {

template <class W>
constexpr W kInf = sssp::WeightTraits<W>::infinity();

}  // namespace

TEMPLATE_TEST_CASE("joins: pivots joined by a corridor, radial copies stay local", "",
                   std::uint64_t, double) {
  using W = TestType;
  // Pivots 0 and 1 (b = 2); cells {0, 2}, {1, 3}; edge 0 is the corridor (macro), 3 a near
  // radial edge between the two positive vertices (r = 100).
  const std::vector<Edge<W>> edges{{0, 1, W{1}}, {0, 2, W{100}}, {1, 3, W{100}}, {2, 3, W{200}}};
  const auto g = CsrGraph<W>::from_edges(4, edges);
  const auto od = sssp::paper::owner_decomposition(g, 2);
  const auto rs = sssp::paper::local_radial_states(g, od);
  REQUIRE(rs.num_states() == 10);
  const StateId in0 = rs.native_inward[0], out0 = rs.native_outward[0];
  const StateId in1 = rs.native_inward[1], out1 = rs.native_outward[1];
  const StateId in2 = rs.native_inward[2], out2 = rs.native_outward[2], up2 = rs.upper_outward[2];
  const StateId in3 = rs.native_inward[3], out3 = rs.native_outward[3], up3 = rs.upper_outward[3];

  sssp::Diagnostics diag;
  const auto j = compute_joins(od, rs, &diag);
  const W inf = kInf<W>;
  const auto at = [&](const std::vector<W>& v, StateId s) { return v[s]; };
  // I0: reachable from the head resets in0, in1 (weight 1).
  for (const StateId s : {in0, out0, in1, out1, out2, out3}) CHECK(at(j.i0, s) == W{1});
  for (const StateId s : {in2, in3, up2, up3}) CHECK(at(j.i0, s) == inf);
  // O0: reaching a tail copy of the corridor (copies of 0 and 1).
  for (const StateId s : {in0, out0, in1, out1, in2, in3}) CHECK(at(j.o0, s) == W{1});
  for (const StateId s : {out2, out3, up2, up3}) CHECK(at(j.o0, s) == inf);
  // Bottlenecks: out2, out3 are reached from the core through the pivot arcs (100).
  for (const StateId s : {in0, out0, in1, out1}) CHECK(at(j.i_star, s) == W{1});
  CHECK(at(j.i_star, out2) == W{100});
  CHECK(at(j.i_star, out3) == W{100});
  CHECK(at(j.i_star, in2) == inf);
  CHECK(at(j.o_star, in2) == W{100});
  CHECK(at(j.o_star, in3) == W{100});
  CHECK(at(j.o_star, out2) == inf);
  for (const StateId s : {in0, out0, in1, out1}) {
    CHECK(at(j.theta, s) == W{1});
    CHECK(at(j.b0, s) == W{1});
    CHECK(j.b0_edge[s] == 0);
    CHECK(j.b0_key[s] == od.corridor_max(0, 1));
  }
  for (const StateId s : {in2, out2, up2, in3, out3, up3}) {
    CHECK(at(j.theta, s) == inf);
    CHECK(at(j.b0, s) == inf);
    CHECK(j.b0_edge[s] == kNoEdge);
  }

  // Join records: 3 block MST edges (weight 200), 4 owner joins (weight 1), 1 corridor join.
  CHECK(j.num_join_vertices() == 12);
  CHECK(j.aux_vertex(0) == 10);
  CHECK(j.aux_vertex(1) == 11);
  std::size_t mst = 0, own = 0, cor = 0;
  for (const auto& r : j.joins) {
    switch (r.kind) {
      case JoinKind::kBlockMst:
        ++mst;
        CHECK(r.weight == W{200});
        CHECK(rs.block[r.a] == rs.block[r.b]);
        break;
      case JoinKind::kOwner:
        ++own;
        CHECK(r.weight == W{1});
        CHECK(r.b == j.aux_vertex(od.owner[rs.projection[r.a]]));
        break;
      case JoinKind::kCorridor:
        ++cor;
        CHECK(r.weight == W{1});
        CHECK(((r.a == 10 && r.b == 11) || (r.a == 11 && r.b == 10)));
        break;
    }
  }
  CHECK(mst == 3);
  CHECK(own == 4);
  CHECK(cor == 1);
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["threshold_joins_calls"] == 1);
  CHECK(diag.counters["threshold_joins_records"] == 8);
#endif
}

TEMPLATE_TEST_CASE("joins: a nontree macro edge gives positive copies finite theta", "",
                   std::uint64_t, double) {
  using W = TestType;
  // Path 0 -1- 1 -5- 2 -3- 3 plus nontree edge 1-3 (6 <= D_T = 8). b = 2: pivots 0, 2; the
  // corridor 0-1-2 is cut at edge 1 (weight 5); edge 3 is intercell with 6 <= N^3 * 5: macro.
  const std::vector<Edge<W>> edges{{0, 1, W{1}}, {1, 2, W{5}}, {2, 3, W{3}}, {1, 3, W{6}}};
  const auto g = CsrGraph<W>::from_edges(4, edges);
  const auto od = sssp::paper::owner_decomposition(g, 2);
  REQUIRE(od.pivots == std::vector<Vertex>{0, 2});
  const auto rs = sssp::paper::local_radial_states(g, od);
  REQUIRE(rs.num_states() == 10);
  REQUIRE(rs.num_blocks == 10);
  const StateId in0 = rs.native_inward[0], out0 = rs.native_outward[0];
  const StateId in1 = rs.native_inward[1], out1 = rs.native_outward[1], up1 = rs.upper_outward[1];
  const StateId in2 = rs.native_inward[2], out2 = rs.native_outward[2];
  const StateId in3 = rs.native_inward[3], out3 = rs.native_outward[3], up3 = rs.upper_outward[3];
  const auto j = compute_joins(od, rs);
  const W inf = kInf<W>;

  const std::vector<StateId> order{in0, out0, in1, out1, up1, in2, out2, in3, out3, up3};
  const auto pick = [&](const std::vector<W>& v) {
    std::vector<W> out;
    for (const StateId s : order) out.push_back(v[s]);
    return out;
  };
  const W f{5}, s{6};
  CHECK(pick(j.i0) == std::vector<W>{f, f, f, f, inf, f, f, s, f, inf});
  CHECK(pick(j.o0) == std::vector<W>{f, f, f, f, f, f, f, f, s, s});
  CHECK(pick(j.i_star) == std::vector<W>{f, f, f, f, inf, f, f, s, f, inf});
  CHECK(pick(j.o_star) == std::vector<W>{f, f, f, f, f, f, f, f, s, s});
  CHECK(pick(j.theta) == std::vector<W>{f, f, f, f, inf, f, f, s, s, inf});
  CHECK(pick(j.b0) == std::vector<W>{f, f, f, f, inf, f, f, s, s, inf});
  // Attaining macro edges: (5, edge 1) wins ties with (6, edge 3) by weight; B0 at in3 and out3
  // is attained by edge 3.
  CHECK(j.i0_edge[in0] == 1);
  CHECK(j.i0_edge[in3] == 3);
  CHECK(j.o0_edge[out3] == 3);
  CHECK(j.b0_edge[in3] == 3);
  CHECK(j.b0_edge[out3] == 3);
  CHECK(j.b0_edge[in1] == 1);
  CHECK(j.i0_edge[up1] == kNoEdge);

  // 8 owner joins (finite theta) and one corridor join; no block MST edges.
  CHECK(j.joins.size() == 9);
}

TEST_CASE("joins reject a state graph of a different decomposition") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(3, {{0, 1, 1}, {1, 2, 2}});
  const auto g2 = CsrGraph<W>::from_edges(2, {{0, 1, 1}});
  const auto od = sssp::paper::owner_decomposition(g, 1);
  const auto rs2 = sssp::paper::local_radial_states(g2, sssp::paper::owner_decomposition(g2, 1));
  CHECK_THROWS_AS(compute_joins(od, rs2), std::invalid_argument);
}
