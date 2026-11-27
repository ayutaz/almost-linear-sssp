// Tests for the radial state representation of Sec 5.2 (Def 5.4, Lemma 5.7) on hand-built
// inputs. The generic builder takes the coordinate, the shell ratio K = 2^k and the
// radial/macro classification as parameters (Def 5.4, l.1057-1060), so the example below uses
// K = 2 to keep the shells small; the local rule (K >= 8h, macro iff w <= N^3 sigma) is tested
// through classify_local_edges / local_radial_states.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
using sssp::paper::build_radial_states;
using sssp::paper::classify_local_edges;
using sssp::paper::EdgeClass;
using sssp::paper::EdgeClassification;
using sssp::paper::kNoCorridor;
using sssp::paper::kNoNearComponent;
using sssp::paper::kNoState;
using sssp::paper::local_radial_states;
using sssp::paper::local_shell_ratio_log2;
using sssp::paper::Phase;
using sssp::paper::RadialArcKind;
using sssp::paper::StateId;

namespace {

using ArcTuple = std::tuple<StateId, StateId, std::uint64_t, RadialArcKind>;

template <class W>
std::vector<ArcTuple> arc_tuples(const sssp::paper::RadialStateGraph<W>& rs) {
  std::vector<ArcTuple> out;
  for (std::size_t i = 0; i < rs.graph.arcs.size(); ++i) {
    const auto& a = rs.graph.arcs[i];
    out.emplace_back(a.tail, a.head, static_cast<std::uint64_t>(a.length), rs.arc_kind[i]);
  }
  std::sort(out.begin(), out.end());
  return out;
}

}  // namespace

TEST_CASE("local shell ratio K is the least power of two >= 8h") {
  CHECK(local_shell_ratio_log2(1) == 3);
  CHECK(local_shell_ratio_log2(2) == 4);
  CHECK(local_shell_ratio_log2(3) == 5);  // 24 -> 32
  CHECK(local_shell_ratio_log2(4) == 5);
  CHECK(local_shell_ratio_log2(104977) == 20);  // h = L^4 + 1 at the first local size
}

TEMPLATE_TEST_CASE("radial states of a hand example with K = 2", "", std::uint64_t, double) {
  using W = TestType;
  // Vertex 0 has coordinate 0 (a pivot); 1 and 2 are near (ratio 2 <= K^2 = 4), 2 and 3 far
  // (ratio 8 > 4). Edge 3 is a macro edge.
  const std::vector<Edge<W>> edges{{0, 1, W{1}}, {1, 2, W{2}}, {2, 3, W{16}}, {0, 3, W{20}}};
  const auto g = CsrGraph<W>::from_edges(4, edges);
  const std::vector<W> coord{W{0}, W{1}, W{2}, W{16}};
  EdgeClassification cls{
      {EdgeClass::kRadial, EdgeClass::kRadial, EdgeClass::kRadial, EdgeClass::kMacro},
      {kNoCorridor, kNoCorridor, kNoCorridor, 7}};
  const auto rs = build_radial_states(g, coord, cls, 1);

  // Near components {1, 2} (anchor 1) and {3} (anchor 16); shells i(1) = 0, i(2) = 1, i(3) = 0.
  CHECK(rs.shell == std::vector<std::uint32_t>{0, 0, 1, 0});
  CHECK(rs.near_component[0] == kNoNearComponent);
  CHECK(rs.near_component[1] == rs.near_component[2]);
  CHECK(rs.near_component[3] != rs.near_component[1]);
  CHECK(rs.component_anchor[rs.near_component[1]] == W{1});
  CHECK(rs.component_anchor[rs.near_component[3]] == W{16});

  // Copies: pivot in/out; shell-0 vertices have no lower inward copy.
  CHECK(rs.graph.num_states == 12);
  CHECK(rs.lower_inward[0] == kNoState);
  CHECK(rs.upper_outward[0] == kNoState);
  CHECK(rs.lower_inward[1] == kNoState);
  CHECK(rs.lower_inward[3] == kNoState);
  const StateId I0 = rs.native_inward[0], O0 = rs.native_outward[0];
  const StateId in1 = rs.native_inward[1], out1 = rs.native_outward[1], up1 = rs.upper_outward[1];
  const StateId in2 = rs.native_inward[2], lo2 = rs.lower_inward[2], out2 = rs.native_outward[2],
                up2 = rs.upper_outward[2];
  const StateId in3 = rs.native_inward[3], out3 = rs.native_outward[3], up3 = rs.upper_outward[3];
  for (const StateId s : {I0, in1, in2, lo2, in3}) CHECK(rs.phase[s] == Phase::kInward);
  for (const StateId s : {O0, out1, up1, out2, up2, out3, up3})
    CHECK(rs.phase[s] == Phase::kOutward);
  CHECK(rs.window[lo2] == 0);
  CHECK(rs.window[in2] == 1);
  CHECK(rs.window[up2] == 2);
  CHECK(rs.projection[up2] == 2);
  CHECK(rs.projection[I0] == 0);

  using K = RadialArcKind;
  std::vector<ArcTuple> want{
      // pivot radial edge 0-1 (Q-32: from every inward copy of 1)
      {in1, I0, 1, K::kPivot},
      {O0, out1, 1, K::kPivot},
      // near edge 1-2
      {in1, lo2, 2, K::kWithinWindow},
      {lo2, in1, 2, K::kWithinWindow},
      {up1, out2, 2, K::kWithinWindow},
      {out2, up1, 2, K::kWithinWindow},
      {out1, out2, 2, K::kNearTransfer},
      {in2, in1, 2, K::kNearTransfer},
      // far edge 2-3: inward copies of the higher end, outward copies of the lower end
      {in3, in2, 16, K::kFarTransfer},
      {out2, out3, 16, K::kFarTransfer},
      {up2, out3, 16, K::kFarTransfer},
      // macro edge 0-3: every copy of the tail to the native inward copy of the head
      {I0, in3, 20, K::kMacro},
      {O0, in3, 20, K::kMacro},
      {in3, I0, 20, K::kMacro},
      {out3, I0, 20, K::kMacro},
      {up3, I0, 20, K::kMacro},
      // phase switches
      {I0, O0, 0, K::kSwitch},
      {in1, out1, 0, K::kSwitch},
      {in2, out2, 0, K::kSwitch},
      {lo2, out2, 0, K::kSwitch},
      {in3, out3, 0, K::kSwitch},
  };
  std::sort(want.begin(), want.end());
  CHECK(arc_tuples(rs) == want);

  // Blocks: {in1, lo2} (inward window 0), {up1, out2} (outward window 1), the rest singletons.
  CHECK(rs.num_blocks == 10);
  CHECK(rs.block[in1] == rs.block[lo2]);
  CHECK(rs.block[up1] == rs.block[out2]);
  CHECK(rs.block[in1] != rs.block[up1]);
  CHECK(rs.block[out1] != rs.block[up1]);

  // The macro key and arc edges are carried along.
  CHECK(rs.classification.macro_key[3] == 7);
  for (std::size_t i = 0; i < rs.graph.arcs.size(); ++i) {
    if (rs.arc_kind[i] == K::kSwitch) {
      CHECK(rs.arc_edge[i] == sssp::kNoEdge);
    } else {
      CHECK(edges[rs.arc_edge[i]].weight == rs.graph.arcs[i].length);
    }
  }
}

TEST_CASE("shell powers saturate instead of overflowing") {
  using W = std::uint64_t;
  // Radii 2^62 and 2^63 are near for K = 2 (ratio 2 <= 4); 2^62 * K = 2^63 <= 2^63 < 2^64.
  const W a = W{1} << 62;
  const auto g = CsrGraph<W>::from_edges(3, {{0, 1, a}, {1, 2, a}});
  const std::vector<W> coord{0, a, W{1} << 63};
  EdgeClassification cls{{EdgeClass::kRadial, EdgeClass::kRadial}, {kNoCorridor, kNoCorridor}};
  const auto rs = build_radial_states(g, coord, cls, 1);
  CHECK(rs.shell == std::vector<std::uint32_t>{0, 0, 1});
  CHECK(rs.near_component[1] == rs.near_component[2]);
}

TEST_CASE("omitted edges produce no arcs") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(2, {{0, 1, 3}, {0, 1, 9}, {1, 1, 2}});
  EdgeClassification cls{{EdgeClass::kMacro, EdgeClass::kOmitted, EdgeClass::kOmitted},
                         {0, kNoCorridor, kNoCorridor}};
  const auto rs = build_radial_states(g, std::vector<W>{0, 0}, cls, 3);
  CHECK(rs.graph.num_states == 4);
  std::size_t macro = 0;
  for (std::size_t i = 0; i < rs.graph.arcs.size(); ++i) {
    CHECK(rs.arc_edge[i] != 1);
    CHECK(rs.arc_edge[i] != 2);
    macro += rs.arc_kind[i] == RadialArcKind::kMacro ? 1 : 0;
  }
  CHECK(macro == 4);  // two copies per tail, both directions
  CHECK(rs.graph.arcs.size() == 6);
}

TEST_CASE("radial state builder rejects inputs outside Def 5.4") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(2, {{0, 1, 3}});
  const EdgeClassification radial{{EdgeClass::kRadial}, {kNoCorridor}};
  // A radial edge between two zero coordinates (no positive radial edge joins two pivots).
  CHECK_THROWS_AS(build_radial_states(g, std::vector<W>{0, 0}, radial, 3), std::invalid_argument);
  CHECK_THROWS_AS(build_radial_states(g, std::vector<W>{0}, radial, 3), std::invalid_argument);
  CHECK_THROWS_AS(build_radial_states(g, std::vector<W>{0, 1}, EdgeClassification{}, 3),
                  std::invalid_argument);
  CHECK_THROWS_AS(build_radial_states(g, std::vector<W>{0, 1}, radial, 0), std::invalid_argument);
}

TEST_CASE("local classification: corridor edges are macro, heavy intercell edges radial") {
  using W = std::uint64_t;
  // Pivots 0 and 1 (b = 2), cells {0, 2} and {1, 3}; r_2 = r_3 = 100, sigma(0, 1) = 1.
  // Edge 3 (2-3, w = 200) is intercell with 200 > N^3 sigma = 64: radial.
  const std::vector<Edge<W>> edges{{0, 1, 1}, {0, 2, 100}, {1, 3, 100}, {2, 3, 200}};
  const auto g = CsrGraph<W>::from_edges(4, edges);
  const auto od = sssp::paper::owner_decomposition(g, 2);
  REQUIRE(od.pivots == std::vector<Vertex>{0, 1});
  const auto cls = classify_local_edges(g, od);
  CHECK(cls.edge_class == std::vector<EdgeClass>{EdgeClass::kMacro, EdgeClass::kRadial,
                                                 EdgeClass::kRadial, EdgeClass::kRadial});
  CHECK(cls.macro_key[0] == od.corridor_max(0, 1));
  CHECK(cls.macro_key[3] == kNoCorridor);
  // With macro exponent 4 (N^4 sigma = 256 >= 200) the same edge would be macro.
  CHECK(classify_local_edges(g, od, 4).edge_class[3] == EdgeClass::kMacro);

  const auto rs = local_radial_states(g, od);
  CHECK(rs.shell_ratio_log2 == 5);  // h = 3, 8h = 24 -> K = 32
  CHECK(rs.coordinate == od.radius);
  CHECK(rs.graph.num_states == 2 + 2 + 3 + 3);
}

TEST_CASE("local classification omits edges outside the reduced graph") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(3, {{0, 1, 1}, {1, 2, 1}, {0, 2, 5}, {2, 2, 1}});
  const auto od = sssp::paper::owner_decomposition(g, 4);
  const auto cls = classify_local_edges(g, od);
  CHECK(cls.edge_class[2] == EdgeClass::kOmitted);
  CHECK(cls.edge_class[3] == EdgeClass::kOmitted);
  CHECK(cls.edge_class[0] == EdgeClass::kRadial);
}
