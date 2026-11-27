// Sec 7.3 (Lemmas 7.3-7.6) of arXiv 2611.34911v1 on a hand-computed auxiliary graph: tree and
// graph pieces from the overlap of [r/4, 4R], the dyadic grid of each graph piece (7 octaves
// below the origin, at least 2 above the largest R), widths 32q <= r < 64q, upper markers
// R < U+ <= 2R, the (C, q) width calendars merged into one marker chain per graph piece, the
// interpiece marker arcs U+_u -> 16 q_v, and the atoms of the edges w' < 2 max(q_u, q_v).
//
// Calendar values are stored after scaling every length by 2^7 (Lemma 7.3 l.1723-1724: "the
// lowest few boundaries can equivalently be represented after scaling all quantities by a fixed
// power of two"); the expectations below are written in those scaled units (real value x 128).
//
// The auxiliary graph is written by hand (it is a valid output of Sec 7.2: T' is an MST of G',
// rooted at a = 0, with r, R the root bottleneck and distance in T', depth <= h and R <= h r).

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
#include "paper/owner/calendars.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::owner::AuxGraph;
using sssp::paper::owner::build_calendars;
using sssp::paper::owner::calendar_scale;
using sssp::paper::owner::check_atoms;
using sssp::paper::owner::check_calendars;
using sssp::paper::owner::find_atoms;
using sssp::paper::owner::kAuxRoot;
using sssp::paper::owner::kCalendarScaleLog2;
using sssp::paper::owner::kNoAtom;
using sssp::paper::owner::kNoCalendar;
using sssp::paper::owner::kNoPiece;

namespace {

/// G' with root 0 and the given T' edge ids; T' is rooted by a BFS over the tree arcs in CSR
/// order, r / R are the root bottleneck / distance along it (Sec 7.3 l.1693).
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

// The hand example (b = 4, h = 5). Vertices 0 (= a), 1..8.
//   id: 0 (0,1) 64      1 (1,2) 2      2 (2,3) 3      3 (1,4) 320     4 (0,5) 100000
//       5 (5,6) 7       6 (6,7) 100    7 (1,3) 6      8 (3,5) 100001  9 (7,4) 150000
//       10 (2,4) 330    11 (4,8) 5
// T' = {0,1,2,3,4,5,6,11}; BFS order 0, 1, 5, 2, 4, 6, 3, 8, 7.
//   r: 1,2,3 -> 64; 4,8 -> 320; 5,6,7 -> 100000.  R: 1 64, 2 66, 3 69, 4 384, 8 389,
//   5 100000, 6 100007, 7 100107.
// Overlap r_v <= 16 R_u and r_u <= 16 R_v holds on every tree edge off a, on edges 7 and 10, and
// fails on edges 8 (r_5 = 100000 > 16 R_3 = 1104) and 9 (r_7 = 100000 > 16 R_4 = 6144).
// Tree pieces = graph pieces: {1,2,3,4,8} rooted at 1 and {5,6,7} rooted at 5.
// Grid of C0: origin 64, scaled grid[j] = 64 * 2^j, j = 0..12 (grid[12] = 262144 >= 4 * 128 * 389
// = 199168 > grid[11]). q: 1,2,3 -> grid[2] = 256 (real 2), 4,8 -> grid[4] = 1024 (real 8).
// U+: 1,2,3 -> grid[8] = 16384 (real 128); 4,8 -> grid[10] = 65536 (real 512).
// Grid of C1: origin 100000, grid[j] = 100000 * 2^j, j = 0..10 (grid[10] = 102400000 >=
// 4 * 128 * 100107 = 51254784 > grid[9]). q = grid[2] = 400000, U+ = grid[8] = 25600000.
// Calendars, coarse to fine: 0 = (C0, j 4) 64 cells of 1024; 1 = (C0, j 2) 64 cells of 256;
// 2 = (C1, j 2) 64 cells of 400000. Chain of C0: 65 + 65 - 17 shared = 113 markers (ids 0..112);
// chain of C1: 65 markers (ids 113..177).
// Interpiece arcs: edge 8 tail 3 head 5: U+_3 = 16384 (id 64) -> 16 q_5 = 6400000 (id 129);
// edge 9 tail 4 head 7: U+_4 = 65536 (id 112) -> 16 q_7 (id 129).
// Atoms (w' * 128 < 2 max(q)): edges 1, 2 (512), 11 (2048), 5, 6 (800000); edge 7 (768) is not.
// Atoms in representative order: {1,2,3} L 5; {5,6,7} L 107; {4,8} L 5.
template <class W>
AuxGraph<W> hand_example() {
  return make_aux<W>(9,
                     {{0, 1, W{64}},
                      {1, 2, W{2}},
                      {2, 3, W{3}},
                      {1, 4, W{320}},
                      {0, 5, W{100000}},
                      {5, 6, W{7}},
                      {6, 7, W{100}},
                      {1, 3, W{6}},
                      {3, 5, W{100001}},
                      {7, 4, W{150000}},
                      {2, 4, W{330}},
                      {4, 8, W{5}}},
                     {0, 1, 2, 3, 4, 5, 6, 11}, 5);
}

}  // namespace

TEMPLATE_TEST_CASE("calendar_scale multiplies by 2^7", "[unit][calendars]", std::uint64_t, double) {
  using W = TestType;
  STATIC_REQUIRE(kCalendarScaleLog2 == 7);
  CHECK(calendar_scale(W{0}) == W{0});
  CHECK(calendar_scale(W{3}) == W{384});
  if constexpr (std::is_floating_point_v<W>) CHECK(calendar_scale(W{0.25}) == W{32});
  CHECK_THROWS_AS(calendar_scale(std::numeric_limits<W>::max()), std::overflow_error);
}

TEMPLATE_TEST_CASE("Sec 7.3 pieces, grids, widths and upper markers on the hand example",
                   "[unit][calendars]", std::uint64_t, double) {
  using W = TestType;
  const auto aux = hand_example<W>();
  Diagnostics diag;
  const auto cal = build_calendars(aux, &diag);

  const std::vector<std::uint32_t> tp{kNoPiece, 0, 0, 0, 0, 1, 1, 1, 0};
  CHECK(cal.tree_piece == tp);
  CHECK(cal.graph_piece == tp);
  CHECK(cal.tree_piece_root == std::vector<Vertex>{1, 5});
  CHECK(cal.tree_piece_graph_piece == std::vector<std::uint32_t>{0, 1});
  const std::vector<std::uint8_t> overlap{0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1};
  CHECK(cal.overlap == overlap);

  REQUIRE(cal.graph_pieces.size() == 2);
  const auto& c0 = cal.graph_pieces[0];
  CHECK(c0.origin == W{64});
  CHECK(c0.grid_root == 1);
  REQUIRE(c0.grid.size() == 13);
  for (std::size_t j = 0; j < 13; ++j) CHECK(c0.grid[j] == W{64} * static_cast<W>(1ULL << j));
  const auto& c1 = cal.graph_pieces[1];
  CHECK(c1.origin == W{100000});
  CHECK(c1.grid_root == 5);
  REQUIRE(c1.grid.size() == 11);
  for (std::size_t j = 0; j < 11; ++j) {
    CHECK(c1.grid[j] == W{100000} * static_cast<W>(1ULL << j));
  }
  CHECK(cal.num_grid_points() == 24);

  const std::uint32_t nr = kNoPiece;
  CHECK(cal.radius_index == std::vector<std::uint32_t>{nr, 7, 7, 7, 9, 7, 7, 7, 9});
  CHECK(cal.width_index == std::vector<std::uint32_t>{nr, 2, 2, 2, 4, 2, 2, 2, 4});
  CHECK(cal.upper_index == std::vector<std::uint32_t>{nr, 8, 8, 8, 10, 8, 8, 8, 10});
  CHECK(cal.width == std::vector<W>{W{0}, W{256}, W{256}, W{256}, W{1024}, W{400000}, W{400000},
                                    W{400000}, W{1024}});
  CHECK(cal.upper == std::vector<W>{W{0}, W{16384}, W{16384}, W{16384}, W{65536}, W{25600000},
                                    W{25600000}, W{25600000}, W{65536}});
  CHECK(cal.calendar == std::vector<std::uint32_t>{kNoCalendar, 1, 1, 1, 0, 2, 2, 2, 0});
  // Every vertex: 32 q <= 2^7 r < 64 q and 2^7 R < U+ <= 2^8 R.
  for (Vertex v = 1; v < aux.num_vertices(); ++v) {
    const W sr = calendar_scale(aux.radius[v]);
    const W sR = calendar_scale(aux.owner_distance[v]);
    CHECK(W{32} * cal.width[v] <= sr);
    CHECK(sr < W{64} * cal.width[v]);
    CHECK(sR < cal.upper[v]);
    CHECK(cal.upper[v] <= W{2} * sR);
  }

#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["calendars.tree_pieces"] == 2);
  CHECK(diag.counters["calendars.graph_pieces"] == 2);
  CHECK(diag.counters["calendars.grid_points"] == 24);
  CHECK(diag.counters["calendars.overlap_edges"] == 8);
  CHECK(diag.counters["calendars.grid_search_steps"] > 0);
#endif
}

TEMPLATE_TEST_CASE("Sec 7.3 / Lemma 7.4: width calendars and merged marker chains",
                   "[unit][calendars]", std::uint64_t, double) {
  using W = TestType;
  const auto aux = hand_example<W>();
  Diagnostics diag;
  const auto cal = build_calendars(aux, &diag);

  REQUIRE(cal.calendars.size() == 3);
  const W widths[] = {W{1024}, W{256}, W{400000}};
  const std::uint32_t pieces[] = {0, 0, 1};
  const std::uint32_t grid_index[] = {4, 2, 2};
  for (std::size_t c = 0; c < 3; ++c) {
    INFO("calendar " << c);
    const auto& wc = cal.calendars[c];
    CHECK(wc.graph_piece == pieces[c]);
    CHECK(wc.grid_index == grid_index[c]);
    CHECK(wc.width == widths[c]);
    REQUIRE(wc.num_cells() == 64);
    REQUIRE(wc.markers.size() == 65);
    for (std::size_t i = 0; i <= 64; ++i) {
      CHECK(cal.marker_value[wc.markers[i]] == static_cast<W>(i) * widths[c]);
    }
  }
  CHECK(cal.graph_pieces[0].calendars == std::vector<std::uint32_t>{0, 1});
  CHECK(cal.graph_pieces[1].calendars == std::vector<std::uint32_t>{2});
  CHECK(cal.num_cells() == 192);

  // Chains: ids [first, end) per graph piece, strictly increasing values.
  REQUIRE(cal.marker_value.size() == 178);
  CHECK(cal.graph_pieces[0].first_marker == 0);
  CHECK(cal.graph_pieces[0].end_marker == 113);
  CHECK(cal.graph_pieces[1].first_marker == 113);
  CHECK(cal.graph_pieces[1].end_marker == 178);
  for (std::uint32_t m = 0; m < 178; ++m) CHECK(cal.marker_graph_piece[m] == (m < 113 ? 0U : 1U));
  // C0: 0, 256, ..., 16384 (ids 0..64), then 17408, ..., 65536 (ids 65..112).
  for (std::uint32_t m = 0; m <= 64; ++m) CHECK(cal.marker_value[m] == W{256} * static_cast<W>(m));
  for (std::uint32_t m = 65; m < 113; ++m) {
    CHECK(cal.marker_value[m] == W{16384} + W{1024} * static_cast<W>(m - 64));
  }
  // The coarse calendar shares the fine markers below 16384: 1024 i is id 4 i.
  for (std::uint32_t i = 0; i <= 16; ++i) CHECK(cal.calendars[0].markers[i] == 4 * i);
  for (std::uint32_t i = 17; i <= 64; ++i) CHECK(cal.calendars[0].markers[i] == 48 + i);
  for (std::uint32_t i = 0; i <= 64; ++i) CHECK(cal.calendars[1].markers[i] == i);
  for (std::uint32_t i = 0; i <= 64; ++i) CHECK(cal.calendars[2].markers[i] == 113 + i);

#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["calendars.calendars"] == 3);
  CHECK(diag.counters["calendars.cells"] == 192);
  CHECK(diag.counters["calendars.markers"] == 178);
#endif
}

TEMPLATE_TEST_CASE("Lemma 7.6: interpiece marker arcs and the topological marker order",
                   "[unit][calendars]", std::uint64_t, double) {
  using W = TestType;
  const auto aux = hand_example<W>();
  Diagnostics diag;
  const auto cal = build_calendars(aux, &diag);

  REQUIRE(cal.interpiece_arcs.size() == 2);
  const auto& a8 = cal.interpiece_arcs[0];
  CHECK(a8.aux_edge == 8);
  CHECK(a8.tail == 3);
  CHECK(a8.head == 5);
  CHECK(a8.from == 64);
  CHECK(a8.to == 129);
  const auto& a9 = cal.interpiece_arcs[1];
  CHECK(a9.aux_edge == 9);
  CHECK(a9.tail == 4);
  CHECK(a9.head == 7);
  CHECK(a9.from == 112);
  CHECK(a9.to == 129);
  CHECK(cal.marker_value[129] == W{6400000});
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["calendars.interpiece_arcs"] == 2);
#endif

  // marker_order is a permutation in which every chain arc and every interpiece arc points
  // forward.
  REQUIRE(cal.marker_order.size() == 178);
  std::vector<std::size_t> pos(178, 178);
  for (std::size_t k = 0; k < 178; ++k) {
    REQUIRE(cal.marker_order[k] < 178);
    REQUIRE(pos[cal.marker_order[k]] == 178);
    pos[cal.marker_order[k]] = k;
  }
  for (std::uint32_t m = 0; m + 1 < 178; ++m) {
    if (m + 1 != 113) CHECK(pos[m] < pos[m + 1]);
  }
  for (const auto& a : cal.interpiece_arcs) {
    CHECK(pos[a.from] < pos[a.to]);
    CHECK(cal.marker_value[a.from] < cal.marker_value[a.to]);
  }
  CHECK(check_calendars(aux, cal).empty());
}

TEMPLATE_TEST_CASE("Lemma 7.5: atoms of the light auxiliary edges", "[unit][atoms]", std::uint64_t,
                   double) {
  using W = TestType;
  const auto aux = hand_example<W>();
  const auto cal = build_calendars(aux);
  Diagnostics diag;
  const auto atoms = find_atoms(aux, cal, &diag);

  REQUIRE(atoms.num_atoms() == 3);
  CHECK(atoms.atom == std::vector<std::uint32_t>{kNoAtom, 0, 0, 0, 2, 1, 1, 1, 2});
  CHECK(atoms.offsets == std::vector<std::uint32_t>{0, 3, 6, 8});
  CHECK(atoms.vertices == std::vector<Vertex>{1, 2, 3, 5, 6, 7, 4, 8});
  CHECK(atoms.representative == std::vector<Vertex>{1, 5, 4});
  CHECK(atoms.tree_offsets == std::vector<std::uint32_t>{0, 2, 4, 5});
  CHECK(atoms.tree_edges == std::vector<EdgeId>{1, 2, 5, 6, 11});
  CHECK(atoms.length == std::vector<W>{W{5}, W{107}, W{5}});
  CHECK(atoms.calendar == std::vector<std::uint32_t>{1, 2, 0});
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["atoms.atoms"] == 3);
  CHECK(diag.counters["atoms.defining_edges"] == 5);
  CHECK(diag.counters["atoms.max_size"] == 3);
#endif
  CHECK(check_atoms(aux, cal, atoms).empty());
}

TEMPLATE_TEST_CASE("a defining chord that closes a cycle stays out of the discovery tree",
                   "[unit][atoms]", std::uint64_t, double) {
  using W = TestType;
  // Root 0 -- 1 (64), path 1-2-3 (2, 3) and the chord 1-3 of weight 3 (2^7 * 3 = 384 < 512):
  // all three edges are defining, the tree uses the first two (edge order), L = 5.
  const auto aux =
      make_aux<W>(4, {{0, 1, W{64}}, {1, 2, W{2}}, {2, 3, W{3}}, {1, 3, W{3}}}, {0, 1, 2}, 5);
  const auto cal = build_calendars(aux);
  Diagnostics diag;
  const auto atoms = find_atoms(aux, cal, &diag);
  REQUIRE(atoms.num_atoms() == 1);
  CHECK(atoms.tree_edges == std::vector<EdgeId>{1, 2});
  CHECK(atoms.length == std::vector<W>{W{5}});
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["atoms.defining_edges"] == 3);
#endif
  CHECK(check_atoms(aux, cal, atoms).empty());
}

TEMPLATE_TEST_CASE("degenerate auxiliary graphs", "[unit][calendars]", std::uint64_t, double) {
  using W = TestType;
  SECTION("only the root") {
    const auto aux = make_aux<W>(1, {}, {}, 2);
    const auto cal = build_calendars(aux);
    CHECK(cal.graph_pieces.empty());
    CHECK(cal.calendars.empty());
    CHECK(cal.marker_value.empty());
    CHECK(cal.marker_order.empty());
    CHECK(check_calendars(aux, cal).empty());
    const auto atoms = find_atoms(aux, cal);
    CHECK(atoms.num_atoms() == 0);
    CHECK(check_atoms(aux, cal, atoms).empty());
  }
  SECTION("one child of the root") {
    // r = R = 3: grid[j] = 3 * 2^j, j = 0..9 (grid[9] = 1536 = 4 * 128 * 3);
    // q = grid[2] = 12 (real 3/32), U+ = grid[8] = 768 (real 6 = 2R): 64 cells.
    // h = 3 (b = 2): one nonroot vertex needs b >= 2 (components of T - A have < b vertices).
    const auto aux = make_aux<W>(2, {{1, 0, W{3}}}, {0}, 3);
    const auto cal = build_calendars(aux);
    REQUIRE(cal.graph_pieces.size() == 1);
    CHECK(cal.graph_pieces[0].grid.size() == 10);
    CHECK(cal.width[1] == W{12});
    CHECK(cal.upper[1] == W{768});
    REQUIRE(cal.calendars.size() == 1);
    CHECK(cal.calendars[0].num_cells() == 64);
    CHECK(cal.marker_value.size() == 65);
    CHECK(cal.interpiece_arcs.empty());
    CHECK(check_calendars(aux, cal).empty());
    const auto atoms = find_atoms(aux, cal);
    CHECK(atoms.num_atoms() == 1);
    CHECK(atoms.length == std::vector<W>{W{0}});
    CHECK(check_atoms(aux, cal, atoms).empty());
  }
}

TEMPLATE_TEST_CASE("build_calendars rejects malformed input", "[unit][calendars]", std::uint64_t,
                   double) {
  using W = TestType;
  SECTION("sizes do not match") {
    auto aux = hand_example<W>();
    aux.radius.pop_back();
    CHECK_THROWS_AS(build_calendars(aux), std::invalid_argument);
  }
  SECTION("a zero radius off the root (the input is positive, Sec 7.1 l.1642)") {
    auto aux = make_aux<W>(3, {{0, 1, W{0}}, {1, 2, W{0}}}, {0, 1}, 3);
    CHECK_THROWS_AS(build_calendars(aux), std::invalid_argument);
  }
  SECTION("2^10 R does not fit the type") {
    if constexpr (std::is_integral_v<W>) {
      const auto aux = make_aux<W>(2, {{0, 1, W{1} << 55}}, {0}, 2);
      CHECK_THROWS_AS(build_calendars(aux), std::overflow_error);
    } else {
      const auto aux = make_aux<W>(2, {{0, 1, std::numeric_limits<W>::max() / 4}}, {0}, 2);
      CHECK_THROWS_AS(build_calendars(aux), std::overflow_error);
    }
  }
  SECTION("find_atoms with calendars of another graph") {
    const auto aux = hand_example<W>();
    const auto other = make_aux<W>(2, {{1, 0, W{3}}}, {0}, 2);
    const auto cal = build_calendars(other);
    CHECK_THROWS_AS(find_atoms(aux, cal), std::invalid_argument);
  }
}

TEMPLATE_TEST_CASE("check_calendars and check_atoms report corrupted outputs", "[unit][calendars]",
                   std::uint64_t, double) {
  using W = TestType;
  const auto aux = hand_example<W>();
  const auto good = build_calendars(aux);
  const auto good_atoms = find_atoms(aux, good);
  REQUIRE(check_calendars(aux, good).empty());
  REQUIRE(check_atoms(aux, good, good_atoms).empty());

  SECTION("width off by one octave") {
    auto cal = good;
    cal.width[2] = cal.width[2] + cal.width[2];
    cal.width_index[2] += 1;
    CHECK_FALSE(check_calendars(aux, cal).empty());
  }
  SECTION("upper marker one octave too high") {
    auto cal = good;
    cal.upper[3] = cal.upper[3] + cal.upper[3];
    cal.upper_index[3] += 1;
    CHECK_FALSE(check_calendars(aux, cal).empty());
  }
  SECTION("grid shifted to 6 octaves below the origin") {
    auto cal = good;
    cal.graph_pieces[0].grid.erase(cal.graph_pieces[0].grid.begin());
    CHECK_FALSE(check_calendars(aux, cal).empty());
  }
  SECTION("two tree pieces merged") {
    auto cal = good;
    cal.tree_piece[5] = cal.tree_piece[6] = cal.tree_piece[7] = 0;
    CHECK_FALSE(check_calendars(aux, cal).empty());
  }
  SECTION("marker chain not increasing") {
    auto cal = good;
    std::swap(cal.marker_value[10], cal.marker_value[11]);
    CHECK_FALSE(check_calendars(aux, cal).empty());
  }
  SECTION("missing interpiece arc") {
    auto cal = good;
    cal.interpiece_arcs.pop_back();
    CHECK_FALSE(check_calendars(aux, cal).empty());
  }
  SECTION("interpiece arc against its orientation") {
    auto cal = good;
    std::swap(cal.interpiece_arcs[0].from, cal.interpiece_arcs[0].to);
    std::swap(cal.interpiece_arcs[0].tail, cal.interpiece_arcs[0].head);
    CHECK_FALSE(check_calendars(aux, cal).empty());
  }
  SECTION("marker order not topological") {
    auto cal = good;
    std::swap(cal.marker_order.front(), cal.marker_order.back());
    CHECK_FALSE(check_calendars(aux, cal).empty());
  }
  SECTION("calendar with a cell missing at the end") {
    auto cal = good;
    cal.calendars[1].markers.pop_back();
    CHECK_FALSE(check_calendars(aux, cal).empty());
  }
  SECTION("atom split") {
    auto atoms = good_atoms;
    atoms.atom[3] = 2;
    CHECK_FALSE(check_atoms(aux, good, atoms).empty());
  }
  SECTION("atom length wrong") {
    auto atoms = good_atoms;
    atoms.length[0] = W{4};
    CHECK_FALSE(check_atoms(aux, good, atoms).empty());
  }
  SECTION("atom larger than b - 1 = 3") {
    auto big = aux;
    big.hop_bound = 3;  // b = 2: the atoms {1,2,3}, {5,6,7}, {4,8} are too large
    CHECK_FALSE(check_atoms(big, good, good_atoms).empty());
  }
  SECTION("representative not the top of T'[Y]") {
    auto atoms = good_atoms;
    atoms.representative[0] = 2;
    CHECK_FALSE(check_atoms(aux, good, atoms).empty());
  }
}

TEMPLATE_TEST_CASE("paper outer separator: one pivot, calendars of the whole graph",
                   "[unit][calendars]", std::uint64_t, double) {
  using W = TestType;
  // A path with a heavy middle edge: with b = L_N^64 there is one pivot (vertex 0), G' is the
  // whole graph rooted at 0 and the heavy edge splits two tree pieces.
  const std::vector<Edge<W>> edges{
      {0, 1, W{4}}, {1, 2, W{1}}, {2, 3, W{5000}}, {3, 4, W{2}}, {0, 4, W{9000}}};
  const auto g = CsrGraph<W>::from_edges(5, edges);
  const auto od = sssp::paper::owner_decomposition(
      g, sssp::paper::owner::outer_separator_size(g.num_vertices()));
  REQUIRE(od.pivots.size() == 1);
  REQUIRE(od.pivots[0] == 0);  // Q-31: the unmarked case uses vertex 0
  sssp::paper::model::WordRng rng(7);
  const auto cls = sssp::paper::owner::classify_interowner_edges(g, od, rng);
  const auto aux = sssp::paper::owner::build_aux_graph(g, od, cls);
  Diagnostics diag;
  const auto cal = build_calendars(aux, &diag);
  CHECK(check_calendars(aux, cal).empty());
  CHECK(cal.tree_piece_root.size() == 2);
  CHECK(diag.fallbacks.empty());
  const auto atoms = find_atoms(aux, cal, &diag);
  CHECK(check_atoms(aux, cal, atoms).empty());
}
