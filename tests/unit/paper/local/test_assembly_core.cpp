// Unit tests of the shared assembly core (task P5-LOCAL-UNIFY, Q-59; arXiv 2611.34911v1
// Sec 5.4 l.1182-1243 and the proof of Lemma 5.10 l.1265-1300):
//   * LengthScale applies its integer factors and doublings in order and saturates (Q-38);
//   * Rule::paper(N) carries the paper's constants: c = 10 (l.1206) and N^16 (l.1216);
//   * the support grid entry point (key_supports + build_support_grid), which Sec 8.2 will call
//     with its own unit and number of octaves: hand examples (separate components, overlap,
//     a class stopped at the next component l.1279, lowered start, cap), the same grid as the
//     local assembly when given the local unit and N^16, and other units / octave counts;
//   * doubling_points (separate grids, l.1209-1213) and aligned_by_medians (l.1214-1218);
//   * ClassMerger (simultaneous merges, Lemma 3.6; l.1281) and link_and_order (l.1282-1293) on
//     a hand example with an auxiliary join vertex;
//   * assemble() with the paper's rule equals paper::assemble_local_hierarchy on every field,
//     and other injected constants change the alignment exactly as their test says.
// Expected values are computed by hand from the cited lines, not copied from a run.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/assemble.hpp"
#include "paper/local/assembly_core.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/saturating.hpp"
#include "paper/local/threshold_joins.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/prng.hpp"
#include "unit/paper/local/assembly_hash.hpp"
#include "unit/paper/local/clustered_graphs.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::WeightTraits;
using sssp::paper::BlockScaleKind;
using sssp::paper::GridSupport;
using sssp::paper::HierarchyNode;
using sssp::paper::kNoNode;
using sssp::paper::NodeId;
using sssp::paper::NodeOrigin;
using sssp::paper::StateArc;
namespace core = sssp::paper::assembly;
namespace sat = sssp::paper::sat;
namespace lt = sssp::test::local;

namespace {

template <class W>
struct Instance {
  sssp::paper::OwnerDecomposition<W> od;
  sssp::paper::RadialStateGraph<W> rs;
  sssp::paper::ThresholdJoins<W> joins;
};

template <class W>
Instance<W> build(const CsrGraph<W>& g, std::size_t b) {
  auto od = sssp::paper::owner_decomposition(g, b);
  auto rs = sssp::paper::local_radial_states(g, od);
  auto joins = sssp::paper::compute_joins(od, rs);
  return {std::move(od), std::move(rs), std::move(joins)};
}

template <class W>
std::uint64_t hash_of(const sssp::paper::LocalAssembly<W>& a) {
  lt::Fnv1a f;
  lt::mix_assembly(f, a);
  return f.digest();
}

template <class W>
std::vector<W> lowers(const std::vector<GridSupport<W>>& s) {
  std::vector<W> out;
  for (const auto& x : s) out.push_back(x.lower);
  return out;
}

template <class W>
std::vector<W> uppers(const std::vector<GridSupport<W>>& s) {
  std::vector<W> out;
  for (const auto& x : s) out.push_back(x.upper);
  return out;
}

template <class W>
std::vector<std::pair<std::size_t, std::size_t>> indices(const std::vector<GridSupport<W>>& s) {
  std::vector<std::pair<std::size_t, std::size_t>> out;
  for (const auto& x : s) out.emplace_back(x.lower_index, x.upper_index);
  return out;
}

using Pairs = std::vector<std::pair<std::size_t, std::size_t>>;

}  // namespace

TEMPLATE_TEST_CASE("LengthScale: integer factors and doublings in order, saturating (Q-38)", "",
                   std::uint64_t, double) {
  using W = TestType;
  const W inf = WeightTraits<W>::infinity();
  CHECK(core::LengthScale{}.apply(W{7}) == W{7});
  CHECK(core::LengthScale::int_power(3, 2).apply(W{5}) == W{45});
  CHECK(core::LengthScale::octaves(3).apply(W{5}) == W{40});
  CHECK(core::LengthScale{}.times(2).times(7).times(3).doubled(4).apply(W{1}) == W{672});
  CHECK(core::LengthScale::int_power(6, 0).apply(W{9}) == W{9});
  CHECK(core::LengthScale::saturated().apply(W{1}) == inf);
  CHECK(core::LengthScale::saturated().times(3).apply(W{1}) == inf);
  // The same operations as paper/local/saturating.hpp, hence the same saturation.
  for (const std::uint64_t n : {2U, 32U, 41U, 1000U}) {
    for (const unsigned e : {1U, 3U, 10U, 16U}) {
      for (const W x : {W{1}, W{3}, W{1000}, static_cast<W>(std::uint64_t{1} << 40)}) {
        CHECK(core::LengthScale::int_power(n, e).apply(x) == sat::mul_int_pow(x, n, e));
        CHECK(core::LengthScale::octaves(e).apply(x) == sat::mul_pow2(x, e));
      }
    }
  }
  if constexpr (std::is_same_v<W, std::uint64_t>) {
    CHECK(core::LengthScale::int_power(32, 16).apply(W{1}) == inf);  // 2^80 does not fit
    CHECK(core::LengthScale::octaves(64).apply(W{1}) == inf);
  }
}

TEMPLATE_TEST_CASE("Rule::paper carries c = 10 and N^16 (Sec 5.4 l.1206, l.1216)", "",
                   std::uint64_t, double) {
  using W = TestType;
  for (const std::size_t n : {32U, 41U, 100U}) {
    const auto rule = core::Rule::paper(n);
    for (const W x : {W{1}, W{5}, W{123456}}) {
      CHECK(rule.align.apply(x) == sat::mul_int_pow(x, n, 10));
      CHECK(rule.support_end.apply(x) == sat::mul_int_pow(x, n, 16));
    }
  }
  CHECK(sssp::paper::kLocalAlignExponent == 10);
  CHECK(sssp::paper::kLocalSupportExponent == 16);
}

TEST_CASE("Support grid: hand examples (Sec 5.4 l.1214-1224, l.1279)") {
  using W = std::uint64_t;
  const auto grid_of = [](std::vector<W> units, unsigned octaves, W cap) {
    return core::build_support_grid(
        core::key_supports<W>(units, core::LengthScale::octaves(octaves)), cap);
  };
  SECTION("separate components") {
    // Supports [8, 32] and [1000, 4000]: two doubling grids.
    const auto g = grid_of({8, 1000}, 2, 5000);
    CHECK(g.grid == std::vector<W>{8, 16, 32, 1000, 2000, 4000});
    CHECK(g.grid_upper == std::vector<W>{16, 32, 64, 2000, 4000, 8000});
    CHECK(lowers(g.supports) == std::vector<W>{8, 1000});
    CHECK(uppers(g.supports) == std::vector<W>{32, 4000});
    CHECK(indices(g.supports) == Pairs{{0, 3}, {3, 6}});
    CHECK(g.classify(8, 0) == 0);
    CHECK(g.classify(20, 0) == 1);
    CHECK(g.classify(63, 0) == 2);
    CHECK_FALSE(g.classify(64, 0).has_value());    // above the last class of key 0
    CHECK_FALSE(g.classify(7, 0).has_value());     // below the grid of key 0
    CHECK_FALSE(g.classify(1500, 0).has_value());  // a class of another key
    CHECK(g.classify(3000, 1) == 4);
  }
  SECTION("overlapping supports share one grid") {
    // [8, 32] and [20, 80] overlap: one component [8, 80].
    const auto g = grid_of({8, 20}, 2, 5000);
    CHECK(g.grid == std::vector<W>{8, 16, 32, 64});
    CHECK(g.grid_upper == std::vector<W>{16, 32, 64, 128});
    CHECK(indices(g.supports) == Pairs{{0, 3}, {1, 4}});
    CHECK(g.classify(70, 1) == 3);
  }
  SECTION("keys in any order; touching supports merge") {
    // [16, 64] and [4, 16] touch at 16: one component [4, 64].
    const auto g = grid_of({16, 4}, 2, 5000);
    CHECK(g.grid == std::vector<W>{4, 8, 16, 32, 64});
    CHECK(indices(g.supports) == Pairs{{2, 5}, {0, 3}});
  }
  SECTION("the last class of a component stops at the next component (l.1279)") {
    // [8, 32] and [40, 160]: the class starting at 32 ends at 40, not at 64.
    const auto g = grid_of({8, 40}, 2, 5000);
    CHECK(g.grid == std::vector<W>{8, 16, 32, 40, 80, 160});
    CHECK(g.grid_upper == std::vector<W>{16, 32, 40, 80, 160, 320});
    CHECK(g.classify(39, 0) == 2);
    CHECK_FALSE(g.classify(40, 0).has_value());
    CHECK(g.classify(40, 1) == 3);
  }
  SECTION("no grid point above the cap") {
    const auto g = grid_of({8, 40}, 2, 20);
    CHECK(g.grid == std::vector<W>{8, 16, 40});
    CHECK(g.grid_upper == std::vector<W>{16, 32, 80});
    CHECK(indices(g.supports) == Pairs{{0, 2}, {2, 3}});
  }
  SECTION("a lowered start moves the grid origin") {
    auto supports = core::key_supports<W>(std::vector<W>{8}, core::LengthScale::octaves(2));
    supports[0].lower = 3;  // an aligned block with s_U = 3 (l.1219)
    const auto g = core::build_support_grid(std::move(supports), W{5000});
    CHECK(g.grid == std::vector<W>{3, 6, 12, 24});
    CHECK(g.grid_upper == std::vector<W>{6, 12, 24, 48});
    CHECK(indices(g.supports) == Pairs{{0, 4}});
  }
  SECTION("no keys") {
    const auto g = grid_of({}, 2, 5000);
    CHECK(g.grid.empty());
    CHECK(g.grid_upper.empty());
    CHECK(g.supports.empty());
  }
  SECTION("a saturated support with an infinite cap terminates") {
    const W inf = WeightTraits<W>::infinity();
    const W big = W{1} << 62;
    const auto g = grid_of({big}, 5, inf);
    CHECK(g.grid == std::vector<W>{big, big + big});
    CHECK(g.grid_upper == std::vector<W>{big + big, inf});
    CHECK(uppers(g.supports) == std::vector<W>{inf});
  }
}

TEST_CASE("Support grid: the unit and the number of octaves are the caller's (Sec 8.2 entry)") {
  using W = std::uint64_t;
  const std::vector<W> units{8, 40, 3000};
  const auto base = core::build_support_grid(
      core::key_supports<W>(units, core::LengthScale::octaves(2)), W{20000});
  SECTION("scaling every unit and the cap scales the grid") {
    // Doubling grids are homogeneous: unit * 4 gives every grid point * 4, same indices.
    std::vector<W> scaled;
    for (const W u : units) scaled.push_back(4 * u);
    const auto g = core::build_support_grid(
        core::key_supports<W>(scaled, core::LengthScale::octaves(2)), W{80000});
    REQUIRE(g.grid.size() == base.grid.size());
    for (std::size_t i = 0; i < g.grid.size(); ++i) {
      CHECK(g.grid[i] == 4 * base.grid[i]);
      CHECK(g.grid_upper[i] == 4 * base.grid_upper[i]);
    }
    CHECK(indices(g.supports) == indices(base.supports));
  }
  SECTION("more octaves extend each support upwards") {
    // 2 octaves: [8, 32], [40, 160], [3000, 12000]. 5 octaves: [8, 256] and [40, 1280] merge.
    CHECK(base.grid == std::vector<W>{8, 16, 32, 40, 80, 160, 3000, 6000, 12000});
    const auto g = core::build_support_grid(
        core::key_supports<W>(units, core::LengthScale::octaves(5)), W{200000});
    CHECK(uppers(g.supports) == std::vector<W>{256, 1280, 96000});
    CHECK(g.grid == std::vector<W>{8, 16, 32, 64, 128, 256, 512, 1024, 3000, 6000, 12000, 24000,
                                   48000, 96000});
    CHECK(g.grid_upper.back() == 192000);
    CHECK(g.grid_upper[7] == 2048);  // 2048 < 3000: not cut by the next component
    CHECK(indices(g.supports) == Pairs{{0, 6}, {2, 8}, {8, 14}});
  }
  SECTION("an integer-power end, as in Sec 5.4") {
    const auto s = core::key_supports<W>(units, core::LengthScale::int_power(3, 2));
    CHECK(lowers(s) == units);
    CHECK(uppers(s) == std::vector<W>{72, 360, 27000});
  }
}

TEMPLATE_TEST_CASE("Support grid entry reproduces the grid of the local assembly", "",
                   std::uint64_t, double) {
  using W = TestType;
  std::size_t checked = 0;
  std::size_t lowered = 0;
  const auto check = [&](const CsrGraph<W>& g, std::size_t b) {
    const auto in = build(g, b);
    const std::size_t n = g.num_vertices();
    const auto a = sssp::paper::assemble_local_hierarchy(in.od, in.rs, in.joins);
    // The local unit is the corridor weight (named key, Q-37), the end N^16 sigma (l.1216).
    std::vector<W> units;
    for (const auto& c : in.od.corridors) units.push_back(c.weight);
    auto supports = core::key_supports<W>(units, core::LengthScale::int_power(n, 16));
    for (const auto& bs : a.blocks) {  // l.1219: an aligned block lowers the start to s_U
      if (bs.kind != BlockScaleKind::kAligned) continue;
      if (bs.s < supports[bs.sigma_key].lower) ++lowered;
      supports[bs.sigma_key].lower = std::min(supports[bs.sigma_key].lower, bs.s);
    }
    W cap{0};  // the largest join weight
    for (const auto& j : in.joins.joins) cap = std::max(cap, j.weight);
    const auto sg = core::build_support_grid(std::move(supports), cap);
    CHECK(sg.grid == a.grid);
    CHECK(sg.grid_upper == a.grid_upper);
    CHECK(lowers(sg.supports) == lowers(a.supports));
    CHECK(uppers(sg.supports) == uppers(a.supports));
    CHECK(indices(sg.supports) == indices(a.supports));
    // Every node of a global class has its height on the grid.
    const auto& h = a.result.hierarchy;
    for (std::size_t v = 0; v < h.nodes.size(); ++v) {
      if (h.is_leaf(static_cast<NodeId>(v)) || a.node_origin[v] != NodeOrigin::kGlobalGrid) {
        continue;
      }
      CHECK(std::binary_search(sg.grid.begin(), sg.grid.end(), h.nodes[v].height));
    }
    ++checked;
  };
  for (const auto f :
       {lt::Family::kPath, lt::Family::kRandomTree, lt::Family::kSparse, lt::Family::kGrid}) {
    for (const auto wd : {lt::Weights::kTies, lt::Weights::kUniform, lt::Weights::kExponential}) {
      for (const std::size_t b : {2U, 4U}) {
        INFO(lt::name(f) << " " << lt::name(wd) << " b=" << b);
        check(CsrGraph<W>::from_edges(36, lt::make_graph<W>(f, wd, 36, 77 + b)), b);
      }
    }
  }
  if constexpr (std::is_same_v<W, std::uint64_t>) {
    for (const std::size_t b : {2U, 3U, 5U}) {
      INFO("clustered b=" << b);
      check(CsrGraph<W>::from_edges(40, lt::clustered(40, 40 + b)), b);
    }
  }
  CHECK(checked > 0);
  CHECK(lowered > 0);
}

TEST_CASE("Separate grid points double from the unit up to the largest value (l.1209-1213)") {
  using W = std::uint64_t;
  CHECK(core::doubling_points<W>(5, 0) == std::vector<W>{5});
  CHECK(core::doubling_points<W>(5, 9) == std::vector<W>{5});
  CHECK(core::doubling_points<W>(5, 10) == std::vector<W>{5, 10});
  CHECK(core::doubling_points<W>(5, 39) == std::vector<W>{5, 10, 20});
  CHECK(core::doubling_points<W>(5, 40) == std::vector<W>{5, 10, 20, 40});
  // Another unit: the points of the unit 7.
  CHECK(core::doubling_points<W>(7, 40) == std::vector<W>{7, 14, 28});
  CHECK(core::doubling_points<double>(0.25, 1.0) == std::vector<double>{0.25, 0.5, 1.0});
  // floor_index: the class of a value.
  const std::vector<W> pts{5, 10, 20};
  CHECK(core::floor_index<W>(pts, 0, 3, 4) == 3);  // below every point: hi
  CHECK(core::floor_index<W>(pts, 0, 3, 5) == 0);
  CHECK(core::floor_index<W>(pts, 0, 3, 19) == 1);
  CHECK(core::floor_index<W>(pts, 0, 3, 1000) == 2);
  CHECK(core::floor_index<W>(pts, 1, 3, 7) == 3);
}

TEST_CASE("Alignment by successive medians equals the per-item test (l.1206, l.1214-1218)") {
  using W = std::uint64_t;
  const auto scale = core::LengthScale::octaves(4);  // aligned iff 16 s >= sigma
  using Flags = std::vector<std::uint8_t>;
  CHECK(core::aligned_by_medians<W>(std::vector<W>{1, 2, 4, 8, 16}, 100, scale) ==
        Flags{0, 0, 0, 1, 1});
  CHECK(core::aligned_by_medians<W>(std::vector<W>{16, 1, 8, 2, 4}, 100, scale) ==
        Flags{1, 0, 1, 0, 0});
  CHECK(core::aligned_by_medians<W>(std::vector<W>{7, 7, 6, 6, 7}, 100, scale) ==
        Flags{1, 1, 0, 0, 1});  // 112 >= 100 > 96
  CHECK(core::aligned_by_medians<W>(std::vector<W>{1, 2, 3}, 100, scale) == Flags{0, 0, 0});
  CHECK(core::aligned_by_medians<W>(std::vector<W>{9, 8, 7}, 100, scale) == Flags{1, 1, 1});
  CHECK(core::aligned_by_medians<W>(std::vector<W>{}, 100, scale).empty());
  // Against the definition on pseudo-random keys.
  sssp::test::SplitMix64 rng(11);
  for (int round = 0; round < 50; ++round) {
    std::vector<W> s(1 + rng.below(40));
    for (auto& x : s) x = rng.between(1, 64);
    const W sigma = rng.between(1, 1500);
    const auto got = core::aligned_by_medians<W>(s, sigma, scale);
    REQUIRE(got.size() == s.size());
    for (std::size_t i = 0; i < s.size(); ++i) CHECK((got[i] != 0) == (16 * s[i] >= sigma));
  }
}

TEST_CASE("ClassMerger and link_and_order: hand example with an auxiliary vertex") {
  using W = std::uint64_t;
  // States 0, 1, 2 (leaves 0, 1, 2) and one auxiliary join vertex 3 without a node.
  std::vector<HierarchyNode<W>> nodes(3);
  for (std::uint32_t x = 0; x < 3; ++x) nodes[x].state = x;
  std::vector<W> upper(3, 0);
  std::vector<std::uint8_t> inclusive(3, 0);
  const std::vector<std::pair<std::uint32_t, std::uint32_t>> joins{{0, 1}, {1, 3}, {3, 2}, {0, 3}};
  core::ClassMerger merger(
      std::vector<NodeId>{0, 1, 2, kNoNode}, [&](std::uint32_t j) { return joins[j]; },
      [&](W height, std::vector<NodeId> children, W up, bool incl, NodeOrigin) {
        const auto id = static_cast<NodeId>(nodes.size());
        for (const NodeId c : children) nodes[c].parent = id;
        HierarchyNode<W> node;
        node.children = std::move(children);
        node.height = height;
        nodes.push_back(std::move(node));
        upper.push_back(up);
        inclusive.push_back(incl ? 1 : 0);
        return id;
      });
  const std::vector<std::uint32_t> first{0};
  merger.merge_class(first, W{1}, W{2}, false, NodeOrigin::kGlobalGrid);
  REQUIRE(nodes.size() == 4);
  CHECK(nodes[3].children == std::vector<NodeId>{0, 1});
  // {0, 1} with the auxiliary vertex: one nonempty part, no new node (unary compression).
  const std::vector<std::uint32_t> unary{3};
  merger.merge_class(unary, W{1}, W{2}, false, NodeOrigin::kGlobalGrid);
  CHECK(nodes.size() == 4);
  // The class {1-3, 3-2}: the auxiliary vertex connects node 3 and leaf 2 (l.1281).
  const std::vector<std::uint32_t> second{1, 2};
  merger.merge_class(second, W{2}, W{4}, false, NodeOrigin::kGlobalGrid);
  REQUIRE(nodes.size() == 5);
  CHECK(nodes[4].children == std::vector<NodeId>{3, 2});
  CHECK(nodes[4].height == 2);
  CHECK(upper == std::vector<W>{0, 0, 0, 2, 4});
  CHECK(nodes[3].parent == 4);
  CHECK(nodes[2].parent == 4);
  merger.merge_class({}, W{9}, W{18}, false, NodeOrigin::kGlobalGrid);  // empty class: nothing
  CHECK(nodes.size() == 5);

  // Arcs: 0 <-> 1 of length 1, 2 -> 1 of length 1 (lighter than the height 2 of node 4, so
  // leaf 2 must precede node 3 among its children), 1 -> 2 of length 3.
  const std::vector<StateArc<W>> arcs{{0, 1, 1}, {1, 0, 1}, {1, 2, 3}, {2, 1, 1}};
  const std::vector<NodeId> leaf_of_state{0, 1, 2};
  const auto lo = core::link_and_order<W>(nodes, arcs, leaf_of_state, upper, inclusive, "test");
  CHECK(lo.links == 2);
  CHECK(lo.roots == std::vector<NodeId>{4});
  REQUIRE(nodes[3].links.size() == 1);
  CHECK(nodes[3].links[0].x == 0);
  CHECK(nodes[3].links[0].y == 1);
  CHECK(nodes[3].links[0].length == 1);
  CHECK(nodes[3].links[0].arc == 0);
  REQUIRE(nodes[4].links.size() == 1);
  CHECK(nodes[4].links[0].arc == 2);  // the first arc of length < 4 between the two children
  CHECK(nodes[4].links[0].length == 3);
  CHECK(nodes[4].children == std::vector<NodeId>{2, 3});  // light arc 2 -> 1
  CHECK(nodes[3].children == std::vector<NodeId>{0, 1});
}

TEMPLATE_TEST_CASE("assemble() with the paper's rule is paper::assemble_local_hierarchy", "",
                   std::uint64_t, double) {
  using W = TestType;
  for (const auto f : {lt::Family::kSparse, lt::Family::kGrid}) {
    for (const std::size_t b : {2U, 5U}) {
      INFO(lt::name(f) << " b=" << b);
      const auto g =
          CsrGraph<W>::from_edges(40, lt::make_graph<W>(f, lt::Weights::kExponential, 40, 5 + b));
      const auto in = build(g, b);
      Diagnostics want_diag;
      Diagnostics got_diag;
      const auto want = sssp::paper::assemble_local_hierarchy(in.od, in.rs, in.joins, &want_diag);
      const auto got = core::assemble(in.od, in.rs, in.joins, core::Rule::paper(40), &got_diag);
      CHECK(hash_of(got) == hash_of(want));
      CHECK(got.result.hierarchy.c0 == 2);  // Lemma 5.10
      CHECK(got_diag.counters == want_diag.counters);
    }
  }
}

TEST_CASE("assemble(): injected constants decide the alignment and the support end") {
  using W = std::uint64_t;
  // Clustered graphs: light clusters (edges 1..3) joined by edges of about 2^54. With N = 32
  // the paper's c = 10 gives N^10 s = 2^50 s < sigma for such blocks: unaligned. An injected
  // exponent 12 gives 2^60 s >= sigma: the same blocks are aligned.
  std::size_t flipped = 0;
  for (const std::size_t b : {2U, 3U, 5U, 8U}) {
    INFO("clustered b=" << b);
    const std::size_t n = 32;
    const auto g = CsrGraph<W>::from_edges(n, lt::clustered(n, n + b));
    const auto in = build(g, b);
    const auto paper = core::assemble(in.od, in.rs, in.joins, core::Rule::paper(n));
    core::Rule rule = core::Rule::paper(n);
    rule.align = core::LengthScale::int_power(n, 12);
    rule.support_end = core::LengthScale::int_power(n, 18);
    const auto wide = core::assemble(in.od, in.rs, in.joins, rule);
    REQUIRE(wide.blocks.size() == paper.blocks.size());
    for (std::size_t u = 0; u < wide.blocks.size(); ++u) {
      const auto& bs = wide.blocks[u];
      if (bs.kind != BlockScaleKind::kAligned && bs.kind != BlockScaleKind::kUnaligned) {
        CHECK(bs.kind == paper.blocks[u].kind);  // singleton / unattached: no alignment test
        continue;
      }
      const W sigma = in.od.corridors[bs.sigma_key].weight;
      CHECK((bs.kind == BlockScaleKind::kAligned) == !(sat::mul_int_pow(bs.s, n, 12) < sigma));
      const auto& ps = paper.blocks[u];
      CHECK((ps.kind == BlockScaleKind::kAligned) == !(sat::mul_int_pow(ps.s, n, 10) < sigma));
      if (bs.kind != ps.kind) ++flipped;
    }
    for (std::size_t c = 0; c < wide.supports.size(); ++c) {
      CHECK(wide.supports[c].upper == sat::mul_int_pow(in.od.corridors[c].weight, n, 18));
      CHECK(paper.supports[c].upper == sat::mul_int_pow(in.od.corridors[c].weight, n, 16));
    }
    const auto report = sssp::paper::validate_ordered_hierarchy(wide.result.hierarchy, g);
    INFO(report.summary());
    CHECK(report.ok());
  }
  CHECK(flipped > 0);
}
