// Unit tests for the Sec 5.4 assembly of the local hierarchy (txt l.1182-1243) and the
// post-processing in the proof of Lemma 5.10 (l.1265-1300): preconditions taken from l.1206
// (N >= 32, h <= N), a hand example derived from the rules of Sec 5.1-5.4, and the grid and
// block-scale invariants on a few fixed inputs. The seeded sweeps with the threshold-SCC oracle
// are in tests/property/paper/test_local_assemble_property.cpp.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/assemble.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::Vertex;
using sssp::paper::assemble_local_hierarchy;
using sssp::paper::BlockScaleKind;
using sssp::paper::compute_joins;
using sssp::paper::local_assembly_applies;
using sssp::paper::local_radial_states;
using sssp::paper::NodeOrigin;
using sssp::paper::owner_decomposition;
using sssp::paper::validate_ordered_hierarchy;
namespace lt = sssp::test::local;

namespace {

template <class W>
CsrGraph<W> unit_path(std::size_t n) {
  std::vector<Edge<W>> e;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    e.push_back({static_cast<Vertex>(i), static_cast<Vertex>(i + 1), W{1}});
  }
  return CsrGraph<W>::from_edges(n, std::move(e));
}

template <class W>
auto run(const CsrGraph<W>& g, std::size_t b, Diagnostics* diag = nullptr) {
  const auto od = owner_decomposition(g, b);
  const auto rs = local_radial_states(g, od);
  const auto j = compute_joins(od, rs);
  return assemble_local_hierarchy(od, rs, j, diag);
}

}  // namespace

TEST_CASE("Sec 5.4 preconditions N >= 32 and h <= N (l.1206)") {
  // In the paper these follow from N > 2 L_N^4; with other parameters the caller must meet them.
  CHECK_FALSE(local_assembly_applies(31, 1));
  CHECK(local_assembly_applies(32, 1));
  CHECK(local_assembly_applies(32, 31));        // h = 32 <= N
  CHECK_FALSE(local_assembly_applies(32, 32));  // h = 33 > N
  CHECK_FALSE(local_assembly_applies(0, 0));
  // The paper's own range: N = 209953, b = L_N^4 = 18^4.
  CHECK(local_assembly_applies(209953, 104976));

  const auto small = unit_path<std::uint64_t>(31);
  const auto od = owner_decomposition(small, 2);
  const auto rs = local_radial_states(small, od);
  const auto j = compute_joins(od, rs);
  CHECK_THROWS_AS(assemble_local_hierarchy(od, rs, j), std::invalid_argument);

  const auto g = unit_path<std::uint64_t>(32);
  const auto od2 = owner_decomposition(g, 32);  // h = 33 > N
  const auto rs2 = local_radial_states(g, od2);
  const auto j2 = compute_joins(od2, rs2);
  CHECK_THROWS_AS(assemble_local_hierarchy(od2, rs2, j2), std::invalid_argument);
  // Joins of a different state graph.
  const auto od3 = owner_decomposition(g, 3);
  const auto rs3 = local_radial_states(g, od3);
  CHECK_THROWS_AS(assemble_local_hierarchy(od3, rs3, j2), std::invalid_argument);
}

TEMPLATE_TEST_CASE("Sec 5.4 hand example: one pivot, no macro, every block is its own root", "",
                   std::uint64_t, double) {
  using W = TestType;
  // Path 0-1-...-31 with unit weights and b = 31 (h = 32 = N). Derivation from the paper:
  //  * Sec 5.1: the postorder residual of vertex 1 reaches 31, so A = {1}; r_v = 1 for v != 1.
  //  * Sec 5.2: one cell, so no macro edge; 0-1 and 1-2 are pivot edges; the near component
  //    {2..31} (and {0} alone) has anchor 1 and every vertex is in shell 0. Copies: in_0, out_0,
  //    out_1 per positive vertex, in/out at the pivot: 31 * 3 + 2 = 95 states. Blocks: the
  //    inward window 0, outward window 0 and outward window 1 of {2..31} (30 states each) and
  //    five singletons.
  //  * Sec 5.3: theta is infinite everywhere (no macro), there are no corridors, so the joins
  //    are the 29 unit block-MST edges of each large block.
  //  * Sec 5.4: each large block has no finite attachment, so it uses its own grid from
  //    s_U = 1 with t_U = s_U (l.1213, l.1236); no Y below t_U has two states; the class [1, 2)
  //    merges all 30 states at height 1. With no macro the state graph's SCCs are the blocks
  //    (Lemma 5.7: the macro-free block quotient is a DAG), so there are 8 roots.
  const auto g = unit_path<W>(32);
  Diagnostics diag;
  const auto a = run(g, 31, &diag);
  const auto& h = a.result.hierarchy;
  REQUIRE(h.num_states() == 95);
  CHECK(h.c0 == 2);
  const auto report = validate_ordered_hierarchy(h, g);
  INFO(report.summary());
  CHECK(report.ok());

  CHECK(h.roots.size() == 8);
  std::size_t internal = 0;
  for (std::size_t v = 0; v < h.nodes.size(); ++v) {
    if (h.is_leaf(static_cast<sssp::paper::NodeId>(v))) continue;
    ++internal;
    CHECK(h.nodes[v].parent == sssp::paper::kNoNode);
    CHECK(h.nodes[v].children.size() == 30);
    CHECK(h.nodes[v].height == W{1});
    CHECK(h.nodes[v].links.size() == 29);
    for (const auto& l : h.nodes[v].links) CHECK(l.length == W{1});
    CHECK(a.node_origin[v] == NodeOrigin::kSeparateGrid);
    CHECK(a.node_upper[v] == W{2});
    CHECK(a.node_upper_inclusive[v] == 0);
  }
  CHECK(internal == 3);

  std::size_t unattached = 0;
  std::size_t singletons = 0;
  for (const auto& bs : a.blocks) {
    if (bs.kind == BlockScaleKind::kSingleton) {
      ++singletons;
      continue;
    }
    CHECK(bs.kind == BlockScaleKind::kUnattached);
    CHECK(bs.size == 30);
    CHECK(bs.s == W{1});
    CHECK(bs.t == W{1});
    CHECK(bs.z == sssp::paper::kNoState);
    ++unattached;
  }
  CHECK(unattached == 3);
  CHECK(singletons == 5);
  CHECK(a.grid.empty());  // no named corridor key: no global grid
  CHECK(a.result.native_entry.size() == 32);
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["local_assemble_calls"] == 1);
  CHECK(diag.counters["local_assemble_unattached_blocks"] == 3);
#endif
  CHECK(diag.fallbacks.empty());
}

TEMPLATE_TEST_CASE("Sec 5.4 global grid: doubling classes over the named corridor keys", "",
                   std::uint64_t, double) {
  using W = TestType;
  for (const auto f : {lt::Family::kSparse, lt::Family::kGrid, lt::Family::kRandomTree}) {
    for (const auto wd : {lt::Weights::kUniform, lt::Weights::kExponential}) {
      for (const std::size_t b : {1U, 2U, 4U}) {
        const auto g = CsrGraph<W>::from_edges(48, lt::make_graph<W>(f, wd, 48, 11 + b));
        const auto od = owner_decomposition(g, b);
        const auto rs = local_radial_states(g, od);
        const auto j = compute_joins(od, rs);
        const auto a = assemble_local_hierarchy(od, rs, j);
        INFO(lt::name(f) << " " << lt::name(wd) << " b=" << b);
        // Q-37: one support per corridor (named key), each starting at or below its key.
        REQUIRE(a.supports.size() == od.corridors.size());
        for (std::size_t c = 0; c < od.corridors.size(); ++c) {
          CHECK(a.supports[c].lower <= od.corridors[c].weight);
          CHECK(od.corridors[c].weight <= a.supports[c].upper);
          REQUIRE(a.supports[c].lower_index < a.grid.size());
          CHECK(a.grid[a.supports[c].lower_index] <= a.supports[c].lower);
        }
        // Classes [g_i, upper_i) partition an increasing range; upper_i <= 2 g_i.
        REQUIRE(a.grid.size() == a.grid_upper.size());
        for (std::size_t i = 0; i < a.grid.size(); ++i) {
          CHECK(a.grid[i] > W{0});
          CHECK(a.grid[i] < a.grid_upper[i]);
          CHECK(a.grid_upper[i] <= a.grid[i] + a.grid[i]);
          if (i + 1 < a.grid.size()) CHECK(a.grid_upper[i] <= a.grid[i + 1]);
        }
        CHECK(validate_ordered_hierarchy(a.result.hierarchy, g).ok());
      }
    }
  }
}
