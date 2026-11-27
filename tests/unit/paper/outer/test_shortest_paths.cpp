// The outer construction, Algorithm 8 l.8-15 of arXiv 2611.34911v1 (Sec 8.3, txt l.2167-2196):
// run_outer_construction / outer_construction in src/paper/outer/shortest_paths.hpp.
//
// Acceptance tests of docs/research/04_outer_scales_pipeline.md Sec 6.2 covered here:
//   PL-2  direct call with the paper's constants (single pivot, Q-72)
//   PL-3  the hand examples L0-L4 of Sec 3.5 (expected values written by hand from the note)
//   PL-9  source at a pivot, at a nonpivot outside every atom, inside an atom (Alg 8 l.14)
// plus the seed (l.14), the per-vertex minimum (l.15, Q-73), F from the exact h (Q-71), the
// fallback records and their shares (Q-69), and the argument checks.
// Oracle: Bellman-Ford on the same graph (shares no code with the construction).
#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "paper/local/local_hierarchy.hpp"
#include "paper/model/random.hpp"
#include "paper/outer/scales.hpp"
#include "paper/outer/scales_numeric.hpp"
#include "paper/outer/shortest_paths.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "unit/paper/outer/shortest_paths_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::kNoAtom;
using sssp::paper::kNoState;
using sssp::paper::StateId;
using sssp::paper::model::WordRng;
using sssp::paper::outer::fallback_shares;
using sssp::paper::outer::OuterBlockKind;
using sssp::paper::outer::OuterConstructionParams;
using sssp::paper::outer::OuterSeparator;
using sssp::paper::outer::run_outer_construction;
using sssp::paper::owner::OuterEdgeClass;
namespace tp = sssp::test::pipeline;

namespace {

constexpr std::uint64_t kSeed = 20261127;

/// Alg 8 l.14-15 as stated, re-derived from the stages: one offer of value 0 at the native copy
/// of the native inward state of s; dist[v] = least value over the states projecting to v; no
/// state below the exact distance (Lemma 4.6).
template <class W>
void require_seed_and_minimum(const sssp::paper::outer::OuterConstruction<W>& r,
                              const CsrGraph<W>& g, Vertex s) {
  const std::size_t n = g.num_vertices();
  REQUIRE(r.offers.size() == 1);
  CHECK(r.offers[0].state == r.sub.entry[r.os.states.native_inward[s]]);
  CHECK(r.offers[0].value == W{0});

  const auto oracle = sssp::baseline::bellman_ford(g, s);
  const W inf = WeightTraits<W>::infinity();
  std::vector<W> least(n, inf);
  REQUIRE(r.state_value.size() == r.sub.hierarchy.num_states());
  for (StateId x = 0; x < r.state_value.size(); ++x) {
    const Vertex v = r.os.states.projection[r.sub.old_state[x]];
    REQUIRE(v < n);
    CHECK(r.state_value[x] >= oracle.dist[v]);
    least[v] = std::min(least[v], r.state_value[x]);
  }
  REQUIRE(r.dist == least);
  for (Vertex v = 0; v < n; ++v) REQUIRE(r.dist[v] != inf);
}

}  // namespace

TEMPLATE_TEST_CASE("PL-3: forced outer construction on the levels L0-L4",
                   "[paper][outer][pipeline]", std::uint64_t, double) {
  using W = TestType;
  for (const auto& lv : tp::levels<W>()) {
    INFO(lv.name);
    const auto g = CsrGraph<W>::from_edges(lv.n, lv.edges);
    for (Vertex s = 0; s < lv.n; ++s) {
      WordRng rng(kSeed);
      Diagnostics diag;
      const auto r = run_outer_construction(g, s, rng, tp::constant_separator(lv.b), &diag);
      if (s == 0) REQUIRE(r.dist == lv.dist_from_0);
      tp::require_oracle(g, s, r.dist);
      require_seed_and_minimum(r, g, s);
      // No stage is replaced by a search: only Alg 4's recorded base case on atoms appears.
      const std::string base = sssp::paper::local_base_case_reason(sssp::paper::kPaperLocalParams);
      for (const auto& fb : diag.fallbacks) CHECK(fb.reason == base);

      // The feature that defines the level (note 04 Sec 3.5).
      const std::size_t macros = tp::count_class(r.cls, OuterEdgeClass::kMacro);
      if (lv.name == "L0") {
        CHECK(r.od.pivots == std::vector<Vertex>{0});
        CHECK(macros == 0);
      } else if (lv.name == "L1") {
        CHECK(r.od.pivots == std::vector<Vertex>{0, 1});
        CHECK(r.od.corridors.size() == 1);
        CHECK(macros == 1);
      } else if (lv.name == "L2") {
        CHECK(r.od.pivots == std::vector<Vertex>{0, 2});
        CHECK(macros == 1);
        CHECK(r.os.special_inward[1] != kNoState);  // a nonpivot: owner-tree states exist
        CHECK(r.os.special_inward[3] != kNoState);
        CHECK(r.jr.scaled_distance[r.aux.aux_vertex[3]] == W{128 * 3});  // J_3 = 3
      } else if (lv.name == "L3") {
        CHECK(r.od.pivots == std::vector<Vertex>{0, 2});
        CHECK(r.cls.edge_class[0] == OuterEdgeClass::kMacro);   // 0 - 1
        CHECK(r.cls.edge_class[3] == OuterEdgeClass::kRadial);  // the chord 1 - 3
        // Three nonsingleton blocks {1, 3} (inward window 0, outward windows 0 and 1), s_U = 3.
        // The macro 0 - 1 attaches the inward block and, through the switch, outward window
        // 0; nothing enters outward window 1, so that block has no finite attachment (Q-66).
        std::size_t blocks = 0;
        std::size_t attached = 0;
        for (const auto& b : r.sc.blocks) {
          if (b.size < 2) continue;
          CHECK(b.size == 2);
          CHECK(b.s == W{3});
          ++blocks;
          attached += b.kind == OuterBlockKind::kUnattached ? 0 : 1;
        }
        CHECK(blocks == 3);
        CHECK(attached == 2);
        CHECK(r.sc.blocks[r.os.states.block[r.os.states.native_inward[1]]].kind ==
              OuterBlockKind::kAligned);  // F^2 3 >= 10: the alignment test ran
        CHECK(r.sc.atoms.empty());
      } else {
        CHECK(r.od.pivots == std::vector<Vertex>{0, 3});
        CHECK(macros == 1);
        CHECK(r.cls.edge_class[2] == OuterEdgeClass::kMacro);  // 0 - 3
        // The pair {4, 5} is an atom in three blocks (edge of weight 1 below s_U / 2 = 4); the
        // pair {1, 2} (edge 4 = s_U) is not.
        REQUIRE(r.sc.atoms.size() == 3);
        for (const auto& atom : r.sc.atoms) {
          std::set<Vertex> image;
          for (const StateId x : atom.states) image.insert(r.os.states.projection[x]);
          CHECK(image == std::set<Vertex>{4, 5});
          CHECK(atom.cutoff > W{4});
          CHECK(atom.cutoff <= W{8});
        }
        // F: b = 3, t = 2 pivots, the paper's P = 2^17: log2 F = 48 (note 04 Sec 3.4).
        CHECK(r.cls.p_log2 == 17);
        CHECK(r.factor.log2 == 48);
        // One Alg 4 base case per Sec 8.2 atom (2 vertices, 1 edge each) besides Alg 7's.
        std::size_t pairs = 0;
        for (const auto& fb : diag.fallbacks) pairs += (fb.vertices == 2 && fb.edges == 1) ? 1 : 0;
        CHECK(pairs >= 3);
      }
    }
  }
}

TEMPLATE_TEST_CASE("PL-9: the source may be a pivot, a plain nonpivot or inside an atom",
                   "[paper][outer][pipeline]", std::uint64_t, double) {
  using W = TestType;
  const auto lv = tp::levels<W>().back();  // L4
  const auto g = CsrGraph<W>::from_edges(lv.n, lv.edges);
  std::size_t pivot = 0;
  std::size_t plain = 0;
  std::size_t in_atom = 0;
  for (Vertex s = 0; s < lv.n; ++s) {
    WordRng rng(kSeed);
    const auto r = run_outer_construction(g, s, rng, tp::constant_separator(lv.b));
    tp::require_oracle(g, s, r.dist);
    const StateId native = r.os.states.native_inward[s];
    if (r.od.is_pivot(s)) {
      ++pivot;
    } else if (r.sc.state_atom[native] == kNoAtom) {
      ++plain;
      CHECK(r.offers[0].state == r.sub.entry[native]);
    } else {
      ++in_atom;
      // l.14 inside an atom: the offer sits at a state of the atom's capped hierarchy whose
      // underlying state is the native inward copy of s.
      CHECK(r.sub.old_state[r.offers[0].state] == native);
    }
  }
  CHECK(pivot == 2);
  CHECK(plain == 2);
  CHECK(in_atom == 2);
}

TEMPLATE_TEST_CASE("PL-2: outer construction with the paper's constants has one pivot",
                   "[paper][outer][pipeline]", std::uint64_t, double) {
  using W = TestType;
  const std::vector<tp::Input<W>> inputs{
      {"single", 1, {}},
      {"edge", 2, {{0, 1, W{5}}}},
      {"path", 9, sssp::test::path_edges<W>(9, [](std::size_t i) { return W(1 + i % 3); })},
      {"grid", 30, sssp::test::grid_edges<W>(5, 6, W{2})},
      {"sparse", 60, tp::connected_edges<W>(60, 40, 1, 50, 7)},
      {"wide", 120, tp::connected_edges<W>(120, 200, 1, 1000000, 8)},
      {"dense", 40, tp::connected_edges<W>(40, 300, 1, 9, 9)},
      {"large", 300, tp::connected_edges<W>(300, 500, 1, 4096, 10)},
  };
  for (const auto& in : inputs) {
    INFO(in.name);
    const auto g = CsrGraph<W>::from_edges(in.n, in.edges);
    for (const Vertex s :
         {Vertex{0}, static_cast<Vertex>(in.n / 2), static_cast<Vertex>(in.n - 1)}) {
      WordRng rng(kSeed);
      Diagnostics diag;
      const auto r =
          run_outer_construction(g, s, rng, sssp::paper::outer::kPaperOuterParams, &diag);
      tp::require_oracle(g, s, r.dist);
      // Sec 7.1 l.1644: b = L_N^64 exceeds N, so the marking leaves one pivot (Q-72): no
      // interowner edge, no macro, every nonsingleton block unattached.
      CHECK(r.separator == OuterSeparator{sssp::paper::model::L(in.n), 64});
      CHECK(r.od.pivots.size() == 1);
      CHECK(tp::count_class(r.cls, OuterEdgeClass::kMacro) == 0);
      CHECK(tp::count_blocks(r.sc, OuterBlockKind::kAligned) == 0);
      CHECK(tp::count_blocks(r.sc, OuterBlockKind::kUnaligned) == 0);
      // Q-71: F from the exact multiword b, not from the saturated hop bound.
      CHECK(r.factor.log2 == sssp::paper::outer::paper_outer_scale_factor(r.cls.p_log2, in.n).log2);
      // Not a fallback: nothing with an outer reason is recorded.
      for (const auto& fb : diag.fallbacks) CHECK(fb.reason.rfind("outer:", 0) != 0);
#if SSSP_DIAGNOSTICS
      CHECK(diag.counters.at("outer.constructions") == 1);
      CHECK(diag.counters.at("outer.pivots") == 1);
      CHECK(diag.counters.at("outer.vertices") == in.n);
#endif
    }
  }
}

TEST_CASE("outer construction: F of a constant separator comes from the exact h = b + 1",
          "[paper][outer][pipeline]") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(40, tp::connected_edges<W>(40, 60, 1, 30, 3));
  for (const std::uint64_t b : {1u, 2u, 3u, 4u, 8u, 64u}) {
    WordRng rng(kSeed);
    const auto r = run_outer_construction(g, Vertex{0}, rng, tp::constant_separator(b));
    CHECK(r.od.separator == b);
    CHECK(r.factor.log2 == sssp::paper::outer::outer_scale_factor(
                               r.cls.p_log2, b + 1, sssp::paper::outer::kOuterShellRatioLog2)
                               .log2);
    tp::require_oracle(g, Vertex{0}, r.dist);
  }
}

TEST_CASE("outer_construction returns the distances of run_outer_construction",
          "[paper][outer][pipeline]") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(25, tp::connected_edges<W>(25, 30, 1, 20, 5));
  WordRng a(kSeed);
  WordRng b(kSeed);
  Diagnostics da;
  Diagnostics db;
  const auto p = tp::constant_separator(2);
  const auto dist = sssp::paper::outer::outer_construction(g, Vertex{4}, a, p, &da);
  const auto r = run_outer_construction(g, Vertex{4}, b, p, &db);
  CHECK(dist == r.dist);
  CHECK(da.fallbacks.size() == db.fallbacks.size());
  // The generator is the caller's and is advanced identically.
  CHECK(a.next_word() == b.next_word());
}

TEST_CASE("outer construction rejects invalid arguments instead of searching",
          "[paper][outer][pipeline]") {
  using W = std::uint64_t;
  const auto g = CsrGraph<W>::from_edges(3, {{0, 1, 2}, {1, 2, 3}});
  WordRng rng(kSeed);
  Diagnostics diag;
  CHECK_THROWS_AS(run_outer_construction(g, Vertex{3}, rng, tp::constant_separator(1), &diag),
                  std::out_of_range);
  OuterConstructionParams zero;
  zero.separator = OuterSeparator{0, 1};
  CHECK_THROWS_AS(run_outer_construction(g, Vertex{0}, rng, zero, &diag), std::invalid_argument);
  zero.separator = OuterSeparator{2, 0};
  CHECK_THROWS_AS(run_outer_construction(g, Vertex{0}, rng, zero, &diag), std::invalid_argument);
  // Not connected (Alg 8 l.2 retains the source component): rejected, not answered.
  const auto split = CsrGraph<W>::from_edges(4, {{0, 1, 2}, {2, 3, 3}});
  CHECK_THROWS_AS(run_outer_construction(split, Vertex{0}, rng, tp::constant_separator(1), &diag),
                  std::invalid_argument);
  // A zero weight (not contracted): rejected.
  const auto zero_edge = CsrGraph<W>::from_edges(2, {{0, 1, 0}});
  CHECK_THROWS_AS(
      run_outer_construction(zero_edge, Vertex{0}, rng, tp::constant_separator(1), &diag),
      std::invalid_argument);
  CHECK(diag.fallbacks.empty());
}

TEST_CASE("Q-69: fallback shares are summed per reason and never across reasons",
          "[paper][outer][pipeline][fallback]") {
  Diagnostics d;
  d.fallbacks.push_back({"local: a", 2, 1});
  d.fallbacks.push_back({"name_routes: b", 0, 7});
  d.fallbacks.push_back({"local: a", 3, 4});
  d.fallbacks.push_back({"local: a", 1, 0});
  const auto shares = fallback_shares(d);
  REQUIRE(shares.size() == 2);
  CHECK(shares[0].reason == "local: a");
  CHECK(shares[0].records == 3);
  CHECK(shares[0].vertices == 6);
  CHECK(shares[0].edges == 5);
  CHECK(shares[1].reason == "name_routes: b");
  CHECK(shares[1].records == 1);
  CHECK(shares[1].vertices == 0);
  CHECK(shares[1].edges == 7);
  CHECK(fallback_shares(Diagnostics{}).empty());

  // On L4 the vertices 4 and 5 lie in three Sec 8.2 atoms and in Alg 7's atoms: the vertex sum
  // of the single local reason exceeds what one cover of the graph would give.
  using W = std::uint64_t;
  const auto lv = tp::levels<W>().back();
  const auto g = CsrGraph<W>::from_edges(lv.n, lv.edges);
  WordRng rng(kSeed);
  Diagnostics diag;
  (void)run_outer_construction(g, Vertex{0}, rng, tp::constant_separator(lv.b), &diag);
  const auto l4 = fallback_shares(diag);
  REQUIRE(l4.size() == 1);
  CHECK(l4[0].reason == sssp::paper::local_base_case_reason(sssp::paper::kPaperLocalParams));
  CHECK(l4[0].records == diag.fallbacks.size());
  CHECK(l4[0].vertices >= 3 * 2);
}
