// Numerical scales and atoms of the outer construction: arXiv 2611.34911v1, Sec 8.2 (txt
// l.2076-2104) and Lemma 8.4 (l.2105-2142). Stages S0-S7 of
// docs/research/04_outer_scales_pipeline.md Sec 2: new names, F from the exact h, block scales,
// alignment, supports and the global grid, the cutoffs t_U, atoms and the classes of the joins.
//
// Expected values are derived by hand from the paper's rules (see the comments at each case and
// in outer_scales_test_support.hpp), never copied from a run.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <string>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/atom_replacement.hpp"
#include "paper/local/threshold_joins.hpp"
#include "paper/model/params.hpp"
#include "paper/outer/outer_states.hpp"
#include "paper/outer/scales.hpp"
#include "paper/outer/scales_numeric.hpp"
#include "paper/owner/aux_metric.hpp"
#include "sssp/graph.hpp"
#include "unit/paper/outer/outer_scales_test_support.hpp"
#include "unit/paper/outer/outer_states_test_support.hpp"

using sssp::Edge;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::JoinKind;
using sssp::paper::kNoAtom;
using sssp::paper::StateId;
using sssp::paper::outer::outer_alignment_by_median_pruning;
using sssp::paper::outer::outer_block_aligned;
using sssp::paper::outer::outer_name_weights;
using sssp::paper::outer::outer_scale_factor;
using sssp::paper::outer::outer_scale_factor_power_separator;
using sssp::paper::outer::OuterBlockKind;
using sssp::paper::outer::OuterEventGrid;
using sssp::paper::outer::OuterScaleFactor;
using sssp::paper::outer::paper_outer_scale_factor;
namespace to = sssp::test::outer;

namespace {

// ---- a second multiword routine for SC-1 (little-endian 32-bit limbs), written here.
using Big = std::vector<std::uint32_t>;

Big big_from(std::uint64_t x) {
  Big r;
  for (; x != 0; x >>= 32) r.push_back(static_cast<std::uint32_t>(x));
  return r;
}

Big big_add(const Big& a, const Big& b) {
  Big r;
  std::uint64_t carry = 0;
  for (std::size_t i = 0; i < std::max(a.size(), b.size()) || carry != 0; ++i) {
    std::uint64_t cur = carry;
    if (i < a.size()) cur += a[i];
    if (i < b.size()) cur += b[i];
    r.push_back(static_cast<std::uint32_t>(cur));
    carry = cur >> 32;
  }
  return r;
}

Big big_mul(const Big& a, const Big& b) {
  Big r(a.size() + b.size(), 0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    std::uint64_t carry = 0;
    for (std::size_t j = 0; j < b.size() || carry != 0; ++j) {
      std::uint64_t cur = r[i + j] + carry;
      if (j < b.size()) cur += static_cast<std::uint64_t>(a[i]) * b[j];
      r[i + j] = static_cast<std::uint32_t>(cur);
      carry = cur >> 32;
    }
  }
  while (!r.empty() && r.back() == 0) r.pop_back();
  return r;
}

bool big_less(const Big& a, const Big& b) {
  if (a.size() != b.size()) return a.size() < b.size();
  for (std::size_t i = a.size(); i-- > 0;) {
    if (a[i] != b[i]) return a[i] < b[i];
  }
  return false;
}

Big big_pow2(std::uint64_t e) {
  Big r(e / 32 + 1, 0);
  r.back() = std::uint32_t{1} << (e % 32);
  return r;
}

/// X = 2 P (2h+1) (8 K^2 h)^2 with P = 2^p, K = 16, h = b + 1 (Sec 8.2 l.2080).
Big scale_product(unsigned p, const Big& b) {
  const Big h = big_add(b, big_from(1));
  const Big inner = big_mul(big_from(8 * 16 * 16), h);
  return big_mul(big_mul(big_pow2(p + 1), big_add(big_add(h, h), big_from(1))),
                 big_mul(inner, inner));
}

Big big_power(std::uint64_t base, unsigned e) {
  Big r = big_from(1);
  for (unsigned i = 0; i < e; ++i) r = big_mul(r, big_from(base));
  return r;
}

/// F = 2^f is the least power of two strictly above X: 2^(f-1) <= X < 2^f.
void check_least_power(std::uint64_t f, const Big& x) {
  CHECK(big_less(x, big_pow2(f)));
  CHECK_FALSE(big_less(x, big_pow2(f - 1)));
}

template <class W>
constexpr W kInf = WeightTraits<W>::infinity();

/// The join records of kind `kind` whose first endpoint lies in block `u`.
template <class W>
std::vector<std::uint32_t> joins_of_block(const to::ScalesCase<W>& r, JoinKind kind,
                                          std::uint32_t u) {
  std::vector<std::uint32_t> out;
  const auto& jr = r.os.joins.joins;
  for (std::uint32_t i = 0; i < jr.size(); ++i) {
    if (jr[i].kind == kind && r.os.states.block[jr[i].a] == u) out.push_back(i);
  }
  return out;
}

template <class W>
std::size_t count_kind(const to::ScalesCase<W>& r, OuterBlockKind k) {
  std::size_t c = 0;
  for (const auto& b : r.sc.blocks) c += b.kind == k ? 1 : 0;
  return c;
}

}  // namespace

// ------------------------------------------------------------------------------ S1: F (SC-1, SC-2)

TEST_CASE("outer scales: scale_factor_log2 matches big-integer values (SC-1)",
          "[paper][outer][scales]") {
  // Note 04 Sec 1.4. b = 589^64 is the paper's separator at the least N that takes the outer
  // construction (t = 1, p = 17); b = 20^64 the paper's b at N = 10^6.
  CHECK(outer_scale_factor_power_separator(17, 589, 64, 4).log2 == 1808);
  check_least_power(1808, scale_product(17, big_power(589, 64)));
  CHECK(outer_scale_factor_power_separator(22, 589, 64, 4).log2 == 1813);
  CHECK(outer_scale_factor_power_separator(17, 590, 64, 4).log2 == 1809);
  check_least_power(1809, scale_product(17, big_power(590, 64)));
  CHECK(outer_scale_factor_power_separator(17, 10, 64, 4).log2 == 679);
  CHECK(outer_scale_factor_power_separator(17, 20, 64, 4).log2 == 871);
  check_least_power(871, scale_product(17, big_power(20, 64)));
  CHECK(outer_scale_factor_power_separator(17, 64, 64, 4).log2 == 1194);
  check_least_power(1194, scale_product(17, big_power(64, 64)));
  // One-word separators agree with the one-word routine: b = 3, p = 17 (48) and b = 160000,
  // p = 1 (77).
  CHECK(outer_scale_factor_power_separator(17, 3, 1, 4).log2 == 48);
  CHECK(outer_scale_factor(17, 4, 4).log2 == 48);
  check_least_power(48, scale_product(17, big_from(3)));
  CHECK(outer_scale_factor_power_separator(1, 160000, 1, 4).log2 == 77);
  CHECK(outer_scale_factor(1, 160001, 4).log2 == 77);
  check_least_power(77, scale_product(1, big_from(160000)));
  CHECK(outer_scale_factor_power_separator(1, 400, 2, 4).log2 == 77);  // 400^2 = 160000
  // The value saturates, the exponent does not (Q-06).
  const auto f = outer_scale_factor_power_separator(17, 589, 64, 4);
  CHECK(f.value == std::numeric_limits<std::uint64_t>::max());
  CHECK_FALSE(f.value_exact());
  CHECK(outer_scale_factor_power_separator(0, 3, 1, 4).value == (std::uint64_t{1} << 31));

  std::mt19937_64 rng(8204);
  for (int it = 0; it < 200; ++it) {
    const auto p = static_cast<unsigned>(rng() % 30);
    const std::uint64_t base = 1 + rng() % 700;
    const auto e = static_cast<unsigned>(rng() % 70);
    CAPTURE(p, base, e);
    const auto got = outer_scale_factor_power_separator(p, base, e, 4);
    check_least_power(got.log2, scale_product(p, big_power(base, e)));
  }
  CHECK_THROWS_AS(outer_scale_factor_power_separator(17, 0, 3, 4), std::invalid_argument);
}

TEST_CASE("outer scales: the scale factor is not derived from a saturated h (SC-2)",
          "[paper][outer][scales]") {
  // N = 10^6: L_N = 20, b = 20^64 (Sec 7.1 l.1644), which the decomposition holds saturated.
  REQUIRE(sssp::paper::model::L(1000000) == 20);
  CHECK(paper_outer_scale_factor(17, 1000000).log2 == 871);
  // What the saturated word would give (note 04 Sec 1.3): a different F.
  const std::uint64_t saturated = sssp::paper::owner::outer_separator_size(1000000);
  REQUIRE(saturated == std::numeric_limits<std::uint64_t>::max());
  CHECK(outer_scale_factor(17, saturated, 4).log2 == 233);
  // The least N of the outer branch does not fit a word; the largest word still has L = 65.
  CHECK(paper_outer_scale_factor(17, std::numeric_limits<std::uint64_t>::max()).log2 ==
        outer_scale_factor_power_separator(17, 65, 64, 4).log2);
  check_least_power(paper_outer_scale_factor(17, 6).log2, scale_product(17, big_power(3, 64)));
}

// ----------------------------------------------------------------------- S3: alignment (SC-3,
// SC-13)

TEST_CASE("outer scales: F^2 s >= a is decided exactly around saturation (SC-3)",
          "[paper][outer][scales]") {
  __extension__ typedef unsigned __int128 Wide;
  std::mt19937_64 rng(77);
  for (int it = 0; it < 4000; ++it) {
    const std::uint64_t log2_f = rng() % 40;
    const OuterScaleFactor f{log2_f, std::uint64_t{1} << log2_f};
    const std::uint64_t s = 1 + (rng() >> (rng() % 64));
    std::uint64_t a = 1 + (rng() >> (rng() % 64));
    if (it % 3 == 0 && 2 * log2_f < 64 && s <= (kInf<std::uint64_t> >> (2 * log2_f)) - 1) {
      a = (s << (2 * log2_f)) + (rng() % 3) - 1;  // the boundary and its neighbours
    }
    if (a == 0 || a == kInf<std::uint64_t> || s == kInf<std::uint64_t>) continue;
    // Exact: F^2 s >= a iff s >= ceil(a / F^2); 2 f < 128 here.
    const bool exact = 2 * log2_f >= 127 || ((static_cast<Wide>(s) << (2 * log2_f)) >= a);
    CAPTURE(log2_f, s, a);
    REQUIRE(outer_block_aligned<std::uint64_t>(f, s, a) == exact);
  }
  // Dyadic doubles: 2^i 2^(2f) >= 2^j iff i + 2f >= j, also beyond the range of double.
  for (int i = -1000; i <= 1000; i += 37) {
    for (int j = -1000; j <= 1000; j += 41) {
      for (const std::uint64_t log2_f :
           {std::uint64_t{28}, std::uint64_t{48}, std::uint64_t{871}, std::uint64_t{1808}}) {
        const OuterScaleFactor f{log2_f, std::numeric_limits<std::uint64_t>::max()};
        const bool exact = static_cast<std::int64_t>(i) + 2 * static_cast<std::int64_t>(log2_f) >=
                           static_cast<std::int64_t>(j);
        REQUIRE(outer_block_aligned<double>(f, std::ldexp(1.0, i), std::ldexp(1.0, j)) == exact);
      }
    }
  }
}

TEMPLATE_TEST_CASE("outer scales: alignment by median pruning equals the per-block test (SC-13)",
                   "[paper][outer][scales]", std::uint64_t, double) {
  using W = TestType;
  std::mt19937_64 rng(1313);
  for (int it = 0; it < 600; ++it) {
    const std::uint64_t log2_f = 1 + rng() % 20;
    const OuterScaleFactor f{log2_f, std::uint64_t{1} << log2_f};
    const std::size_t n = rng() % 40;
    const auto a = static_cast<W>(std::uint64_t{1} << (10 + rng() % 40));
    std::vector<W> s;
    for (std::size_t i = 0; i < n; ++i) {
      // Scales around a / F^2, with ties.
      s.push_back(static_cast<W>(1 + (rng() % 9) * (std::uint64_t{1} << (rng() % 30))));
    }
    std::uint64_t tests = 0;
    const auto got = outer_alignment_by_median_pruning<W>(f, a, std::span<const W>(s), &tests);
    REQUIRE(got.size() == n);
    for (std::size_t i = 0; i < n; ++i) {
      CAPTURE(log2_f, n, i);
      REQUIRE((got[i] != 0) == outer_block_aligned<W>(f, s[i], a));
    }
    // l.2130 / Sec 6.1: each tested median at least halves the group.
    std::uint64_t bound = 0;
    for (std::size_t k = n; k != 0; k >>= 1) ++bound;  // floor(log2 n) + 1; 0 for n = 0
    CHECK(tests <= bound);
  }
}

// --------------------------------------------------------------------------- S0: new names (SC-4)

TEMPLATE_TEST_CASE("outer scales: new names are the least macro weight (SC-4)",
                   "[paper][outer][scales]", std::uint64_t, double) {
  using W = TestType;
  using sssp::paper::owner::kNoName;
  using sssp::paper::owner::OuterEdgeClass;
  // Six edges on five vertices; four macros with two names, a radial interowner edge that
  // carries the name 0 with a lighter weight (not a macro: ignored, Q-65), and an intraowner
  // edge. Dictionary index 2 names no macro.
  const std::vector<Edge<W>> edges{{0, 1, W{40}}, {1, 2, W{24}}, {2, 3, W{90}},
                                   {3, 4, W{70}}, {0, 4, W{5}},  {1, 3, W{2}}};
  const auto g = sssp::CsrGraph<W>::from_edges(5, edges);
  sssp::paper::owner::OuterEdgeClassification<W> cls;
  cls.dictionary = {W{30}, W{80}, W{200}};
  cls.edge_class = {OuterEdgeClass::kMacro, OuterEdgeClass::kMacro,  OuterEdgeClass::kMacro,
                    OuterEdgeClass::kMacro, OuterEdgeClass::kRadial, OuterEdgeClass::kRadial};
  cls.name = {0, 0, 1, 1, 0, kNoName};
  const auto a = outer_name_weights(g, cls);
  REQUIRE(a.size() == 3);
  CHECK(a[0] == W{24});  // min(40, 24); the radial edge of weight 5 does not count
  CHECK(a[1] == W{70});  // min(90, 70)
  CHECK(a[2] == kInf<W>);
  cls.name[3] = 7;  // outside the dictionary
  CHECK_THROWS_AS(outer_name_weights(g, cls), std::invalid_argument);
}

// ------------------------------------------------------------------------ hand example L3 (SC-5)

TEMPLATE_TEST_CASE("outer scales: hand example L3, aligned blocks on the global grid (SC-5)",
                   "[paper][outer][scales]", std::uint64_t, double) {
  using W = TestType;
  for (const auto names : {to::Names::kExactP1, to::Names::kNameRoutesPaper}) {
    auto r = to::run_scales(to::make_outer_case<W>(4, to::level3_edges<W>(W{10}), 2, names, 5));
    const auto& st = r.os.states;
    const auto& sc = r.sc;
    REQUIRE(r.c.od.pivots == std::vector<Vertex>{0, 2});
    REQUIRE(r.c.od.radius[1] == W{3});
    REQUIRE(r.c.od.radius[3] == W{3});
    REQUIRE(r.c.cls.edge_class[0] == sssp::paper::owner::OuterEdgeClass::kMacro);
    // b = 2, h = 3: f = 29 with P = 1 and 46 with the paper's P = 2^17 (note 04 Sec 1.4).
    CHECK(sc.factor.log2 == (names == to::Names::kExactP1 ? 29U : 46U));
    CHECK(r.diag.fallbacks.empty());  // S0-S7 never substitute a search

    // S0: one macro (0-1, weight 10), so its name becomes 10.
    const std::uint32_t name = r.c.cls.name[0];
    CHECK(sc.name_weight[name] == W{10});
    const std::uint32_t key = sc.name_key[name];
    REQUIRE(key < sc.num_name_keys);
    REQUIRE(sc.key_weight.size() == sc.num_name_keys + 1);  // one corridor (0 .. 2)
    CHECK(sc.key_weight[key] == W{10});
    CHECK(sc.key_weight[sc.num_name_keys] == W{10});  // the corridor weight b_T(0, 2)

    // Both vertices 1 and 3 have J = 3, shell 0: copies in_0, out_0, out_1, hence three blocks
    // {1, 3} (inward window 0, outward windows 0 and 1), each with the chord of weight 5.
    const std::uint32_t bin = st.block[st.native_inward[1]];
    const std::uint32_t bout = st.block[st.native_outward[1]];
    const std::uint32_t bup = st.block[st.upper_outward[1]];
    REQUIRE(bin == st.block[st.native_inward[3]]);
    REQUIRE(bout == st.block[st.native_outward[3]]);
    REQUIRE(bup == st.block[st.upper_outward[3]]);
    // The macro 0-1 resets to in_0(1) and leaves from every copy of 1: the inward block and
    // (through the switch) the outward window 0 are attached, both at theta = 10 (the macro is
    // the heaviest edge). Nothing enters outward window 1 (no neighbour of a lower shell, no
    // macro head), so I0 is infinite there: that block has no finite attachment.
    for (const std::uint32_t u : {bin, bout}) {
      const auto& b = sc.blocks[u];
      CHECK(b.kind == OuterBlockKind::kAligned);  // F^2 3 >= 10
      CHECK(b.size == 2);
      CHECK(b.s == W{3});
      CHECK(b.s_vertex == 1);  // least id attaining s_U
      CHECK(b.name == name);
      CHECK(b.key == key);
      REQUIRE(b.z != sssp::paper::kNoState);
      CHECK(st.block[b.z] == u);
      CHECK(r.os.joins.theta[b.z] == W{10});
      // The support of the name starts at min(10, s_U) = 3, so the grid is 3, 6, 12, ... and
      // t_U = 3: s_U / 2 < t_U <= s_U.
      CHECK(b.t == W{3});
      const auto mst = joins_of_block(r, JoinKind::kBlockMst, u);
      REQUIRE(mst.size() == 1);
      CHECK(r.os.joins.joins[mst[0]].weight == W{5});
      const auto& ev = sc.events[mst[0]];
      CHECK(ev.grid == OuterEventGrid::kGlobal);  // 5 >= t_U: significant, class [3, 6)
      CHECK(ev.lower == W{3});
      CHECK(ev.upper == W{6});
    }
    {
      const auto& b = sc.blocks[bup];
      CHECK(b.kind == OuterBlockKind::kUnattached);
      CHECK(b.s == W{3});
      CHECK(b.z == sssp::paper::kNoState);
      // Separate grid (Q-66): the piece of vertex 1 has origin 3 = r_1, so t_U = 3.
      CHECK(b.t == W{3});
      const auto mst = joins_of_block(r, JoinKind::kBlockMst, bup);
      REQUIRE(mst.size() == 1);
      const auto& ev = sc.events[mst[0]];
      CHECK(ev.grid == OuterEventGrid::kSeparate);
      CHECK(ev.index == 0);
      CHECK(ev.lower == W{3});
      CHECK(ev.upper == W{6});
    }
    CHECK(count_kind(r, OuterBlockKind::kAligned) == 2);
    CHECK(count_kind(r, OuterBlockKind::kUnattached) == 1);
    CHECK(count_kind(r, OuterBlockKind::kUnaligned) == 0);
    CHECK(count_kind(r, OuterBlockKind::kSingleton) == st.num_blocks - 3);
    CHECK(sc.atoms.empty());
    for (const auto a : sc.state_atom) CHECK(a == kNoAtom);

    // S4: the two supports [3, F^3 10] and [10, F^3 10] overlap: one component from 3.
    REQUIRE(sc.grid.component_begin == std::vector<std::size_t>{0});
    REQUIRE(sc.grid.points.size() >= 3);
    CHECK(sc.grid.points[0] == W{3});
    CHECK(sc.grid.points[1] == W{6});
    CHECK(sc.grid.points[2] == W{12});
    CHECK(sc.grid.supports[key].lower == W{3});
    CHECK(sc.grid.supports[sc.num_name_keys].lower == W{10});

    // S7: every finite attachment is 10 (all bottlenecks are at most the macro weight), class
    // [6, 12); so is the corridor join.
    std::size_t owner_joins = 0;
    std::size_t corridor_joins = 0;
    const auto& jr = r.os.joins.joins;
    for (std::size_t i = 0; i < jr.size(); ++i) {
      if (jr[i].kind == JoinKind::kBlockMst) continue;
      (jr[i].kind == JoinKind::kOwner ? owner_joins : corridor_joins) += 1;
      CHECK(jr[i].weight == W{10});
      CHECK(sc.events[i].grid == OuterEventGrid::kGlobal);
      CHECK(sc.events[i].index == 1);
      CHECK(sc.events[i].lower == W{6});
      CHECK(sc.events[i].upper == W{12});
    }
    CHECK(corridor_joins == 1);
    CHECK(owner_joins >= 4 + 4);  // at least the four pivot copies and the two attached blocks
    CHECK(sssp::paper::outer::check_outer_scales(r.c.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.calendars,
                                                 r.os, sc)
              .empty());
  }
}

// ------------------------------------------------------------------------ hand example L4 (SC-6)

TEMPLATE_TEST_CASE("outer scales: hand example L4, atoms below the cutoff (SC-6, S0-S7)",
                   "[paper][outer][scales]", std::uint64_t, double) {
  using W = TestType;
  for (const auto names : {to::Names::kExactP1, to::Names::kNameRoutesPaper}) {
    auto r = to::run_scales(to::make_outer_case<W>(6, to::level4_edges<W>(), 3, names, 9));
    const auto& st = r.os.states;
    const auto& sc = r.sc;
    REQUIRE(r.c.od.pivots == std::vector<Vertex>{0, 3});
    REQUIRE(r.c.od.radius[4] == W{8});
    REQUIRE(r.c.od.radius[5] == W{8});
    REQUIRE(r.c.od.radius[1] == W{4});
    REQUIRE(r.c.od.radius[2] == W{4});
    CHECK(sc.name_weight[r.c.cls.name[2]] == W{100});

    // The only macro joins the two pivots, so no macro head lies in an inward block of 4, 5 or
    // of 1, 2 and no outward copy can return to a macro tail: every nonsingleton block is
    // unattached and takes its floor on the separate grid of its piece (Q-66).
    CHECK(count_kind(r, OuterBlockKind::kUnattached) == 6);
    CHECK(count_kind(r, OuterBlockKind::kAligned) == 0);
    CHECK(count_kind(r, OuterBlockKind::kUnaligned) == 0);

    // Blocks {4, 5}: s_U = 8 (r_4 = r_5 = 8); the piece origin is 8, so t_U = 8, in (4, 8]. The
    // edge 4-5 of weight 1 is below t_U: one atom per block, H = that edge.
    REQUIRE(sc.atoms.size() == 3);
    std::vector<std::uint32_t> seen;
    for (const auto& [x4, x5] : {std::pair{st.native_inward[4], st.native_inward[5]},
                                 std::pair{st.native_outward[4], st.native_outward[5]},
                                 std::pair{st.upper_outward[4], st.upper_outward[5]}}) {
      const std::uint32_t u = st.block[x4];
      REQUIRE(u == st.block[x5]);
      const auto& b = sc.blocks[u];
      CHECK(b.kind == OuterBlockKind::kUnattached);
      CHECK(b.size == 2);
      CHECK(b.s == W{8});
      CHECK(b.s_vertex == 4);
      CHECK(b.t == W{8});
      const std::uint32_t a = sc.state_atom[x4];
      REQUIRE(a != kNoAtom);
      REQUIRE(a == sc.state_atom[x5]);
      seen.push_back(a);
      const auto& atom = sc.atoms[a];
      CHECK(sc.atom_block[a] == u);
      CHECK(atom.states == std::vector<StateId>{std::min(x4, x5), std::max(x4, x5)});
      REQUIRE(atom.edges.size() == 1);
      CHECK(atom.edges[0].weight == W{1});
      CHECK(atom.edges[0].u + atom.edges[0].v == 1);  // local vertices 0 and 1
      CHECK(atom.spanning_tree == std::vector<sssp::EdgeId>{0});
      CHECK(atom.cutoff == W{8});
      const auto mst = joins_of_block(r, JoinKind::kBlockMst, u);
      REQUIRE(mst.size() == 1);
      CHECK(sc.events[mst[0]].grid == OuterEventGrid::kAtom);
      CHECK(sc.events[mst[0]].index == a);
    }
    std::ranges::sort(seen);
    CHECK(seen == std::vector<std::uint32_t>{0, 1, 2});

    // Blocks {1, 2}: s_U = t_U = 4 and the edge 1-2 of weight 4 is significant (not < t_U):
    // no atom, separate class [4, 8).
    for (const auto& [x1, x2] : {std::pair{st.native_inward[1], st.native_inward[2]},
                                 std::pair{st.native_outward[1], st.native_outward[2]},
                                 std::pair{st.upper_outward[1], st.upper_outward[2]}}) {
      const std::uint32_t u = st.block[x1];
      REQUIRE(u == st.block[x2]);
      CHECK(sc.blocks[u].s == W{4});
      CHECK(sc.blocks[u].t == W{4});
      CHECK(sc.state_atom[x1] == kNoAtom);
      CHECK(sc.state_atom[x2] == kNoAtom);
      const auto mst = joins_of_block(r, JoinKind::kBlockMst, u);
      REQUIRE(mst.size() == 1);
      CHECK(sc.events[mst[0]].grid == OuterEventGrid::kSeparate);
      CHECK(sc.events[mst[0]].lower == W{4});
      CHECK(sc.events[mst[0]].upper == W{8});
    }
    CHECK(sssp::paper::outer::check_outer_scales(r.c.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.calendars,
                                                 r.os, sc)
              .empty());
  }
}

// ----------------------------------------------------------------- unaligned blocks (SC-7, SC-8)

TEST_CASE("outer scales: unaligned block on uint64 (SC-7)", "[paper][outer][scales]") {
  using W = std::uint64_t;
  // L3 with the corridor edge 0-1 at 2^60 and P = 1: f = 29, F^2 s_U = 3 2^58 < 2^60 = a_U.
  const W big = W{1} << 60;
  auto r = to::run_scales(
      to::make_outer_case<W>(4, to::level3_edges<W>(big), 2, to::Names::kExactP1, 3));
  const auto& st = r.os.states;
  const auto& sc = r.sc;
  REQUIRE(r.c.od.pivots == std::vector<Vertex>{0, 2});
  REQUIRE(sc.factor.log2 == 29);
  const std::uint32_t name = r.c.cls.name[0];
  CHECK(sc.name_weight[name] == big);
  for (const StateId x : {st.native_inward[1], st.native_outward[1]}) {
    const std::uint32_t u = st.block[x];
    const auto& b = sc.blocks[u];
    CHECK(b.kind == OuterBlockKind::kUnaligned);
    CHECK(b.s == W{3});
    CHECK(b.name == name);
    // l.2085-2086: every within-window edge is below a_U / 2.
    const auto mst = joins_of_block(r, JoinKind::kBlockMst, u);
    REQUIRE(mst.size() == 1);
    CHECK(r.os.joins.joins[mst[0]].weight == W{5});
    CHECK(W{5} + W{5} < sc.name_weight[b.name]);
    // Lemma 8.4 proof (1): equal attachment threshold on the block.
    for (StateId y = 0; y < st.num_states(); ++y) {
      if (st.block[y] == u) CHECK(r.os.joins.theta[y] == big);
    }
    // t_U on the piece grid: the piece of vertex 1 has origin 3.
    CHECK(b.t == W{3});
    const auto& ev = sc.events[mst[0]];
    CHECK(ev.grid == OuterEventGrid::kSeparate);
    CHECK(ev.lower == W{3});
    CHECK(ev.upper == W{6});
    // l.2093: the block is complete before its first global event: 2 p_top <= g.
    const auto own = joins_of_block(r, JoinKind::kOwner, u);
    REQUIRE(own.size() == 2);
    for (const auto j : own) {
      CHECK(sc.events[j].grid == OuterEventGrid::kGlobal);
      CHECK(sc.events[j].lower == big);  // the grid starts at the name 2^60
      CHECK(ev.lower + ev.lower <= sc.events[j].lower);
    }
  }
  CHECK(sc.blocks[st.block[st.upper_outward[1]]].kind == OuterBlockKind::kUnattached);
  CHECK(sc.atoms.empty());
  // No aligned block lowers a support: the grid starts at 2^60.
  REQUIRE(!sc.grid.points.empty());
  CHECK(sc.grid.points[0] == big);
  CHECK(sssp::paper::outer::check_outer_scales(r.c.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.calendars,
                                               r.os, sc)
            .empty());
}

TEST_CASE("outer scales: unaligned block with an atom (SC-8)", "[paper][outer][scales]") {
  using W = std::uint64_t;
  // The L4 tree with 0-3 at 2^63 + 1, 3-4 at 2, 4-5 at 1, plus the nontree edge 2-5 of weight
  // 2^63 + 2 (kept: it is not longer than its tree path), b = 3, P = 1: f = 31.
  // [Correction to note 04 Sec 6.1, SC-8] Without the edge 2-5 the only macro joins the two
  // pivots and no nonsingleton block has a finite attachment (see SC-6); 2-5 is a macro with a
  // nonpivot head in each cell. Both macros get the name 2^63 + 1 (the least instance weight
  // that is at least their pivot distance), so the new name is a = 2^63 + 1.
  const W corridor = (W{1} << 63) + 1;
  const std::vector<Edge<W>> edges{{0, 1, W{4}}, {1, 2, W{4}}, {0, 3, corridor},
                                   {3, 4, W{2}}, {4, 5, W{1}}, {2, 5, corridor + 1}};
  auto r = to::run_scales(to::make_outer_case_named<W>(6, edges, 3, {corridor}, {0, 0}));
  const auto& st = r.os.states;
  const auto& sc = r.sc;
  REQUIRE(r.c.od.pivots == std::vector<Vertex>{0, 3});
  REQUIRE(r.c.cls.edge_class[2] == sssp::paper::owner::OuterEdgeClass::kMacro);
  REQUIRE(r.c.cls.edge_class[5] == sssp::paper::owner::OuterEdgeClass::kMacro);
  REQUIRE(r.c.od.radius[4] == W{2});
  REQUIRE(r.c.od.radius[5] == W{2});
  REQUIRE(sc.factor.log2 == 31);
  CHECK(sc.name_weight[0] == corridor);

  // Blocks {4, 5}: s_U = 2 and F^2 2 = 2^63 < a_U. The inward block (macro head in_0(5)) and
  // outward window 0 (through the switch) are attached, hence unaligned; outward window 1 is
  // unattached. In all three t_U = 2 (piece origin 2) and the edge of weight 1 forms an atom.
  REQUIRE(sc.atoms.size() == 3);
  const std::pair<StateId, OuterBlockKind> blocks45[] = {
      {st.native_inward[4], OuterBlockKind::kUnaligned},
      {st.native_outward[4], OuterBlockKind::kUnaligned},
      {st.upper_outward[4], OuterBlockKind::kUnattached}};
  for (const auto& [x, kind] : blocks45) {
    const auto& b = sc.blocks[st.block[x]];
    CHECK(b.kind == kind);
    CHECK(b.size == 2);
    CHECK(b.s == W{2});
    CHECK(b.t == W{2});
    const std::uint32_t a = sc.state_atom[x];
    REQUIRE(a != kNoAtom);
    CHECK(sc.atoms[a].states.size() == 2);  // the whole block is one atom
    CHECK(sc.atoms[a].cutoff == W{2});
    REQUIRE(sc.atoms[a].edges.size() == 1);
    CHECK(sc.atoms[a].edges[0].weight == W{1});
    if (kind == OuterBlockKind::kUnaligned) {
      CHECK(sc.name_weight[b.name] == corridor);
      // The first global event of the atom is far above its cutoff: class [2^63, infinity).
      for (const auto j : joins_of_block(r, JoinKind::kOwner, st.block[x])) {
        CHECK(sc.events[j].grid == OuterEventGrid::kGlobal);
        CHECK(sc.events[j].lower == W{1} << 63);
        CHECK(sc.events[j].upper == kInf<W>);
      }
    }
  }
  // Blocks {1, 2}: s_U = 4 and F^2 4 = 2^64 >= a_U: the attached ones are aligned, their name's
  // support starts at 4, t_U = 4, and the edge 1-2 of weight 4 is in the global class [4, 8).
  const std::pair<StateId, OuterBlockKind> blocks12[] = {
      {st.native_inward[1], OuterBlockKind::kAligned},
      {st.native_outward[1], OuterBlockKind::kAligned},
      {st.upper_outward[1], OuterBlockKind::kUnattached}};
  for (const auto& [x, kind] : blocks12) {
    const std::uint32_t u = st.block[x];
    const auto& b = sc.blocks[u];
    CHECK(b.kind == kind);
    CHECK(b.s == W{4});
    CHECK(b.t == W{4});
    CHECK(sc.state_atom[x] == kNoAtom);
    const auto mst = joins_of_block(r, JoinKind::kBlockMst, u);
    REQUIRE(mst.size() == 1);
    CHECK(sc.events[mst[0]].grid ==
          (kind == OuterBlockKind::kAligned ? OuterEventGrid::kGlobal : OuterEventGrid::kSeparate));
    CHECK(sc.events[mst[0]].lower == W{4});
    CHECK(sc.events[mst[0]].upper == W{8});
  }
  REQUIRE(!sc.grid.points.empty());
  CHECK(sc.grid.points[0] == W{4});
  CHECK(sc.grid.points.back() == kInf<W>);  // every support end saturates (D-13)
  CHECK(sssp::paper::outer::check_outer_scales(r.c.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.calendars,
                                               r.os, sc)
            .empty());
}

// ------------------------------------------------------------- several grid components (SC-9)

TEST_CASE("outer scales: several grid components with dyadic doubles (SC-9)",
          "[paper][outer][scales]") {
  using W = double;
  // The path 0 - 1 - 2 - 3 with b = 1: every vertex is a pivot, every edge a corridor and a
  // macro named by its own weight. h = 2 and P = 1 give f = 28, so the supports are
  // [1, 2^84], [2^200, 2^284] and [3 2^283, 3 2^367]: the first is more than a factor two
  // below the second (a separate component), the third starts below 2 * 2^284 (merged).
  const W w1 = 0x1p200;
  const W w2 = 0x3p283;
  const std::vector<Edge<W>> edges{{0, 1, W{1}}, {1, 2, w1}, {2, 3, w2}};
  auto r = to::run_scales(to::make_outer_case<W>(4, edges, 1, to::Names::kExactP1, 1));
  const auto& sc = r.sc;
  REQUIRE(r.c.od.pivots.size() == 4);
  REQUIRE(sc.factor.log2 == 28);
  const auto& pts = sc.grid.points;
  // Component 0: 2^0 .. 2^84 and the sentinel 2^85; component 1: 2^200 .. 2^369 (the first
  // point strictly above 3 2^367).
  REQUIRE(sc.grid.component_begin == std::vector<std::size_t>{0, 86});
  REQUIRE(pts.size() == 86 + 170);
  CHECK(pts[0] == W{1});
  CHECK(pts[85] == 0x1p85);
  CHECK(pts[86] == w1);
  CHECK(pts.back() == 0x1p369);
  for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
    if (i + 1 == 86) {
      CHECK(pts[i] + pts[i] <= pts[i + 1]);  // between components: at least a factor two
    } else {
      CHECK(pts[i + 1] == pts[i] + pts[i]);
    }
  }
  // Every key: lower_index is the floor of its start, upper_index the first point above its end.
  REQUIRE(sc.key_weight.size() == sc.num_name_keys + 3);
  for (std::size_t k = 0; k < sc.key_weight.size(); ++k) {
    const auto& sp = sc.grid.supports[k];
    CHECK(sp.lower == sc.key_weight[k]);
    CHECK(pts[sp.lower_index] <= sp.lower);
    CHECK(sp.lower < pts[sp.lower_index + 1]);
    CHECK(sp.upper < pts[sp.upper_index]);
    CHECK(pts[sp.upper_index - 1] <= sp.upper);
    // No class used by this key crosses a component start.
    CHECK((sp.lower_index < 86) == (sp.upper_index <= 86));
  }
  // The joins: the corridor of weight 3 2^283 falls in the class [2^284, 2^285) of the merged
  // component, the corridor of weight 1 in the class [1, 2).
  const auto& jr = r.os.joins.joins;
  std::size_t corridors = 0;
  for (std::size_t i = 0; i < jr.size(); ++i) {
    const auto& ev = sc.events[i];
    REQUIRE(ev.grid == OuterEventGrid::kGlobal);
    CHECK(ev.lower <= jr[i].weight);
    CHECK(jr[i].weight < ev.upper);
    CHECK(ev.upper == ev.lower + ev.lower);
    if (jr[i].kind != JoinKind::kCorridor) continue;
    ++corridors;
    if (jr[i].weight == W{1}) CHECK(ev.index == 0);
    if (jr[i].weight == w1) CHECK(ev.index == 86);
    if (jr[i].weight == w2) CHECK(ev.index == 86 + 84);
  }
  CHECK(corridors == 3);
  CHECK(sssp::paper::outer::check_outer_scales(r.c.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.calendars,
                                               r.os, sc)
            .empty());
}

// ------------------------------------------------------------------ unattached blocks (SC-10)

TEMPLATE_TEST_CASE("outer scales: one pivot, every block unattached (SC-10, S0-S7)",
                   "[paper][outer][scales]", std::uint64_t, double) {
  using W = TestType;
  // Mode A: b = L_6^64 saturates, one pivot, no interowner edge, no name, no corridor.
  auto r = to::run_scales(to::make_outer_case<W>(6, to::level4_edges<W>(),
                                                 sssp::paper::owner::outer_separator_size(6),
                                                 to::Names::kNameRoutesPaper, 2));
  const auto& sc = r.sc;
  REQUIRE(r.c.od.pivots.size() == 1);
  // F from the exact b = 3^64 (L_6 = 3), never from the saturated word (Q-71).
  CHECK(sc.factor.log2 == paper_outer_scale_factor(r.c.cls.p_log2, 6).log2);
  CHECK(sc.grid.points.empty());
  CHECK(sc.key_weight.empty());
  CHECK(count_kind(r, OuterBlockKind::kAligned) == 0);
  CHECK(count_kind(r, OuterBlockKind::kUnaligned) == 0);
  CHECK(count_kind(r, OuterBlockKind::kUnattached) > 0);
  for (const auto& ev : sc.events) CHECK(ev.grid != OuterEventGrid::kGlobal);
  for (const auto& b : sc.blocks) {
    if (b.kind == OuterBlockKind::kSingleton) continue;
    CHECK(b.t <= b.s);
    CHECK(b.s < b.t + b.t);
  }
  CHECK(r.diag.fallbacks.empty());
  CHECK(sssp::paper::outer::check_outer_scales(r.c.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.calendars,
                                               r.os, sc)
            .empty());
}

TEMPLATE_TEST_CASE("outer scales: rejects inputs of different constructions",
                   "[paper][outer][scales]", std::uint64_t, double) {
  using W = TestType;
  auto r = to::run_scales(
      to::make_outer_case<W>(4, to::level3_edges<W>(W{10}), 2, to::Names::kExactP1, 5));
  auto other = to::make_outer_case<W>(6, to::level4_edges<W>(), 3, to::Names::kExactP1, 5);
  using sssp::paper::outer::build_outer_scales;
  CHECK_THROWS_AS(
      build_outer_scales(other.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.calendars, r.os, r.f),
      std::invalid_argument);
  CHECK_THROWS_AS(
      build_outer_scales(r.c.g, r.c.od, r.c.cls, r.c.aux, other.jr.calendars, r.os, r.f),
      std::invalid_argument);
  auto os = r.os;
  os.joins.theta.pop_back();
  CHECK_THROWS_AS(build_outer_scales(r.c.g, r.c.od, r.c.cls, r.c.aux, r.c.jr.calendars, os, r.f),
                  std::invalid_argument);
  // The merge plan reads the events and the atoms of the scales: scales cut short do not belong
  // to these outer states.
  auto sc = r.sc;
  sc.events.pop_back();
  CHECK_THROWS_AS(sssp::paper::outer::build_outer_merge_plan(r.os, sc), std::invalid_argument);
  CHECK_THROWS_AS(sssp::paper::outer::build_outer_partition(r.os, sc), std::invalid_argument);
}

// ------------------------------------------------- S8-S9: partition hierarchy and substitution

namespace {

/// Alg 8 l.14-15 on a substituted hierarchy: one seed at the native inward copy of s, Evaluate,
/// minimum emitted value per vertex.
template <class W>
std::vector<W> evaluate_from(const to::ScalesCase<W>& r,
                             const sssp::paper::AtomSubstitution<W>& sub, Vertex s) {
  const std::vector<sssp::paper::StateOffer<W>> seed{{r.os.states.native_inward[s], W{0}}};
  sssp::paper::SeededOrderedHierarchy<W> seeded{
      sub.hierarchy, {{s, W{0}}}, sssp::paper::redirect_offers<W>(sub, seed)};
  const auto ev = sssp::paper::evaluate(seeded);
  std::vector<W> dist(r.c.n, kInf<W>);
  for (StateId x = 0; x < sub.hierarchy.num_states(); ++x) {
    const Vertex v = r.os.states.projection[sub.old_state[x]];
    dist[v] = std::min(dist[v], ev.value[x]);
  }
  return dist;
}

/// H4 as produced: every arc between two roots goes forward in `roots`.
template <class W>
bool roots_in_condensation_order(const sssp::paper::outer::OuterPartition<W>& part) {
  const auto& nodes = part.hierarchy.nodes;
  std::vector<std::size_t> rank(nodes.size(), 0);
  for (std::size_t i = 0; i < part.roots.size(); ++i) rank[part.roots[i]] = i;
  const auto root_of = [&](StateId x) {
    sssp::paper::NodeId v = part.leaf_of_state[x];
    while (nodes[v].parent != sssp::paper::kNoNode) v = nodes[v].parent;
    return v;
  };
  for (const auto& arc : part.hierarchy.graph.arcs) {
    if (rank[root_of(arc.tail)] > rank[root_of(arc.head)]) return false;
  }
  return true;
}

}  // namespace

TEMPLATE_TEST_CASE("outer scales: hand example L4, partition hierarchy and substitution (SC-6)",
                   "[paper][outer][scales]", std::uint64_t, double) {
  using W = TestType;
  using sssp::paper::kNoNode;
  using sssp::paper::NodeId;
  auto r =
      to::run_scales(to::make_outer_case<W>(6, to::level4_edges<W>(), 3, to::Names::kExactP1, 9));
  const auto& st = r.os.states;
  sssp::Diagnostics diag;
  const auto part = sssp::paper::outer::build_outer_partition(r.os, r.sc, &diag);
  const auto& ph = part.hierarchy;
  CHECK(diag.fallbacks.empty());
  CHECK(ph.c0 == 2);
  // Parts: the three atoms {4, 5} and every other state on its own.
  const std::size_t parts = st.num_states() - 3;
  REQUIRE(ph.atoms.size() == 3);
  // D loses the two arcs of the edge 4-5 inside each atom.
  CHECK(ph.graph.arcs.size() + 6 == st.graph.arcs.size());
  // Internal nodes, by hand:
  //   * the three blocks {1, 2}: separate class [4, 8) of the edge 1-2 (weight 4): height 4,
  //     two leaf children, one link of length 4;
  //   * the macro 0-3 (weight 100, name and corridor 100, grid 100, 200, ...): the four pivot
  //     copies are the only states with a finite theta (the macro head resets to in(0) or
  //     in(3) and leaves from any copy of 0 or 3; no nonpivot copy lies on a macro cycle), all
  //     at theta = 100, and the corridor join has weight 100: one node of the class [100, 200)
  //     with the four pivot copies as children.
  REQUIRE(ph.nodes.size() == parts + 4);
  std::size_t low = 0;
  std::size_t high = 0;
  for (NodeId v = static_cast<NodeId>(parts); v < ph.nodes.size(); ++v) {
    const auto& node = ph.nodes[v];
    CHECK(node.parent == kNoNode);
    if (node.height == W{4}) {
      ++low;
      CHECK(part.node_upper[v] == W{8});
      CHECK(part.node_grid[v] == OuterEventGrid::kSeparate);
      REQUIRE(node.children.size() == 2);
      REQUIRE(node.links.size() == 1);
      CHECK(node.links[0].length == W{4});
      for (const NodeId c : node.children) {
        const Vertex pv = st.projection[ph.nodes[c].state];
        CHECK((pv == 1 || pv == 2));
      }
    } else {
      ++high;
      CHECK(node.height == W{100});
      CHECK(part.node_upper[v] == W{200});
      CHECK(part.node_grid[v] == OuterEventGrid::kGlobal);
      REQUIRE(node.children.size() == 4);
      std::vector<StateId> kids;
      for (const NodeId c : node.children) kids.push_back(ph.nodes[c].state);
      std::ranges::sort(kids);
      std::vector<StateId> pivots{st.native_inward[0], st.native_outward[0], st.native_inward[3],
                                  st.native_outward[3]};
      std::ranges::sort(pivots);
      CHECK(kids == pivots);
      // Links are arcs strictly below 200: the switches (0) and the macro arcs (100).
      for (const auto& link : node.links) CHECK(link.length <= W{100});
    }
  }
  CHECK(low == 3);
  CHECK(high == 1);
  // Atom leaves are roots here: the blocks {4, 5} have no significant edge and no attachment.
  for (std::uint32_t y = 0; y < 3; ++y) {
    const NodeId leaf = part.leaf_of_state[ph.atoms[y].states[0]];
    CHECK(ph.leaf_atom[leaf] == y);
    CHECK(ph.nodes[leaf].parent == kNoNode);
    CHECK(part.leaf_of_state[ph.atoms[y].states[1]] == leaf);
  }
  CHECK(part.roots.size() == parts - 3 - 3);
  CHECK(roots_in_condensation_order(part));

  // S9 (Lemma 5.15): accepted by substitute_atoms; one Alg 4 call per atom, each the recorded
  // base case on 2 vertices and 1 edge; Def 4.1 holds with c0 <= 2.
  sssp::Diagnostics sub_diag;
  const auto sub = sssp::paper::substitute_atoms(ph, &sub_diag);
  REQUIRE(sub_diag.fallbacks.size() == 3);
  for (const auto& fb : sub_diag.fallbacks) {
    CHECK(fb.reason == "local: k <= 2*L_k^4 (paper base case)");
    CHECK(fb.vertices == 2);
    CHECK(fb.edges == 1);
  }
  CHECK(sub.hierarchy.c0 <= 2);
  const auto report = sssp::paper::validate_ordered_hierarchy(sub.hierarchy, r.c.g);
  INFO(report.summary());
  CHECK(report.ok());
  // Distances by hand on the tree 2 -4- 1 -4- 0 -100- 3 -8- 4 -1- 5; the source 5 lies in atoms.
  CHECK(evaluate_from(r, sub, 0) == std::vector<W>{W{0}, W{4}, W{8}, W{100}, W{108}, W{109}});
  CHECK(evaluate_from(r, sub, 5) == std::vector<W>{W{109}, W{113}, W{117}, W{9}, W{1}, W{0}});
}

TEST_CASE("outer scales: unaligned atoms precede their first outer merge (SC-8, S8-S9)",
          "[paper][outer][scales]") {
  using W = std::uint64_t;
  using sssp::paper::kNoNode;
  using sssp::paper::NodeId;
  const W corridor = (W{1} << 63) + 1;
  const std::vector<Edge<W>> edges{{0, 1, W{4}}, {1, 2, W{4}}, {0, 3, corridor},
                                   {3, 4, W{2}}, {4, 5, W{1}}, {2, 5, corridor + 1}};
  auto r = to::run_scales(to::make_outer_case_named<W>(6, edges, 3, {corridor}, {0, 0}));
  const auto& st = r.os.states;
  const auto part = sssp::paper::outer::build_outer_partition(r.os, r.sc);
  const auto& ph = part.hierarchy;
  // l.2093: an unaligned block (here a whole atom, cutoff 2) is one node before its first
  // global event: the parent of its leaf is a node of the class [2^63, infinity).
  for (const StateId x : {st.native_inward[4], st.native_outward[4]}) {
    const NodeId leaf = part.leaf_of_state[x];
    REQUIRE(ph.leaf_atom[leaf] != kNoAtom);
    const NodeId par = ph.nodes[leaf].parent;
    REQUIRE(par != kNoNode);
    CHECK(ph.nodes[par].height == W{1} << 63);
    CHECK(part.node_upper[par] == kInf<W>);
    CHECK(ph.atoms[ph.leaf_atom[leaf]].cutoff <= ph.nodes[par].height);
  }
  // The aligned blocks {1, 2} merge in the global class [4, 8) and join the core later.
  for (const StateId x : {st.native_inward[1], st.native_outward[1]}) {
    const NodeId par = ph.nodes[part.leaf_of_state[x]].parent;
    REQUIRE(par != kNoNode);
    CHECK(ph.nodes[par].height == W{4});
    CHECK(part.node_upper[par] == W{8});
    CHECK(part.node_grid[par] == OuterEventGrid::kGlobal);
    REQUIRE(ph.nodes[par].parent != kNoNode);
    CHECK(ph.nodes[ph.nodes[par].parent].height == W{1} << 63);
  }
  CHECK(roots_in_condensation_order(part));
  // Lemma 5.15 accepts the partition hierarchy (parts, outside arcs, cutoffs).
  sssp::Diagnostics sub_diag;
  const auto sub = sssp::paper::substitute_atoms(ph, &sub_diag);
  CHECK(sub_diag.fallbacks.size() == 3);
  CHECK(sub.hierarchy.c0 <= 2);
}

TEMPLATE_TEST_CASE("outer scales: one pivot, several roots in condensation order (SC-10)",
                   "[paper][outer][scales]", std::uint64_t, double) {
  using W = TestType;
  auto r = to::run_scales(to::make_outer_case<W>(6, to::level4_edges<W>(),
                                                 sssp::paper::owner::outer_separator_size(6),
                                                 to::Names::kNameRoutesPaper, 2));
  REQUIRE(r.c.od.pivots == std::vector<Vertex>{0});
  const auto part = sssp::paper::outer::build_outer_partition(r.os, r.sc);
  const auto& ph = part.hierarchy;
  // Without macro edges every SCC of the state graph lies in one block (the macro-free block
  // quotient is acyclic, Lemma 8.1): separate grids only, and more than one root.
  for (std::size_t v = 0; v < ph.nodes.size(); ++v) {
    if (!ph.nodes[v].children.empty()) CHECK(part.node_grid[v] == OuterEventGrid::kSeparate);
  }
  CHECK(part.roots.size() > 1);
  CHECK(roots_in_condensation_order(part));
  sssp::Diagnostics sub_diag;
  const auto sub = sssp::paper::substitute_atoms(ph, &sub_diag);
  const auto report = sssp::paper::validate_ordered_hierarchy(sub.hierarchy, r.c.g);
  INFO(report.summary());
  CHECK(report.ok());
  CHECK(sub.hierarchy.c0 <= 2);
  // The single-pivot run is the paper's construction (Q-72): distances from the pivot's side.
  CHECK(evaluate_from(r, sub, 0) == std::vector<W>{W{0}, W{4}, W{8}, W{100}, W{108}, W{109}});
  CHECK(evaluate_from(r, sub, 4) == std::vector<W>{W{108}, W{112}, W{116}, W{8}, W{0}, W{1}});
}
