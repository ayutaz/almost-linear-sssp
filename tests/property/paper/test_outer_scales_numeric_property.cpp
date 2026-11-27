// Property tests for the numeric parts of Sec 8.2 and Lemma 8.4 (arXiv 2611.34911v1, txt
// l.2076-2142) on seeded random inputs, against a naive reference written in this file.
//
// The reference works on integers in unsigned __int128, never saturates, and follows the
// definitions literally: F^2 s >= a by one multiplication, supports [min(a, aligned s), F^3 a],
// components by merging any two intervals that overlap or are closer than a factor two until
// nothing changes, grids by listing lo, 2 lo, 4 lo, ... through the first point above the
// component end, floors and classes by linear scans.
//
// Two implementations are compared with it on the same instances:
//   * double with dyadic inputs x * 2^-20 (Q-55 option A), x with at most 53 significant bits.
//     Nothing saturates, so every output must be the scaled reference output exactly.
//   * uint64. F^3 a rarely fits, so this is the saturated regime of D-13: the output must be
//     the reference with every value >= 2^64 - 1 replaced by infinity() and the grid cut after
//     the first such point. This is the claim that saturation changes no decision.
//
// The inequalities of Lemma 8.4 item 1 are checked on the outputs:
//   unaligned (F^2 s < a): every w <= 4 K^2 h s has w < a / 2;
//   aligned with s <= F a: 4 K^2 h s <= F^3 a, s/2 < t <= s, and every significant weight
//   t <= w <= 4 K^2 h s falls in a class [g, g') with g' <= 2g of its name's index range.
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <vector>

#include "paper/outer/scales_numeric.hpp"
#include "sssp/graph.hpp"

using sssp::paper::outer::build_outer_global_grid;
using sssp::paper::outer::first_edge_not_below_half;
using sssp::paper::outer::kNoGridClass;
using sssp::paper::outer::outer_aligned_cutoff;
using sssp::paper::outer::outer_aligned_scale_in_range;
using sssp::paper::outer::outer_block_aligned;
using sssp::paper::outer::outer_block_scale;
using sssp::paper::outer::outer_cutoff_valid;
using sssp::paper::outer::outer_grid_class;
using sssp::paper::outer::outer_scale_factor;
using sssp::paper::outer::outer_separate_cutoff;
using sssp::paper::outer::outer_window_edge_bound;
using sssp::paper::outer::OuterAlignedScale;
using sssp::paper::outer::OuterGlobalGrid;
using sssp::paper::outer::OuterScaleFactor;

namespace {

using U128 = unsigned __int128;
constexpr std::uint64_t kU64Max = std::numeric_limits<std::uint64_t>::max();

struct Params {
  unsigned p;       ///< log2 P
  std::uint64_t h;  ///< hop bound
  unsigned k;       ///< log2 K
};

/// The definition of F (l.2080) in 128 bits; the parameters keep the product far below 2^127.
unsigned ref_log2_f(const Params& q) {
  const U128 big_k = U128{1} << q.k;
  const U128 inner = 8 * big_k * big_k * q.h;
  const U128 product = 2 * (U128{1} << q.p) * (2 * U128{q.h} + 1) * inner * inner;
  unsigned e = 0;
  while ((U128{1} << e) <= product) ++e;
  return e;
}

/// 4 K^2 h s (Lemma 8.4 item 1).
U128 ref_window_bound(const Params& q, std::uint64_t s) {
  return 4 * (U128{1} << (2 * q.k)) * q.h * s;
}

std::uint64_t random_bits(std::mt19937_64& rng, unsigned max_bits) {
  const unsigned bits = 1 + static_cast<unsigned>(rng() % max_bits);
  const std::uint64_t x = bits >= 64 ? rng() : rng() & ((std::uint64_t{1} << bits) - 1);
  return x == 0 ? 1 : x;
}

std::uint64_t random_in(std::mt19937_64& rng, std::uint64_t lo, std::uint64_t hi) {
  return lo + rng() % (hi - lo + 1);
}

struct RefSupport {
  U128 lower;
  U128 upper;
};

struct RefGrid {
  std::vector<RefSupport> supports;
  std::vector<U128> points;
  std::vector<std::size_t> component_begin;
};

RefGrid ref_grid(unsigned log2_f, const std::vector<std::uint64_t>& names,
                 const std::vector<OuterAlignedScale<std::uint64_t>>& blocks) {
  RefGrid g;
  for (const std::uint64_t a : names) g.supports.push_back({a, U128{a} << (3 * log2_f)});
  for (const auto& b : blocks) {
    g.supports[b.name].lower = std::min<U128>(g.supports[b.name].lower, b.s);
  }
  // Components: merge any two intervals that overlap or leave a gap below a factor two.
  std::vector<RefSupport> comps = g.supports;
  for (bool changed = true; changed;) {
    changed = false;
    for (std::size_t i = 0; i < comps.size() && !changed; ++i) {
      for (std::size_t j = i + 1; j < comps.size() && !changed; ++j) {
        if (comps[j].lower < 2 * comps[i].upper && comps[i].lower < 2 * comps[j].upper) {
          comps[i] = {std::min(comps[i].lower, comps[j].lower),
                      std::max(comps[i].upper, comps[j].upper)};
          comps.erase(comps.begin() + static_cast<std::ptrdiff_t>(j));
          changed = true;
        }
      }
    }
  }
  std::sort(comps.begin(), comps.end(),
            [](const RefSupport& x, const RefSupport& y) { return x.lower < y.lower; });
  for (const RefSupport& c : comps) {
    U128 x = c.lower;
    if (g.points.empty() || g.points.back() != x) {
      g.component_begin.push_back(g.points.size());
      g.points.push_back(x);
    } else {
      g.component_begin.push_back(g.points.size() - 1);
    }
    while (x <= c.upper) {  // through the first point strictly above the end
      x *= 2;
      g.points.push_back(x);
    }
  }
  return g;
}

std::size_t ref_floor(const std::vector<U128>& points, U128 x) {
  std::size_t best = points.size();
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (points[i] <= x) best = i;
  }
  return best;
}

std::size_t ref_first_above(const std::vector<U128>& points, U128 x) {
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (points[i] > x) return i;
  }
  return points.size();
}

constexpr Params kSmallParams[] = {{0, 2, 1}, {1, 2, 1}, {0, 3, 1}, {0, 2, 0}, {1, 3, 0}};

double dy(std::uint64_t x) {
  return static_cast<double>(x) * 0x1p-20;
}

/// Exact for values with at most 53 significant bits.
double dy128(U128 x) {
  return static_cast<double>(x) * 0x1p-20;
}

std::uint64_t clamp64(U128 x) {
  return x >= kU64Max ? kU64Max : static_cast<std::uint64_t>(x);
}

}  // namespace

TEST_CASE("outer scales numeric property: alignment equals F^2 s >= a",
          "[paper][outer][scales_numeric][property]") {
  std::mt19937_64 rng(84001);
  std::size_t aligned = 0;
  std::size_t unaligned = 0;
  for (int it = 0; it < 20000; ++it) {
    const Params q{static_cast<unsigned>(rng() % 6), random_in(rng, 2, 40),
                   static_cast<unsigned>(rng() % 4)};
    const unsigned e = ref_log2_f(q);
    REQUIRE(e <= 60);
    const auto f = outer_scale_factor(q.p, q.h, q.k);
    REQUIRE(f.log2 == e);
    REQUIRE(f.value == std::uint64_t{1} << e);
    std::uint64_t s = random_bits(rng, 20);
    std::uint64_t a = random_bits(rng, 63);
    if (it % 4 == 0) {  // near the boundary F^2 s = a when it fits
      const U128 edge = U128{s} << (2 * e);
      if (edge < kU64Max - 2) a = static_cast<std::uint64_t>(edge) + (rng() % 3) - 1;
    }
    if (a == kU64Max) a = kU64Max - 1;
    const bool want = (U128{s} << (2 * e)) >= a;
    CAPTURE(q.p, q.h, q.k, s, a);
    REQUIRE(outer_block_aligned<std::uint64_t>(f, s, a) == want);
    (want ? aligned : unaligned) += 1;

    // Lemma 8.4 item 1, unaligned: every window edge w <= 4 K^2 h s is below a / 2.
    const U128 bound = ref_window_bound(q, s);
    const std::uint64_t got_bound = outer_window_edge_bound<std::uint64_t>(s, q.h, q.k);
    REQUIRE(got_bound == (bound >= kU64Max ? kU64Max : static_cast<std::uint64_t>(bound)));
    if (!want) {
      REQUIRE(2 * bound < a);
      std::vector<std::uint64_t> w;
      for (int i = 0; i < 6; ++i) w.push_back(random_in(rng, 1, got_bound));
      w.push_back(got_bound);
      REQUIRE(first_edge_not_below_half<std::uint64_t>(std::span<const std::uint64_t>(w), a) ==
              w.size());
    }
  }
  CHECK(aligned > 1000);
  CHECK(unaligned > 1000);
}

TEST_CASE("outer scales numeric property: half-name test and s_U against linear scans",
          "[paper][outer][scales_numeric][property]") {
  std::mt19937_64 rng(84002);
  for (int it = 0; it < 5000; ++it) {
    const std::size_t n = 1 + rng() % 12;
    std::vector<std::uint64_t> w;
    for (std::size_t i = 0; i < n; ++i) w.push_back(random_bits(rng, it % 2 == 0 ? 12 : 64));
    for (auto& x : w) x = std::min(x, kU64Max - 1);
    std::uint64_t a = random_bits(rng, it % 2 == 0 ? 14 : 64);
    a = std::min(a, kU64Max - 1);
    std::size_t want = n;
    for (std::size_t i = n; i-- > 0;) {
      if (2 * U128{w[i]} >= a) want = i;
    }
    REQUIRE(first_edge_not_below_half<std::uint64_t>(std::span<const std::uint64_t>(w), a) == want);
    std::uint64_t least = w[0];
    for (const std::uint64_t x : w) least = x < least ? x : least;
    REQUIRE(outer_block_scale<std::uint64_t>(std::span<const std::uint64_t>(w)) == least);
  }
}

TEST_CASE("outer scales numeric property: supports, grid, cutoffs and classes",
          "[paper][outer][scales_numeric][property]") {
  std::mt19937_64 rng(84003);
  std::size_t multi_component = 0;
  std::size_t shared_sentinel = 0;
  std::size_t classified = 0;
  std::size_t saturated_grids = 0;
  for (int it = 0; it < 1500; ++it) {
    const Params q = kSmallParams[rng() % std::size(kSmallParams)];
    const unsigned e = ref_log2_f(q);
    REQUIRE(e <= 17);
    const auto f = outer_scale_factor(q.p, q.h, q.k);
    REQUIRE(f.log2 == e);

    // Names with at most 11 significant bits, below 2^62, in clusters far apart, so that
    // separate components occur (a gap of at least 2 F^3 between clusters).
    const std::size_t nk = 1 + rng() % 5;
    std::vector<std::uint64_t> names;
    for (std::size_t i = 0; i < nk; ++i) {
      const unsigned shift = rng() % 2 == 0 ? 0 : static_cast<unsigned>(rng() % 51);
      names.push_back(random_bits(rng, 11) << shift);
    }
    if (it % 5 == 0 && 3 * e + 2 < 60) {
      // A gap of exactly two: the second name is twice the end of the first support.
      names = {std::uint64_t{1} << (rng() % 4)};
      names.push_back(names[0] << (3 * e + 1));
    }
    const std::size_t nnames = names.size();
    // Aligned blocks with a / F^2 <= s <= F a (Lemma 8.4, l.2118); s has at most 16 significant
    // bits and is below 2^62.
    std::vector<OuterAlignedScale<std::uint64_t>> blocks;
    const std::size_t tries = rng() % 9;
    for (std::size_t i = 0; i < tries; ++i) {
      const auto name = static_cast<std::uint32_t>(rng() % nnames);
      const std::uint64_t a = names[name];
      const std::uint64_t s = random_bits(rng, 16) << (rng() % 47);
      const bool lower_ok = (U128{s} << (2 * e)) >= a;
      const bool upper_ok = U128{s} <= (U128{a} << e);
      REQUIRE(outer_block_aligned<std::uint64_t>(f, s, a) == lower_ok);
      REQUIRE(outer_aligned_scale_in_range<std::uint64_t>(f, s, a) == (lower_ok && upper_ok));
      REQUIRE(outer_block_aligned<double>(f, dy(s), dy(a)) == lower_ok);
      REQUIRE(outer_aligned_scale_in_range<double>(f, dy(s), dy(a)) == (lower_ok && upper_ok));
      if (lower_ok && upper_ok) blocks.push_back({s, name});
    }
    std::vector<double> dnames;
    for (const std::uint64_t a : names) dnames.push_back(dy(a));
    std::vector<OuterAlignedScale<double>> dblocks;
    for (const auto& b : blocks) dblocks.push_back({dy(b.s), b.name});

    const RefGrid want = ref_grid(e, names, blocks);
    const auto got = build_outer_global_grid<std::uint64_t>(
        f, std::span<const std::uint64_t>(names),
        std::span<const OuterAlignedScale<std::uint64_t>>(blocks));
    const auto dgot = build_outer_global_grid<double>(
        f, std::span<const double>(dnames), std::span<const OuterAlignedScale<double>>(dblocks));

    multi_component += want.component_begin.size() > 1 ? 1 : 0;
    for (std::size_t c = 1; c < want.component_begin.size(); ++c) {
      // A shared sentinel: the previous grid's last doubling is this component's start.
      const std::size_t b = want.component_begin[c];
      if (want.points[b] == 2 * want.points[b - 1]) ++shared_sentinel;
    }

    // double: the reference exactly.
    REQUIRE(dgot.points.size() == want.points.size());
    for (std::size_t i = 0; i < want.points.size(); ++i) {
      REQUIRE(dgot.points[i] == dy128(want.points[i]));
    }
    REQUIRE(dgot.component_begin == want.component_begin);
    // uint64: the reference cut after its first point that does not fit (D-13).
    std::size_t cut = want.points.size();
    for (std::size_t i = want.points.size(); i-- > 0;) {
      if (want.points[i] >= kU64Max) cut = i + 1;
    }
    saturated_grids += want.points.back() >= kU64Max ? 1 : 0;
    REQUIRE(got.points.size() == cut);
    for (std::size_t i = 0; i < cut; ++i) REQUIRE(got.points[i] == clamp64(want.points[i]));
    REQUIRE(got.component_begin == want.component_begin);

    REQUIRE(got.supports.size() == nnames);
    REQUIRE(dgot.supports.size() == nnames);
    for (std::size_t c = 0; c < nnames; ++c) {
      const std::size_t lower_index = ref_floor(want.points, want.supports[c].lower);
      const std::size_t upper_index = ref_first_above(want.points, want.supports[c].upper);
      REQUIRE(upper_index < want.points.size());
      REQUIRE(dgot.supports[c].lower == dy128(want.supports[c].lower));
      REQUIRE(dgot.supports[c].upper == dy128(want.supports[c].upper));
      REQUIRE(dgot.supports[c].lower_index == lower_index);
      REQUIRE(dgot.supports[c].upper_index == upper_index);
      REQUIRE(got.supports[c].lower == want.supports[c].lower);
      REQUIRE(got.supports[c].upper == clamp64(want.supports[c].upper));
      REQUIRE(got.supports[c].lower_index == lower_index);
      REQUIRE(got.supports[c].upper_index == std::min(upper_index, cut - 1));
      // l.2081, l.2086: the support holds its name.
      REQUIRE(got.supports[c].lower <= names[c]);
      REQUIRE(names[c] <= got.supports[c].upper);
    }

    for (const auto& b : blocks) {
      const auto& sup = got.supports[b.name];
      const std::uint64_t a = names[b.name];
      REQUIRE(sup.lower <= b.s);
      // l.2090-2092: t_U is the grid floor of s_U and s_U / 2 < t_U <= s_U.
      const U128 want_t = want.points[ref_floor(want.points, b.s)];
      const std::uint64_t t = outer_aligned_cutoff<std::uint64_t>(got, b.name, b.s);
      REQUIRE(U128{t} == want_t);
      REQUIRE(outer_aligned_cutoff<double>(dgot, b.name, dy(b.s)) == dy128(want_t));
      REQUIRE(t <= b.s);
      REQUIRE(U128{b.s} < 2 * U128{t});
      REQUIRE(outer_cutoff_valid<std::uint64_t>(b.s, t));
      // l.2119: the floor of the support start is at or below t_U.
      REQUIRE(got.points[sup.lower_index] <= t);
      // Significant weights t <= w <= 4 K^2 h s are inside the support (in exact arithmetic)
      // and inside the support's index range.
      const U128 bound = ref_window_bound(q, b.s);
      REQUIRE(bound <= (U128{a} << (3 * e)));
      const std::uint64_t top = bound >= kU64Max ? kU64Max - 1 : static_cast<std::uint64_t>(bound);
      for (int i = 0; i < 8; ++i) {
        const std::uint64_t w = i == 0 ? t : (i == 1 ? top : random_in(rng, t, top));
        const std::size_t cls = outer_grid_class<std::uint64_t>(got, b.name, w);
        REQUIRE(cls != kNoGridClass);
        REQUIRE(cls == ref_floor(want.points, w));
        REQUIRE(cls >= sup.lower_index);
        REQUIRE(cls < sup.upper_index);
        REQUIRE(got.points[cls] <= w);
        REQUIRE(w < got.points[cls + 1]);
        REQUIRE(U128{got.points[cls + 1]} <= 2 * U128{got.points[cls]});
        ++classified;
      }
      // Below the designated floor a weight is local (l.2122): no class.
      if (got.points[sup.lower_index] > 1) {
        REQUIRE(outer_grid_class<std::uint64_t>(got, b.name, got.points[sup.lower_index] - 1) ==
                kNoGridClass);
      }
    }
    // Every name is classified on its own index range, identically in both types.
    for (std::size_t c = 0; c < nnames; ++c) {
      const auto name = static_cast<std::uint32_t>(c);
      const std::size_t cls = ref_floor(want.points, names[c]);
      REQUIRE(outer_grid_class<std::uint64_t>(got, name, names[c]) == cls);
      REQUIRE(outer_grid_class<double>(dgot, name, dnames[c]) == cls);
    }
  }
  CHECK(multi_component > 50);
  CHECK(shared_sentinel > 50);
  CHECK(classified > 5000);
  CHECK(saturated_grids > 50);
}

TEST_CASE("outer scales numeric property: double decisions equal the integer decisions",
          "[paper][outer][scales_numeric][property]") {
  std::mt19937_64 rng(84004);
  for (int it = 0; it < 20000; ++it) {
    const Params q = kSmallParams[rng() % std::size(kSmallParams)];
    const auto f = outer_scale_factor(q.p, q.h, q.k);
    const std::uint64_t s = random_bits(rng, 18);
    std::uint64_t a = random_bits(rng, 52);
    if (it % 3 == 0) {
      const U128 edge = U128{s} << (2 * f.log2);
      if (edge < (U128{1} << 52)) a = static_cast<std::uint64_t>(edge) + (rng() % 3) - 1;
    }
    const std::uint64_t t = random_in(rng, 1, 2 * s);
    CAPTURE(s, a, t, f.log2);
    REQUIRE(outer_block_aligned<double>(f, dy(s), dy(a)) ==
            outer_block_aligned<std::uint64_t>(f, s, a));
    REQUIRE(outer_aligned_scale_in_range<double>(f, dy(s), dy(a)) ==
            outer_aligned_scale_in_range<std::uint64_t>(f, s, a));
    REQUIRE(outer_cutoff_valid<double>(dy(s), dy(t)) == outer_cutoff_valid<std::uint64_t>(s, t));
    REQUIRE(outer_cutoff_valid<std::uint64_t>(s, t) == (t <= s && s < 2 * t));
    const std::vector<std::uint64_t> w{random_bits(rng, 50), random_bits(rng, 50)};
    const std::vector<double> dw{dy(w[0]), dy(w[1])};
    REQUIRE(first_edge_not_below_half<double>(std::span<const double>(dw), dy(a)) ==
            first_edge_not_below_half<std::uint64_t>(std::span<const std::uint64_t>(w), a));
    // 4 K^2 h s stays below 2^53 here, so the integer multiple is exact in double as well.
    REQUIRE(outer_window_edge_bound<double>(dy(s), q.h, q.k) ==
            dy(outer_window_edge_bound<std::uint64_t>(s, q.h, q.k)));
  }
}

TEST_CASE("outer scales numeric property: separate-grid cutoff is the floor on a doubling table",
          "[paper][outer][scales_numeric][property]") {
  std::mt19937_64 rng(84005);
  for (int it = 0; it < 5000; ++it) {
    const std::uint64_t origin = random_bits(rng, 20);
    const std::size_t len = 1 + rng() % 30;
    std::vector<std::uint64_t> table;
    for (std::size_t i = 0; i < len; ++i) table.push_back(origin << i);
    const std::uint64_t s = random_in(rng, origin, 2 * table.back() - 1);
    const std::uint64_t t =
        outer_separate_cutoff<std::uint64_t>(std::span<const std::uint64_t>(table), s);
    std::uint64_t want = 0;
    for (const std::uint64_t x : table) {
      if (x <= s) want = x;
    }
    REQUIRE(t == want);
    REQUIRE(t <= s);
    REQUIRE(s < 2 * t);  // s/2 < t_U <= s_U (Q-36)
  }
}

TEST_CASE("outer scales numeric property: saturated supports keep the grid structure (D-13)",
          "[paper][outer][scales_numeric][property]") {
  std::mt19937_64 rng(84006);
  // F = 2^64 and F = 2^96 (thresholds.md Sec 5.2) and a small F with large names.
  const OuterScaleFactor factors[] = {outer_scale_factor(26, std::uint64_t{21}, 4),
                                      outer_scale_factor(20, std::uint64_t{160001}, 4),
                                      outer_scale_factor(0, std::uint64_t{2}, 1)};
  std::size_t saturated = 0;
  for (int it = 0; it < 1500; ++it) {
    const OuterScaleFactor& f = factors[rng() % 3];
    const std::size_t nk = 1 + rng() % 6;
    std::vector<std::uint64_t> names;
    for (std::size_t i = 0; i < nk; ++i) {
      names.push_back(std::min(random_bits(rng, 64), kU64Max - 1));
    }
    std::vector<OuterAlignedScale<std::uint64_t>> blocks;
    const std::size_t nblocks = rng() % 5;
    for (std::size_t i = 0; i < nblocks; ++i) {
      const auto name = static_cast<std::uint32_t>(rng() % nk);
      const std::uint64_t s = std::min(random_bits(rng, 64), kU64Max - 1);
      const U128 a = names[name];
      // Keep only blocks satisfying a / F^2 <= s <= F a, decided in 128 bits by the exponent.
      const bool lower_ok = 2 * f.log2 >= 64 || (U128{s} << (2 * f.log2)) >= a;
      const bool upper_ok = f.log2 >= 64 || U128{s} <= (a << f.log2);
      REQUIRE(outer_block_aligned<std::uint64_t>(f, s, names[name]) == lower_ok);
      REQUIRE(outer_aligned_scale_in_range<std::uint64_t>(f, s, names[name]) ==
              (lower_ok && upper_ok));
      if (lower_ok && upper_ok) blocks.push_back({s, name});
    }
    const auto g = build_outer_global_grid<std::uint64_t>(
        f, std::span<const std::uint64_t>(names),
        std::span<const OuterAlignedScale<std::uint64_t>>(blocks));

    REQUIRE(!g.points.empty());
    for (std::size_t i = 0; i + 1 < g.points.size(); ++i) REQUIRE(g.points[i] < g.points[i + 1]);
    REQUIRE(!g.component_begin.empty());
    REQUIRE(g.component_begin.front() == 0);
    for (std::size_t c = 0; c < g.component_begin.size(); ++c) {
      const std::size_t begin = g.component_begin[c];
      const std::size_t end =
          c + 1 < g.component_begin.size() ? g.component_begin[c + 1] : g.points.size() - 1;
      // Inside a component every step is one doubling; a step that does not fit is the
      // saturated sentinel. The step into the next component's start is a doubling when the
      // sentinel is shared and otherwise leaves this component's sentinel.
      for (std::size_t i = begin; i < end; ++i) {
        const std::uint64_t doubled = clamp64(2 * U128{g.points[i]});
        if (i + 1 == end && c + 1 < g.component_begin.size() && i > begin) {
          REQUIRE(clamp64(2 * U128{g.points[i - 1]}) == g.points[i]);
        } else {
          REQUIRE(g.points[i + 1] == doubled);
        }
      }
    }
    // l.2088: a later component starts at least a factor two above everything before it.
    for (std::size_t c = 0; c < nk; ++c) {
      for (const std::size_t begin : g.component_begin) {
        if (g.points[begin] > g.supports[c].lower) {
          REQUIRE(U128{g.points[begin]} >= 2 * U128{g.supports[c].upper});
        }
      }
    }
    for (std::size_t c = 0; c < nk; ++c) {
      const auto& sup = g.supports[c];
      const U128 exact_upper = 3 * f.log2 >= 128 ? ~U128{0} : U128{names[c]} << (3 * f.log2);
      const bool sat = 3 * f.log2 >= 64 || exact_upper >= kU64Max;
      REQUIRE(sup.upper == (sat ? kU64Max : static_cast<std::uint64_t>(exact_upper)));
      saturated += sat ? 1 : 0;
      REQUIRE(sup.lower <= names[c]);
      REQUIRE(sup.lower_index < sup.upper_index);
      REQUIRE(sup.upper_index < g.points.size());
      REQUIRE(g.points[sup.lower_index] <= sup.lower);
      REQUIRE(sup.lower < g.points[sup.lower_index + 1]);
      // First point strictly above the end, or the saturated sentinel when the end saturates.
      if (sat) {
        REQUIRE(g.points[sup.upper_index] == kU64Max);
      } else {
        REQUIRE(g.points[sup.upper_index] > sup.upper);
        REQUIRE(g.points[sup.upper_index - 1] <= sup.upper);
      }
      const std::size_t cls =
          outer_grid_class<std::uint64_t>(g, static_cast<std::uint32_t>(c), names[c]);
      REQUIRE(cls != kNoGridClass);
      REQUIRE(g.points[cls] <= names[c]);
      REQUIRE(names[c] < g.points[cls + 1]);
    }
    for (const auto& b : blocks) {
      const std::uint64_t t = outer_aligned_cutoff<std::uint64_t>(g, b.name, b.s);
      REQUIRE(t <= b.s);
      REQUIRE(U128{b.s} < 2 * U128{t});
      REQUIRE(g.supports[b.name].lower <= b.s);
    }
  }
  CHECK(saturated > 1000);
}
