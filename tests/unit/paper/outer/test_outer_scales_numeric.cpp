// Sec 8.2 (arXiv 2611.34911v1, txt l.2076-2104) and Lemma 8.4 (l.2105-2142): the numeric parts
// that do not depend on the outer state types. Hand examples, boundaries, saturation (D-13) and
// rejected inputs.
//
// F is compared with a multi-precision reference written in this file (base 2^32 limbs): the
// product 2P(2h+1)(8K^2h)^2 is formed in full and the least power of two strictly above it is
// found by doubling. Expected exponents of the paper-size and Mode B rows are the ones listed in
// docs/research/thresholds.md Sec 5 (computed independently of the implementation).
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/outer/scales_numeric.hpp"
#include "sssp/graph.hpp"

using sssp::WeightTraits;
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
using sssp::paper::outer::OuterScaleFactor;
using sssp::paper::outer::scale_pow2;

namespace {

constexpr std::uint64_t kU64Max = std::numeric_limits<std::uint64_t>::max();
constexpr double kDInf = std::numeric_limits<double>::infinity();

// ---- multi-precision reference (little-endian base 2^32, no leading zero limbs) ----
using Big = std::vector<std::uint32_t>;

void trim(Big& a) {
  while (!a.empty() && a.back() == 0) a.pop_back();
}

Big big_from(std::uint64_t x) {
  Big a{static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(x >> 32)};
  trim(a);
  return a;
}

Big big_add(const Big& a, const Big& b) {
  Big r;
  std::uint64_t carry = 0;
  for (std::size_t i = 0; i < a.size() || i < b.size() || carry != 0; ++i) {
    std::uint64_t s = carry;
    if (i < a.size()) s += a[i];
    if (i < b.size()) s += b[i];
    r.push_back(static_cast<std::uint32_t>(s));
    carry = s >> 32;
  }
  return r;
}

Big big_mul(const Big& a, const Big& b) {
  Big r(a.size() + b.size() + 1, 0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    std::uint64_t carry = 0;
    for (std::size_t j = 0; j < b.size() || carry != 0; ++j) {
      std::uint64_t cur = r[i + j] + carry;
      if (j < b.size()) cur += static_cast<std::uint64_t>(a[i]) * b[j];
      r[i + j] = static_cast<std::uint32_t>(cur);
      carry = cur >> 32;
    }
  }
  trim(r);
  return r;
}

bool big_less(const Big& a, const Big& b) {
  if (a.size() != b.size()) return a.size() < b.size();
  for (std::size_t i = a.size(); i-- > 0;) {
    if (a[i] != b[i]) return a[i] < b[i];
  }
  return false;
}

Big big_pow(const Big& base, unsigned e) {
  Big r = big_from(1);
  for (unsigned i = 0; i < e; ++i) r = big_mul(r, base);
  return r;
}

/// The definition, l.2080: the least e with 2^e > 2 P (2h+1) (8 K^2 h)^2, P = 2^p, K = 2^k.
std::uint64_t reference_log2_f(unsigned p, const Big& h, unsigned k) {
  const Big two = big_from(2);
  const Big big_p = big_pow(two, p);
  const Big big_k = big_pow(two, k);
  const Big two_h_plus_1 = big_add(big_mul(two, h), big_from(1));
  const Big inner = big_mul(big_from(8), big_mul(big_mul(big_k, big_k), h));  // 8 K^2 h
  const Big product = big_mul(big_mul(two, big_p), big_mul(two_h_plus_1, big_mul(inner, inner)));
  Big pow = big_from(1);
  std::uint64_t e = 0;
  while (!big_less(product, pow)) {  // until pow > product
    pow = big_add(pow, pow);
    ++e;
  }
  return e;
}

template <class W>
std::span<const W> sp(const std::vector<W>& v) {
  return std::span<const W>(v);
}

/// F = 2^16: P = 1, h = 2, K = 2 gives 2 * 5 * (8 * 4 * 2)^2 = 40960.
OuterScaleFactor small_f() {
  return outer_scale_factor(0, std::uint64_t{2}, 1);
}

/// An F of the paper's magnitude, 2^1808 (thresholds.md Sec 5.1), reached with a one-word h:
/// 7 + 1780 + 16 + bitlen(20) = 1808.
OuterScaleFactor paper_f() {
  return outer_scale_factor(1780, std::uint64_t{2}, 4);
}

}  // namespace

TEST_CASE("outer scales numeric: F on hand examples", "[paper][outer][scales_numeric]") {
  // 2 * 1 * 5 * (8 * 4 * 2)^2 = 40960 < 65536 = 2^16.
  const auto f = outer_scale_factor(0, std::uint64_t{2}, 1);
  CHECK(f.log2 == 16);
  CHECK(f.value == 65536);
  CHECK(f.value_exact());
  // 2 * 5 * (8 * 256 * 2)^2 = 167772160, between 2^27 and 2^28.
  const auto g = outer_scale_factor(0, std::uint64_t{2}, 4);
  CHECK(g.log2 == 28);
  CHECK(g.value == (std::uint64_t{1} << 28));
  // 2 * 2 * 7 * (8 * 1 * 3)^2 = 16128, between 2^13 and 2^14.
  const auto e = outer_scale_factor(1, std::uint64_t{3}, 0);
  CHECK(e.log2 == 14);
  CHECK(e.value == 16384);
}

TEST_CASE("outer scales numeric: the F formula gives 2^1808 at the paper's size",
          "[paper][outer][scales_numeric]") {
  // h = 589^64 + 1 does not fit a word, so the implementation cannot be called with it. The
  // reference evaluates the definition in full, and the closed form the implementation uses,
  // 7 + p + 4k + bitlen(h^2 (2h + 1)), is evaluated on the same multi-precision h.
  const Big h = big_add(big_pow(big_from(589), 64), big_from(1));
  CHECK(reference_log2_f(17, h, 4) == 1808);  // thresholds.md Sec 5.1
  const Big m = big_mul(big_mul(h, h), big_add(big_add(h, h), big_from(1)));
  Big pow = big_from(1);
  std::uint64_t bits = 0;
  while (!big_less(m, pow)) {
    pow = big_add(pow, pow);
    ++bits;
  }
  CHECK(7 + 17 + 4 * 4 + bits == 1808);
  // A one-word h reaching the same exponent: the value saturates (D-13).
  const auto f = paper_f();
  CHECK(f.log2 == 1808);
  CHECK(f.log2 == reference_log2_f(1780, big_from(2), 4));
  CHECK(f.value == kU64Max);
  CHECK_FALSE(f.value_exact());
}

TEST_CASE("outer scales numeric: F with small separators and the paper's P",
          "[paper][outer][scales_numeric]") {
  // (h = b + 1, log2 F) at p = 17, K = 16: docs/research/04_outer_scales_pipeline.md Sec 1.4.
  const std::pair<std::uint64_t, std::uint64_t> rows[] = {{2, 45}, {3, 46},  {4, 48},  {5, 49},
                                                          {9, 51}, {17, 54}, {34, 57}, {65, 60}};
  for (const auto& [h, log2_f] : rows) {
    CHECK(outer_scale_factor(17, h, 4).log2 == log2_f);
    CHECK(reference_log2_f(17, big_from(h), 4) == log2_f);
  }
  CHECK(outer_scale_factor(0, std::uint64_t{3}, 4).log2 == 29);  // P = 1, b = 2
  CHECK(outer_scale_factor(0, std::uint64_t{4}, 4).log2 == 31);  // P = 1, b = 3
}

TEST_CASE("outer scales numeric: F on the Mode B rows of the threshold table",
          "[paper][outer][scales_numeric]") {
  // (log2 P, h = b + 1, log2 F) from thresholds.md Sec 5.2, K = 16.
  struct Row {
    unsigned p;
    std::uint64_t h;
    std::uint64_t log2_f;
  };
  const Row rows[] = {{26, 21, 64}, {25, 401, 75}, {20, 160001, 96}, {27, 28, 66}};
  for (const Row& r : rows) {
    const auto f = outer_scale_factor(r.p, r.h, 4);
    CHECK(f.log2 == r.log2_f);
    CHECK(f.log2 == reference_log2_f(r.p, big_from(r.h), 4));
    CHECK(f.value == kU64Max);
    CHECK_FALSE(f.value_exact());
  }
}

TEST_CASE("outer scales numeric: the saturated value starts exactly at 2^64",
          "[paper][outer][scales_numeric]") {
  // h = 2: h^2 (2h + 1) = 20, so F = 2^(7 + p + 4k + 5).
  const auto below = outer_scale_factor(51, std::uint64_t{2}, 0);
  CHECK(below.log2 == 63);
  CHECK(below.value == (std::uint64_t{1} << 63));
  CHECK(below.value_exact());
  const auto at = outer_scale_factor(52, std::uint64_t{2}, 0);
  CHECK(at.log2 == 64);
  CHECK(at.value == kU64Max);
  CHECK_FALSE(at.value_exact());
}

TEST_CASE("outer scales numeric: F equals the multi-precision reference",
          "[paper][outer][scales_numeric]") {
  std::mt19937_64 rng(20261127);
  for (int it = 0; it < 400; ++it) {
    const auto p = static_cast<unsigned>(rng() % 41);
    const auto k = static_cast<unsigned>(rng() % 26);
    const unsigned bits = 2 + static_cast<unsigned>(rng() % 63);  // 2..64
    std::uint64_t h = rng();
    if (bits < 64) h &= (std::uint64_t{1} << bits) - 1;
    if (h < 2) h = 2;
    if (it == 0) h = kU64Max;
    if (it == 1) h = std::uint64_t{1} << 63;
    const auto f = outer_scale_factor(p, h, k);
    const std::uint64_t want = reference_log2_f(p, big_from(h), k);
    CAPTURE(p, k, h);
    REQUIRE(f.log2 == want);
    REQUIRE(f.value == (want < 64 ? std::uint64_t{1} << want : kU64Max));
  }
}

TEST_CASE("outer scales numeric: F rejects h below 2", "[paper][outer][scales_numeric]") {
  // h = b + 1 with b >= 1 (Sec 7, l.1644).
  CHECK_THROWS_AS(outer_scale_factor(3, std::uint64_t{0}, 4), std::invalid_argument);
  CHECK_THROWS_AS(outer_scale_factor(3, std::uint64_t{1}, 4), std::invalid_argument);
  CHECK(outer_scale_factor(3, std::uint64_t{2}, 4).log2 == 31);
}

TEST_CASE("outer scales numeric: scale_pow2 doubles and saturates",
          "[paper][outer][scales_numeric]") {
  CHECK(scale_pow2<std::uint64_t>(5, 0) == 5);
  CHECK(scale_pow2<std::uint64_t>(5, 3) == 40);
  CHECK(scale_pow2<std::uint64_t>(1, 63) == (std::uint64_t{1} << 63));
  CHECK(scale_pow2<std::uint64_t>(1, 64) == kU64Max);
  CHECK(scale_pow2<std::uint64_t>(3, 62) == (std::uint64_t{3} << 62));
  CHECK(scale_pow2<std::uint64_t>(3, 63) == kU64Max);
  CHECK(scale_pow2<std::uint64_t>(1, kU64Max) == kU64Max);
  CHECK(scale_pow2<std::uint64_t>(0, kU64Max) == 0);
  CHECK(scale_pow2<std::uint64_t>(kU64Max, 0) == kU64Max);

  CHECK(scale_pow2<double>(0.75, 3) == 6.0);
  CHECK(scale_pow2<double>(1.0, 1023) == 0x1p1023);
  CHECK(scale_pow2<double>(1.0, 1024) == kDInf);
  CHECK(scale_pow2<double>(std::numeric_limits<double>::denorm_min(), 1074) == 1.0);
  CHECK(scale_pow2<double>(0.0, kU64Max) == 0.0);
  CHECK(scale_pow2<double>(0x1p-1000, kU64Max) == kDInf);
}

TEST_CASE("outer scales numeric: s_U is the least radius", "[paper][outer][scales_numeric]") {
  const std::vector<std::uint64_t> r{9, 4, 7, 4, 12};
  CHECK(outer_block_scale(sp(r)) == 4);
  const std::vector<double> rd{2.5, 0.375, 8.0};
  CHECK(outer_block_scale(sp(rd)) == 0.375);
  // A positive ordinary block (l.1999): empty, zero and infinite radii are rejected.
  const std::vector<std::uint64_t> none;
  const std::vector<std::uint64_t> zero{3, 0};
  const std::vector<std::uint64_t> inf{3, kU64Max};
  const std::vector<double> nan{1.0, std::numeric_limits<double>::quiet_NaN()};
  const std::vector<double> dinf{1.0, kDInf};
  CHECK_THROWS_AS(outer_block_scale(sp(none)), std::invalid_argument);
  CHECK_THROWS_AS(outer_block_scale(sp(zero)), std::invalid_argument);
  CHECK_THROWS_AS(outer_block_scale(sp(inf)), std::invalid_argument);
  CHECK_THROWS_AS(outer_block_scale(sp(nan)), std::invalid_argument);
  CHECK_THROWS_AS(outer_block_scale(sp(dinf)), std::invalid_argument);
}

TEST_CASE("outer scales numeric: aligned iff F^2 s >= a, at the boundary",
          "[paper][outer][scales_numeric]") {
  const auto f = small_f();  // F^2 = 2^32
  const std::uint64_t f2 = std::uint64_t{1} << 32;
  CHECK(outer_block_aligned<std::uint64_t>(f, 1, f2));
  CHECK_FALSE(outer_block_aligned<std::uint64_t>(f, 1, f2 + 1));
  CHECK(outer_block_aligned<std::uint64_t>(f, 3, 3 * f2));
  CHECK_FALSE(outer_block_aligned<std::uint64_t>(f, 3, 3 * f2 + 1));
  CHECK(outer_block_aligned<std::uint64_t>(f, 7, 5));
  // F^2 s does not fit: the saturated product still decides exactly.
  CHECK(outer_block_aligned<std::uint64_t>(f, std::uint64_t{1} << 40, kU64Max - 1));
  CHECK(outer_block_aligned<std::uint64_t>(paper_f(), 1, kU64Max - 1));

  CHECK(outer_block_aligned<double>(f, 0x1p-40, 0x1p-8));
  CHECK_FALSE(outer_block_aligned<double>(f, 0x1p-40, 0x1.0000000000001p-8));
  CHECK(outer_block_aligned<double>(f, 0.75, 0.75 * 0x1p32));
  CHECK(outer_block_aligned<double>(paper_f(), std::numeric_limits<double>::denorm_min(),
                                    std::numeric_limits<double>::max()));

  CHECK_THROWS_AS(outer_block_aligned<std::uint64_t>(f, 0, 5), std::invalid_argument);
  CHECK_THROWS_AS(outer_block_aligned<std::uint64_t>(f, 5, 0), std::invalid_argument);
  CHECK_THROWS_AS(outer_block_aligned<std::uint64_t>(f, kU64Max, 5), std::invalid_argument);
  CHECK_THROWS_AS(outer_block_aligned<std::uint64_t>(f, 5, kU64Max), std::invalid_argument);
  CHECK_THROWS_AS(outer_block_aligned<double>(f, 1.0, kDInf), std::invalid_argument);
  CHECK_THROWS_AS(outer_block_aligned<double>(f, -1.0, 1.0), std::invalid_argument);
}

TEST_CASE("outer scales numeric: aligned range a/F^2 <= s <= F a",
          "[paper][outer][scales_numeric]") {
  const auto f = small_f();  // F = 2^16
  CHECK(outer_aligned_scale_in_range<std::uint64_t>(f, 65536, 1));
  CHECK_FALSE(outer_aligned_scale_in_range<std::uint64_t>(f, 65537, 1));
  CHECK(outer_aligned_scale_in_range<std::uint64_t>(f, 1, std::uint64_t{1} << 32));
  CHECK_FALSE(outer_aligned_scale_in_range<std::uint64_t>(f, 1, (std::uint64_t{1} << 32) + 1));
  CHECK(outer_aligned_scale_in_range<double>(f, 0x1p16, 1.0));
  CHECK_FALSE(outer_aligned_scale_in_range<double>(f, 0x1.0000000000001p16, 1.0));
}

TEST_CASE("outer scales numeric: unaligned MST edges are below a_U / 2",
          "[paper][outer][scales_numeric]") {
  const std::vector<std::uint64_t> w{1, 2, 3};
  CHECK(first_edge_not_below_half<std::uint64_t>(sp(w), 7) == 3);  // 2 * 3 = 6 < 7
  CHECK(first_edge_not_below_half<std::uint64_t>(sp(w), 6) == 2);  // 2 * 3 = 6 is not < 6
  CHECK(first_edge_not_below_half<std::uint64_t>(sp(w), 2) == 0);
  const std::vector<std::uint64_t> none;
  CHECK(first_edge_not_below_half<std::uint64_t>(sp(none), 1) == 0);
  // 2w does not fit: w is not below a / 2 for any finite a.
  const std::vector<std::uint64_t> big{5, std::uint64_t{1} << 63};
  CHECK(first_edge_not_below_half<std::uint64_t>(sp(big), kU64Max - 1) == 1);
  const std::vector<double> d{0.25, 0.5};
  CHECK(first_edge_not_below_half<double>(sp(d), 1.25) == 2);
  CHECK(first_edge_not_below_half<double>(sp(d), 1.0) == 1);
  CHECK_THROWS_AS(first_edge_not_below_half<std::uint64_t>(sp(w), kU64Max), std::invalid_argument);
}

TEST_CASE("outer scales numeric: window bound 4 K^2 h s", "[paper][outer][scales_numeric]") {
  CHECK(outer_window_edge_bound<std::uint64_t>(3, 5, 4) == 4 * 256 * 5 * 3);
  CHECK(outer_window_edge_bound<std::uint64_t>(1, 2, 0) == 8);
  CHECK(outer_window_edge_bound<std::uint64_t>(std::uint64_t{1} << 50, 1000, 4) == kU64Max);
  CHECK(outer_window_edge_bound<double>(0.375, 5, 4) == 4.0 * 256 * 5 * 0.375);
  CHECK_THROWS_AS(outer_window_edge_bound<std::uint64_t>(3, 1, 4), std::invalid_argument);
}

TEST_CASE("outer scales numeric: cutoff range s/2 < t <= s", "[paper][outer][scales_numeric]") {
  // Q-36: the outer range is half-open at s/2.
  CHECK(outer_cutoff_valid<std::uint64_t>(8, 8));
  CHECK(outer_cutoff_valid<std::uint64_t>(8, 5));
  CHECK_FALSE(outer_cutoff_valid<std::uint64_t>(8, 4));
  CHECK_FALSE(outer_cutoff_valid<std::uint64_t>(8, 9));
  CHECK(outer_cutoff_valid<std::uint64_t>(9, 5));
  CHECK(outer_cutoff_valid<std::uint64_t>(kU64Max - 1, std::uint64_t{1} << 63));
  CHECK(outer_cutoff_valid<double>(1.0, 0.75));
  CHECK_FALSE(outer_cutoff_valid<double>(1.0, 0.5));
}

TEST_CASE("outer scales numeric: supports, merged grid and sentinels on a hand example",
          "[paper][outer][scales_numeric]") {
  const auto f = small_f();  // F^3 = 2^48
  const std::uint64_t one = 1;
  SECTION("a gap of exactly two is not merged and the sentinel is shared") {
    const std::vector<std::uint64_t> names{one, one << 49};
    const auto g = build_outer_global_grid<std::uint64_t>(f, sp(names), {});
    REQUIRE(g.supports.size() == 2);
    CHECK(g.supports[0].lower == 1);
    CHECK(g.supports[0].upper == one << 48);
    CHECK(g.supports[1].lower == one << 49);
    CHECK(g.supports[1].upper == kU64Max);  // 2^97 saturates (D-13)
    // 1, 2, ..., 2^48, then 2^49 (sentinel of the first grid and start of the second), ...,
    // 2^63, and the saturated sentinel.
    REQUIRE(g.points.size() == 65);
    for (std::size_t i = 0; i < 64; ++i) CHECK(g.points[i] == one << i);
    CHECK(g.points[64] == kU64Max);
    CHECK(g.component_begin == std::vector<std::size_t>{0, 49});
    CHECK(g.supports[0].lower_index == 0);
    CHECK(g.supports[0].upper_index == 49);
    CHECK(g.supports[1].lower_index == 49);
    CHECK(g.supports[1].upper_index == 64);
    CHECK(outer_grid_class<std::uint64_t>(g, 0, 5) == 2);
    CHECK(outer_grid_class<std::uint64_t>(g, 0, (one << 49) - 1) == 48);
    CHECK(outer_grid_class<std::uint64_t>(g, 0, one << 49) == kNoGridClass);
    CHECK(outer_grid_class<std::uint64_t>(g, 1, (one << 49) - 1) == kNoGridClass);
    CHECK(outer_grid_class<std::uint64_t>(g, 1, kU64Max - 1) == 63);
  }
  SECTION("a gap below two is merged into one doubling grid") {
    const std::vector<std::uint64_t> names{(one << 49) - 1, one};
    const auto g = build_outer_global_grid<std::uint64_t>(f, sp(names), {});
    CHECK(g.component_begin == std::vector<std::size_t>{0});
    REQUIRE(g.points.size() == 65);
    CHECK(g.points[49] == one << 49);
    CHECK(g.supports[0].lower_index == 48);  // floor of 2^49 - 1
    CHECK(g.supports[1].lower_index == 0);
    CHECK(g.supports[1].upper_index == 49);
  }
  SECTION("an aligned block lowers the start of its name's interval") {
    const std::vector<std::uint64_t> names{1024, one << 62};
    const std::vector<OuterAlignedScale<std::uint64_t>> blocks{{100, 0}, {3, 0}, {2000, 0}};
    const auto g = build_outer_global_grid<std::uint64_t>(
        f, sp(names), std::span<const OuterAlignedScale<std::uint64_t>>(blocks));
    CHECK(g.supports[0].lower == 3);
    CHECK(g.supports[0].upper == 1024 * (one << 48));
    CHECK(g.points[0] == 3);
    CHECK(g.points[5] == 96);
    // t_U is the floor of s_U on the grid of its name (l.2090).
    CHECK(outer_aligned_cutoff<std::uint64_t>(g, 0, 100) == 96);
    CHECK(outer_aligned_cutoff<std::uint64_t>(g, 0, 3) == 3);
    CHECK(outer_aligned_cutoff<std::uint64_t>(g, 0, 2000) == 1536);
    CHECK(outer_grid_class<std::uint64_t>(g, 0, 2) == kNoGridClass);
    // A scale below the support has no floor: a violated relation of the paper.
    CHECK_THROWS_AS(outer_aligned_cutoff<std::uint64_t>(g, 0, 2), std::logic_error);
    CHECK_THROWS_AS(outer_aligned_cutoff<std::uint64_t>(g, 7, 100), std::invalid_argument);
  }
  SECTION("no names: empty grid") {
    const std::vector<std::uint64_t> names;
    const auto g = build_outer_global_grid<std::uint64_t>(f, sp(names), {});
    CHECK(g.points.empty());
    CHECK(g.supports.empty());
    CHECK(g.component_begin.empty());
  }
}

TEST_CASE("outer scales numeric: the grid builder rejects inputs outside Sec 8.2",
          "[paper][outer][scales_numeric]") {
  const auto f = small_f();
  using Blocks = std::vector<OuterAlignedScale<std::uint64_t>>;
  const auto build = [&](const std::vector<std::uint64_t>& names, const Blocks& blocks) {
    return build_outer_global_grid<std::uint64_t>(
        f, sp(names), std::span<const OuterAlignedScale<std::uint64_t>>(blocks));
  };
  CHECK_THROWS_AS(build({0}, {}), std::invalid_argument);
  CHECK_THROWS_AS(build({kU64Max}, {}), std::invalid_argument);
  CHECK_THROWS_AS(build({8}, {{4, 1}}), std::invalid_argument);        // unknown name
  CHECK_THROWS_AS(build({8}, {{0, 0}}), std::invalid_argument);        // s = 0
  CHECK_THROWS_AS(build({8}, {{kU64Max, 0}}), std::invalid_argument);  // s infinite
  // F^2 s < a: the block is not aligned and must not touch the global supports.
  CHECK_THROWS_AS(build({(std::uint64_t{1} << 32) + 1}, {{1, 0}}), std::invalid_argument);
  // s > F a contradicts Lemma 8.4 (s_U <= F a_U, l.2118).
  CHECK_THROWS_AS(build({1}, {{65537, 0}}), std::logic_error);
  CHECK_NOTHROW(build({1}, {{65536, 0}}));

  const std::vector<double> dn{1.0, std::numeric_limits<double>::quiet_NaN()};
  CHECK_THROWS_AS(build_outer_global_grid<double>(f, sp(dn), {}), std::invalid_argument);
}

TEST_CASE("outer scales numeric: double supports saturate at infinity",
          "[paper][outer][scales_numeric]") {
  const std::vector<double> names{0x1p1000, 0.5};
  const auto g = build_outer_global_grid<double>(paper_f(), sp(names), {});
  CHECK(g.supports[0].upper == kDInf);
  CHECK(g.supports[1].upper == kDInf);
  CHECK(g.component_begin == std::vector<std::size_t>{0});
  // 2^-1, ..., 2^1023 and the saturated sentinel.
  REQUIRE(g.points.size() == 1026);
  CHECK(g.points.front() == 0.5);
  CHECK(g.points[1024] == 0x1p1023);
  CHECK(g.points.back() == kDInf);
  CHECK(g.supports[0].lower_index == 1001);
  CHECK(g.supports[0].upper_index == 1025);
  CHECK(outer_grid_class<double>(g, 1, 3.0) == 2);
  CHECK(outer_grid_class<double>(g, 0, std::numeric_limits<double>::max()) == 1024);
}

TEST_CASE("outer scales numeric: separate-grid cutoff", "[paper][outer][scales_numeric]") {
  const std::vector<std::uint64_t> table{4, 8, 16, 32};
  CHECK(outer_separate_cutoff<std::uint64_t>(sp(table), 20) == 16);
  CHECK(outer_separate_cutoff<std::uint64_t>(sp(table), 4) == 4);
  CHECK(outer_separate_cutoff<std::uint64_t>(sp(table), 63) == 32);
  // No floor, or a floor at most s/2: the table does not cover s (l.2100-2103).
  CHECK_THROWS_AS(outer_separate_cutoff<std::uint64_t>(sp(table), 3), std::logic_error);
  CHECK_THROWS_AS(outer_separate_cutoff<std::uint64_t>(sp(table), 64), std::logic_error);
  const std::vector<double> dt{0.25, 0.5, 1.0};
  CHECK(outer_separate_cutoff<double>(sp(dt), 0.75) == 0.5);
}
