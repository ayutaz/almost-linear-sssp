// Sec 2, "Operations on lengths and indices" (arXiv 2611.34911v1, txt l.292-300): lengths are
// opaque reals (copy, compare, add); infinity is a separate flag; a signed value is a pair of
// nonnegative sums, a-b < c-d is decided by a+d < c+b; multiplication by 2^j is j doublings.
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "paper/model/length.hpp"
#include "sssp/result.hpp"

using sssp::OpCounters;
using sssp::paper::model::EdgeKey;
using sssp::paper::model::Len;
using sssp::paper::model::LenOps;
using sssp::paper::model::SignedPair;

namespace {
constexpr std::uint64_t kU64Max = std::numeric_limits<std::uint64_t>::max();
constexpr bool kCounting = SSSP_DIAGNOSTICS != 0;
__extension__ typedef __int128 i128;

// Deterministic test-local generator (xorshift64*), independent of the library's RNG.
struct TestGen {
  std::uint64_t s;
  std::uint64_t next() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * 2685821657736338717ull;
  }
};
}  // namespace

TEMPLATE_TEST_CASE("infinity flag absorbs addition", "[paper][model][length]", std::uint64_t,
                   double) {
  using W = TestType;
  OpCounters c;
  const LenOps<W> ops(&c);
  const auto inf = Len<W>::infinity();
  const auto five = Len<W>::finite(W{5});
  REQUIRE(inf.is_infinite());
  REQUIRE_FALSE(five.is_infinite());
  REQUIRE(ops.add(inf, five).is_infinite());
  REQUIRE(ops.add(five, inf).is_infinite());
  REQUIRE(ops.add(inf, inf).is_infinite());
  REQUIRE(c.additions == 0);  // no real register is touched when a flag is set
  const auto ten = ops.add(five, five);
  REQUIRE_FALSE(ten.is_infinite());
  REQUIRE(ten.value() == W{10});
  REQUIRE(c.additions == (kCounting ? 1u : 0u));
  REQUIRE_THROWS_AS(inf.value(), std::logic_error);
  REQUIRE(Len<W>::zero().value() == W{0});
  REQUIRE(Len<W>{}.value() == W{0});
}

TEMPLATE_TEST_CASE("comparisons order every finite length below infinity", "[paper][model][length]",
                   std::uint64_t, double) {
  using W = TestType;
  OpCounters c;
  const LenOps<W> ops(&c);
  const auto inf = Len<W>::infinity();
  const auto a = Len<W>::finite(W{3});
  const auto b = Len<W>::finite(W{7});
  REQUIRE(ops.less(a, inf));
  REQUIRE_FALSE(ops.less(inf, a));
  REQUIRE_FALSE(ops.less(inf, inf));
  REQUIRE(ops.equal(inf, inf));
  REQUIRE_FALSE(ops.equal(a, inf));
  REQUIRE(ops.less_equal(inf, inf));
  REQUIRE(ops.less_equal(a, inf));
  REQUIRE_FALSE(ops.less_equal(inf, a));
  REQUIRE(c.comparisons == 0);
  REQUIRE(ops.less(a, b));
  REQUIRE_FALSE(ops.less(b, a));
  REQUIRE_FALSE(ops.less(a, a));
  REQUIRE(ops.less_equal(a, a));
  REQUIRE(ops.equal(a, a));
  REQUIRE_FALSE(ops.equal(a, b));
  REQUIRE(c.comparisons == (kCounting ? 6u : 0u));  // one per call with two finite operands
}

TEST_CASE("integer lengths detect overflow and keep 2^64-1 for the infinity sentinel",
          "[paper][model][length]") {
  using L64 = Len<std::uint64_t>;
  const LenOps<std::uint64_t> ops;
  REQUIRE_THROWS_AS(L64::finite(kU64Max), std::invalid_argument);
  const auto big = L64::finite(kU64Max - 1);
  REQUIRE(ops.add(L64::finite(kU64Max - 2), L64::finite(1)).value() == kU64Max - 1);
  REQUIRE_THROWS_AS(ops.add(big, L64::finite(1)), std::overflow_error);
  REQUIRE_THROWS_AS(ops.add(big, big), std::overflow_error);
  REQUIRE(ops.add(big, L64::zero()).value() == kU64Max - 1);
}

TEST_CASE("double lengths reject invalid values and detect overflow", "[paper][model][length]") {
  using LD = Len<double>;
  const LenOps<double> ops;
  REQUIRE_THROWS_AS(LD::finite(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
  REQUIRE_THROWS_AS(LD::finite(std::numeric_limits<double>::infinity()), std::invalid_argument);
  REQUIRE_THROWS_AS(LD::finite(-1.0), std::invalid_argument);
  REQUIRE_FALSE(std::signbit(LD::finite(-0.0).value()));  // -0.0 is normalised to +0.0
  const auto huge = LD::finite(std::numeric_limits<double>::max());
  REQUIRE_THROWS_AS(ops.add(huge, huge), std::overflow_error);
}

TEMPLATE_TEST_CASE("weights convert to and from the shared infinity sentinel",
                   "[paper][model][length]", std::uint64_t, double) {
  using W = TestType;
  const W sentinel = sssp::WeightTraits<W>::infinity();
  REQUIRE(Len<W>::from_weight(sentinel).is_infinite());
  REQUIRE(Len<W>::infinity().to_weight() == sentinel);
  REQUIRE(Len<W>::from_weight(W{42}).value() == W{42});
  REQUIRE(Len<W>::finite(W{42}).to_weight() == W{42});
}

TEST_CASE("mul_pow2 is exact for double and costs j doublings", "[paper][model][length]") {
  TestGen gen{0x9e3779b97f4a7c15ull};
  for (int trial = 0; trial < 2000; ++trial) {
    const double x = static_cast<double>(gen.next() >> 11) * 0x1.0p-40;  // in [0, 2^13)
    const unsigned j = static_cast<unsigned>(gen.next() % 200);
    OpCounters c;
    const LenOps<double> ops(&c);
    const auto y = ops.mul_pow2(Len<double>::finite(x), j);
    REQUIRE(y.value() == std::ldexp(x, static_cast<int>(j)));
    REQUIRE(c.additions == (kCounting ? j : 0u));
  }
  const LenOps<double> ops;
  REQUIRE(ops.mul_pow2(Len<double>::finite(3.0), 10).value() == 3072.0);
  REQUIRE(ops.mul_pow2(Len<double>::finite(0.1), 3).value() == std::ldexp(0.1, 3));
  REQUIRE(ops.mul_pow2(Len<double>::finite(5.0), 0).value() == 5.0);
  REQUIRE(ops.mul_pow2(Len<double>::infinity(), 7).is_infinite());
  REQUIRE_THROWS_AS(ops.mul_pow2(Len<double>::finite(1.0), 1024), std::overflow_error);
}

TEST_CASE("mul_pow2 on integers is exact and detects overflow", "[paper][model][length]") {
  OpCounters c;
  const LenOps<std::uint64_t> ops(&c);
  using L64 = Len<std::uint64_t>;
  REQUIRE(ops.mul_pow2(L64::finite(1), 63).value() == (std::uint64_t{1} << 63));
  REQUIRE(c.additions == (kCounting ? 63u : 0u));
  REQUIRE(ops.mul_pow2(L64::finite(3), 61).value() == (std::uint64_t{3} << 61));
  REQUIRE(ops.mul_pow2(L64::finite(3), 62).value() == (std::uint64_t{3} << 62));  // < 2^64 - 1
  REQUIRE_THROWS_AS(ops.mul_pow2(L64::finite(1), 64), std::overflow_error);
  REQUIRE_THROWS_AS(ops.mul_pow2(L64::finite(3), 63), std::overflow_error);  // 3 * 2^63 > 2^64
  REQUIRE(ops.mul_pow2(L64::infinity(), 100).is_infinite());
  REQUIRE(ops.mul_pow2(L64::zero(), 100).value() == 0);
}

TEMPLATE_TEST_CASE("signed pair compares a-b vs c-d via a+d vs c+b", "[paper][model][length]",
                   std::uint64_t, double) {
  using W = TestType;
  using SP = SignedPair<W>;
  OpCounters c;
  const LenOps<W> ops(&c);
  const auto len = [](std::uint64_t v) { return Len<W>::finite(static_cast<W>(v)); };
  const SP two(len(5), len(3));     // 5 - 3 = 2
  const SP three(len(10), len(7));  // 10 - 7 = 3
  const SP minus4(len(1), len(5));  // 1 - 5 = -4, not representable as an unsigned difference
  REQUIRE(ops.less(two, three));
  REQUIRE_FALSE(ops.less(three, two));
  REQUIRE(ops.less(minus4, SP::from(len(0))));
  REQUIRE(ops.equal(SP(len(5), len(3)), SP(len(4), len(2))));
  REQUIRE_FALSE(ops.less(SP(len(5), len(3)), SP(len(4), len(2))));
  // each comparison: two additions (a+d and c+b) and one real comparison (Sec 2, l.293-296)
  REQUIRE(c.additions == (kCounting ? 10u : 0u));
  REQUIRE(c.comparisons == (kCounting ? 5u : 0u));

  TestGen gen{12345};
  for (int trial = 0; trial < 20000; ++trial) {
    const std::uint64_t v[4] = {gen.next() >> 24, gen.next() >> 24, gen.next() >> 24,
                                gen.next() >> 24};  // < 2^40, so all sums are exact in double
    const SP x(len(v[0]), len(v[1]));
    const SP y(len(v[2]), len(v[3]));
    const i128 dx = static_cast<i128>(v[0]) - static_cast<i128>(v[1]);
    const i128 dy = static_cast<i128>(v[2]) - static_cast<i128>(v[3]);
    REQUIRE(ops.less(x, y) == (dx < dy));
    REQUIRE(ops.equal(x, y) == (dx == dy));
  }
}

TEMPLATE_TEST_CASE("signed pair addition is componentwise and negation swaps",
                   "[paper][model][length]", std::uint64_t, double) {
  using W = TestType;
  using SP = SignedPair<W>;
  OpCounters c;
  const LenOps<W> ops(&c);
  const auto len = [](std::uint64_t v) { return Len<W>::finite(static_cast<W>(v)); };
  const SP x(len(9), len(2));  // 7
  const SP y(len(1), len(6));  // -5
  const SP s = ops.add(x, y);  // 2
  REQUIRE(s.positive().value() == W{10});
  REQUIRE(s.negative().value() == W{8});
  REQUIRE(c.additions == (kCounting ? 2u : 0u));
  const SP nx = x.neg();  // -7
  REQUIRE(nx.positive().value() == W{2});
  REQUIRE(nx.negative().value() == W{9});
  REQUIRE(ops.equal(ops.add(x, nx), SP::from(len(0))));
  REQUIRE(ops.equal(ops.sub(x, y), SP::from(len(12))));
  REQUIRE(SP{}.positive().value() == W{0});
  REQUIRE(SP{}.negative().value() == W{0});
}

TEMPLATE_TEST_CASE("signed pair rejects infinite components", "[paper][model][length]",
                   std::uint64_t, double) {
  using W = TestType;
  const auto inf = Len<W>::infinity();
  const auto one = Len<W>::finite(W{1});
  REQUIRE_THROWS_AS(SignedPair<W>(inf, one), std::invalid_argument);
  REQUIRE_THROWS_AS(SignedPair<W>(one, inf), std::invalid_argument);
  REQUIRE_THROWS_AS(SignedPair<W>::from(inf), std::invalid_argument);
}

TEMPLATE_TEST_CASE("EdgeKey ties broken by index", "[paper][model][length]", std::uint64_t,
                   double) {
  using W = TestType;
  using K = EdgeKey<W>;
  OpCounters c;
  const LenOps<W> ops(&c);
  REQUIRE(ops.less(K{W{3}, 9}, K{W{4}, 1}));
  REQUIRE(c.comparisons == (kCounting ? 1u : 0u));  // decided by the weight alone
  REQUIRE(ops.less(K{W{3}, 1}, K{W{3}, 2}));
  REQUIRE(c.comparisons == (kCounting ? 3u : 0u));  // equal weights need a second comparison
  REQUIRE_FALSE(ops.less(K{W{3}, 2}, K{W{3}, 2}));
  REQUIRE(K{W{3}, 1} < K{W{3}, 2});
  REQUIRE_FALSE(K{W{3}, 2} < K{W{3}, 1});
  REQUIRE(K{W{3}, 2} == K{W{3}, 2});
  REQUIRE_FALSE(K{W{3}, 2} == K{W{3}, 1});
}

TEMPLATE_TEST_CASE("EdgeKey is a strict total order", "[paper][model][length]", std::uint64_t,
                   double) {
  using W = TestType;
  using K = EdgeKey<W>;
  std::vector<K> keys;
  for (std::uint32_t id = 0; id < 40; ++id) keys.push_back(K{static_cast<W>(id % 5), id});
  if constexpr (std::is_floating_point_v<W>) keys.push_back(K{-0.0, 100});  // equals +0.0
  const LenOps<W> ops;
  for (const auto& a : keys) {
    REQUIRE_FALSE(ops.less(a, a));
    for (const auto& b : keys) {
      const int relations = int{ops.less(a, b)} + int{ops.less(b, a)} + int{a == b};
      REQUIRE(relations == 1);  // trichotomy
      REQUIRE(ops.less(a, b) == (a < b));
      for (const auto& d : keys) {
        if (ops.less(a, b) && ops.less(b, d)) REQUIRE(ops.less(a, d));
      }
    }
  }
}

TEST_CASE("LenOps without a counter sink performs the same operations", "[paper][model][length]") {
  const LenOps<std::uint64_t> silent;
  REQUIRE(silent.counters() == nullptr);
  REQUIRE(silent.add(Len<std::uint64_t>::finite(2), Len<std::uint64_t>::finite(3)).value() == 5);
  REQUIRE(silent.less(Len<std::uint64_t>::finite(2), Len<std::uint64_t>::finite(3)));
}
