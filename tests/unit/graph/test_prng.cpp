#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <vector>

#include "sssp/generators.hpp"

using sssp::gen::SplitMix64;
using sssp::gen::Xoshiro256ss;

// Known-answer vectors come from the reference algorithms (Vigna's splitmix64.c and
// xoshiro256starstar.c), computed independently of this implementation.
TEST_CASE("splitmix64 matches the reference sequence for seed 0", "[graph][prng]") {
  SplitMix64 sm(0);
  REQUIRE(sm.next() == 0xE220A8397B1DCDAFULL);
  REQUIRE(sm.next() == 0x6E789E6AA1B965F4ULL);
  REQUIRE(sm.next() == 0x06C45D188009454FULL);
}

TEST_CASE("xoshiro256** matches the reference sequence for state {1,2,3,4}", "[graph][prng]") {
  auto rng = Xoshiro256ss::from_state({1, 2, 3, 4});
  REQUIRE(rng.next() == 11520ULL);
  REQUIRE(rng.next() == 0ULL);
  REQUIRE(rng.next() == 1509978240ULL);
  REQUIRE(rng.next() == 1215971899390074240ULL);
}

TEST_CASE("xoshiro256** seeded from an integer is deterministic and seed-sensitive",
          "[graph][prng]") {
  Xoshiro256ss a(42), b(42), c(43);
  bool any_diff = false;
  for (int i = 0; i < 16; ++i) {
    const auto x = a.next();
    REQUIRE(x == b.next());
    any_diff = any_diff || (x != c.next());
  }
  REQUIRE(any_diff);
  // seeding goes through splitmix64, so seed 0 does not produce the all-zero state
  Xoshiro256ss z(0);
  REQUIRE((z.next() | z.next() | z.next()) != 0);
}

TEST_CASE("uniform_below stays in range and is roughly uniform", "[graph][prng]") {
  Xoshiro256ss rng(7);
  constexpr std::uint64_t kBuckets = 10;
  constexpr int kDraws = 100000;
  std::array<int, kBuckets> count{};
  for (int i = 0; i < kDraws; ++i) {
    const auto x = rng.uniform_below(kBuckets);
    REQUIRE(x < kBuckets);
    ++count[x];
  }
  // chi-square with 9 degrees of freedom; 40 is far beyond the 99.999% quantile (~33.7)
  double chi2 = 0;
  const double expected = static_cast<double>(kDraws) / kBuckets;
  for (int c : count) chi2 += (c - expected) * (c - expected) / expected;
  REQUIRE(chi2 < 40.0);
  REQUIRE(rng.uniform_below(1) == 0);
  // bound that is not a power of two near 2^64 still terminates and stays in range
  const std::uint64_t big = (1ULL << 63) + 12345;
  for (int i = 0; i < 100; ++i) REQUIRE(rng.uniform_below(big) < big);
}

TEST_CASE("uniform_range is inclusive at both ends", "[graph][prng]") {
  Xoshiro256ss rng(9);
  bool lo = false, hi = false;
  for (int i = 0; i < 1000; ++i) {
    const auto x = rng.uniform_range(3, 5);
    REQUIRE(x >= 3);
    REQUIRE(x <= 5);
    lo = lo || x == 3;
    hi = hi || x == 5;
  }
  REQUIRE(lo);
  REQUIRE(hi);
  REQUIRE(rng.uniform_range(0, ~0ULL) <= ~0ULL);  // full range does not overflow
}

TEST_CASE("uniform01 lies in [0,1) with mean near 1/2", "[graph][prng]") {
  Xoshiro256ss rng(11);
  double sum = 0;
  constexpr int kDraws = 100000;
  for (int i = 0; i < kDraws; ++i) {
    const double u = rng.uniform01();
    REQUIRE(u >= 0.0);
    REQUIRE(u < 1.0);
    sum += u;
  }
  REQUIRE(sum / kDraws == Catch::Approx(0.5).margin(0.01));
}

TEST_CASE("deterministic log/exp agree with libm to near machine precision", "[graph][prng]") {
  for (double x : {1e-300, 1e-10, 0.001, 0.5, 0.7071, 1.0, 1.5, 2.0, 10.0, 12345.678, 1e300}) {
    REQUIRE(sssp::gen::det_log(x) == Catch::Approx(std::log(x)).epsilon(1e-14).margin(1e-15));
  }
  for (double x : {-700.0, -20.0, -1.0, -1e-5, 0.0, 1e-5, 0.5, 1.0, 3.0, 50.0, 700.0}) {
    REQUIRE(sssp::gen::det_exp(x) == Catch::Approx(std::exp(x)).epsilon(1e-14));
  }
  REQUIRE(sssp::gen::det_log(1.0) == 0.0);
  REQUIRE(sssp::gen::det_exp(0.0) == 1.0);
}
