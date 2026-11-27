// Sec 2 (txt l.284-287, l.304): L_k = ceil(log2(k+2)); outer construction needs N > L_N^64,
// local construction needs N > 2 L_N^4; size-test products may be capped at the tested limit.
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "paper/model/params.hpp"

namespace model = sssp::paper::model;
using model::capped_pow;
using model::L;
using model::local_nontrivial;
using model::outer_nontrivial;
using model::PaperParams;
using model::saturating_pow;

namespace {
constexpr std::uint64_t kMax = std::numeric_limits<std::uint64_t>::max();
__extension__ typedef unsigned __int128 u128;

// Independent reference: exact value of coeff * base^e in 128 bits, clamped to 2^100 so that it
// never wraps (every value compared against it is below 2^64).
u128 exact_scaled_pow(std::uint64_t coeff, std::uint64_t base, unsigned e) {
  const u128 clamp = static_cast<u128>(1) << 100;
  u128 acc = coeff;
  for (unsigned i = 0; i < e && acc < clamp; ++i) acc *= base;
  return acc < clamp ? acc : clamp;
}

// Independent reference for L_k: least c with k + 2 <= 2^c, evaluated in 128 bits.
unsigned reference_L(std::uint64_t k) {
  const u128 target = static_cast<u128>(k) + 2;
  unsigned c = 0;
  while ((static_cast<u128>(1) << c) < target) ++c;
  return c;
}
}  // namespace

TEST_CASE("L_k matches ceil(log2(k+2))", "[paper][model][params]") {
  REQUIRE(L(1) == 2);  // log2 3
  REQUIRE(L(2) == 2);  // log2 4
  REQUIRE(L(3) == 3);  // log2 5
  REQUIRE(L(6) == 3);  // log2 8
  REQUIRE(L(7) == 4);  // log2 9
  REQUIRE(L(209952) == 18);
  for (std::uint64_t k = 1; k <= 100000; ++k) REQUIRE(L(k) == reference_L(k));
  for (unsigned b = 1; b < 64; ++b) {
    const std::uint64_t p = std::uint64_t{1} << b;
    for (const std::uint64_t k : {p - 3, p - 2, p - 1, p, p + 1}) {
      if (k >= 1) REQUIRE(L(k) == reference_L(k));
    }
  }
  REQUIRE(L(kMax - 2) == 64);  // k + 2 = 2^64 - 1
  REQUIRE(L(kMax - 1) == 64);  // k + 2 = 2^64
  REQUIRE(L(kMax) == 65);      // k + 2 = 2^64 + 1
}

TEST_CASE("Mode A threshold constants are the paper's and immutable", "[paper][model][params]") {
  STATIC_REQUIRE(model::kOuterExponent == 64);
  STATIC_REQUIRE(model::kLocalExponent == 4);
  STATIC_REQUIRE(model::kLocalCoeff == 2);
  STATIC_REQUIRE(model::kModeAParams.outer_exponent == 64);
  STATIC_REQUIRE(model::kModeAParams.local_exponent == 4);
  STATIC_REQUIRE(model::kModeAParams.local_coeff == 2);
  STATIC_REQUIRE(std::is_const_v<decltype(model::kModeAParams)>);
}

TEST_CASE("saturating power never wraps", "[paper][model][params]") {
  REQUIRE(saturating_pow(2, 63) == (std::uint64_t{1} << 63));
  REQUIRE(saturating_pow(2, 64) == kMax);
  REQUIRE(saturating_pow(65, 64) == kMax);  // L_N^64 for N near 2^64
  REQUIRE(saturating_pow(3, 40) == 12157665459056928801ull);
  REQUIRE(saturating_pow(3, 41) == kMax);
  REQUIRE(saturating_pow(0, 0) == 1);
  REQUIRE(saturating_pow(0, 5) == 0);
  REQUIRE(saturating_pow(1, 1000) == 1);
  REQUIRE(saturating_pow(kMax, 1) == kMax);
  REQUIRE(saturating_pow(kMax, 2) == kMax);
  for (std::uint64_t base = 0; base <= 70; ++base) {
    for (unsigned e = 0; e <= 70; ++e) {
      const u128 exact = exact_scaled_pow(1, base, e);
      const std::uint64_t expected = exact > kMax ? kMax : static_cast<std::uint64_t>(exact);
      REQUIRE(saturating_pow(base, e) == expected);
    }
  }
}

TEST_CASE("capped power stops at the tested limit", "[paper][model][params]") {
  REQUIRE(capped_pow(18, 4, 1000000) == 104976);
  REQUIRE(capped_pow(18, 4, 104976) == 104976);
  REQUIRE(capped_pow(18, 4, 100000) == 100000);
  REQUIRE(capped_pow(65, 64, 12345) == 12345);
  REQUIRE(capped_pow(2, 0, 0) == 0);
  for (std::uint64_t base = 1; base <= 40; ++base) {
    for (unsigned e = 0; e <= 40; ++e) {
      for (const std::uint64_t cap : {std::uint64_t{1}, std::uint64_t{1000}, kMax / 3, kMax}) {
        const u128 exact = exact_scaled_pow(1, base, e);
        const std::uint64_t expected = exact >= cap ? cap : static_cast<std::uint64_t>(exact);
        REQUIRE(capped_pow(base, e, cap) == expected);
      }
    }
  }
}

TEST_CASE("local threshold first holds at 209953", "[paper][model][params]") {
  // L_209952 = L_209953 = 18 and 2 * 18^4 = 209952.
  REQUIRE_FALSE(local_nontrivial(209952));
  REQUIRE(local_nontrivial(209953));
  REQUIRE_FALSE(local_nontrivial(0));
  REQUIRE_FALSE(local_nontrivial(1));
  for (std::uint64_t n = 1; n <= 300000; ++n) {
    const bool expected = static_cast<u128>(n) > exact_scaled_pow(2, reference_L(n), 4);
    REQUIRE(local_nontrivial(n) == expected);
  }
  REQUIRE(local_nontrivial(kMax));
}

TEST_CASE("outer threshold never holds below 2^64", "[paper][model][params]") {
  // For every N < 2^64, L_N >= 2 so L_N^64 >= 2^64 > N.
  REQUIRE_FALSE(outer_nontrivial(0));
  for (std::uint64_t n = 1; n <= 100000; ++n) REQUIRE_FALSE(outer_nontrivial(n));
  for (unsigned b = 1; b < 64; ++b) {
    const std::uint64_t p = std::uint64_t{1} << b;
    REQUIRE_FALSE(outer_nontrivial(p - 1));
    REQUIRE_FALSE(outer_nontrivial(p));
    REQUIRE_FALSE(outer_nontrivial(p + 1));
  }
  REQUIRE_FALSE(outer_nontrivial(kMax));
}

TEST_CASE("Mode A overloads equal the parameterised tests with the paper's values",
          "[paper][model][params]") {
  for (std::uint64_t n = 0; n <= 250000; n += 7) {
    REQUIRE(local_nontrivial(n) == local_nontrivial(n, model::kModeAParams));
    REQUIRE(outer_nontrivial(n) == outer_nontrivial(n, model::kModeAParams));
  }
  REQUIRE(local_nontrivial(209953, model::kModeAParams));
  REQUIRE_FALSE(local_nontrivial(209952, model::kModeAParams));
}

TEST_CASE("Mode B thresholds follow injected parameters exactly", "[paper][model][params]") {
  for (unsigned outer_e = 1; outer_e <= 6; ++outer_e) {
    for (unsigned local_e = 1; local_e <= 6; ++local_e) {
      for (std::uint64_t coeff = 1; coeff <= 3; ++coeff) {
        const PaperParams p{outer_e, local_e, coeff};
        for (std::uint64_t n = 1; n <= 5000; ++n) {
          const u128 nn = n;
          const unsigned l = reference_L(n);
          REQUIRE(outer_nontrivial(n, p) == (nn > exact_scaled_pow(1, l, outer_e)));
          REQUIRE(local_nontrivial(n, p) == (nn > exact_scaled_pow(coeff, l, local_e)));
        }
      }
    }
  }
  // Large coefficients and exponents saturate instead of wrapping into a false "nontrivial".
  const PaperParams huge{64, 64, kMax};
  REQUIRE_FALSE(outer_nontrivial(kMax, huge));
  REQUIRE_FALSE(local_nontrivial(kMax, huge));
  const PaperParams coeff_only{1, 0, 7};  // N > 7 * L^0 = 7
  REQUIRE_FALSE(local_nontrivial(7, coeff_only));
  REQUIRE(local_nontrivial(8, coeff_only));
}
