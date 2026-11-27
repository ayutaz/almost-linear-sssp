#pragma once
// Size parameters and construction thresholds of Sec 2 (arXiv 2611.34911v1, "Notation",
// txt l.284-287).
//
//   L_k = ceil(log2(k + 2))                     (l.284)
//   outer construction assumes N > L_N^64       (l.286; Alg 8 l.4 takes the base case otherwise)
//   local construction assumes N > 2 * L_N^4    (l.286; Alg 4 l.2 takes the base case otherwise)
//
// All of this is integer word arithmetic on sizes (never on real lengths). Powers are evaluated
// with products capped at the tested limit, as Sec 2 permits ("products used solely in integer
// size tests may be capped once they reach the tested limit", l.304), so no product can wrap.
//
// Mode A (src/paper/) uses the fixed constants below and the overloads without parameters.
// Mode B (src/experimental/paper_experimental/) may pass its own PaperParams; any value other
// than kModeAParams is a deviation that must be listed in docs/PAPER_MAPPING.md.

#include <cstdint>
#include <limits>

#include "paper/model/ilog.hpp"

namespace sssp::paper::model {

/// Threshold constants. The exponents and coefficient enter the constant relations used in the
/// proofs of Sec 5-8 (e.g. the bounds stated around txt l.1206), so they are part of the paper's
/// specification, not tuning knobs.
struct PaperParams {
  unsigned outer_exponent;    ///< e in "N > L_N^e" for the outer construction (paper: 64)
  unsigned local_exponent;    ///< e in "N > c * L_N^e" for the local construction (paper: 4)
  std::uint64_t local_coeff;  ///< c in "N > c * L_N^e" for the local construction (paper: 2)
};

inline constexpr unsigned kOuterExponent = 64;
inline constexpr unsigned kLocalExponent = 4;
inline constexpr std::uint64_t kLocalCoeff = 2;
/// The paper's values (Sec 2, l.286). Immutable; Mode A never uses anything else.
inline constexpr PaperParams kModeAParams{kOuterExponent, kLocalExponent, kLocalCoeff};

/// L_k = ceil(log2(k + 2)) (Sec 2, l.284). Defined for every 64-bit k (the paper states k >= 1;
/// the formula also gives L_0 = 1). For k + 2 >= 2, ceil(log2(k + 2)) = floor(log2(k + 1)) + 1,
/// which avoids forming k + 2; only k = 2^64 - 1 needs k + 1 > 2^64 - 1 and gives 65.
constexpr unsigned L(std::uint64_t k) {
  if (k == std::numeric_limits<std::uint64_t>::max()) return 65;
  return floor_log2(k + 1) + 1;
}

/// min(a * b, cap) without wrap-around.
constexpr std::uint64_t capped_mul(std::uint64_t a, std::uint64_t b, std::uint64_t cap) {
  std::uint64_t r = 0;
  if (__builtin_mul_overflow(a, b, &r) || r >= cap) return cap;
  return r;
}

/// min(base^e, cap), with the running product capped once it reaches cap (Sec 2, l.304). At most
/// 64 multiplications are performed for base >= 2.
constexpr std::uint64_t capped_pow(std::uint64_t base, unsigned e, std::uint64_t cap) {
  if (e == 0) return cap < 1 ? cap : 1;
  if (base == 0) return 0;
  if (base == 1) return cap < 1 ? cap : 1;
  std::uint64_t acc = cap < 1 ? cap : 1;
  for (unsigned i = 0; i < e && acc < cap; ++i) acc = capped_mul(acc, base, cap);
  return acc;
}

/// min(base^e, 2^64 - 1).
constexpr std::uint64_t saturating_pow(std::uint64_t base, unsigned e) {
  return capped_pow(base, e, std::numeric_limits<std::uint64_t>::max());
}

/// Exact truth value of N > coeff * base^e. The right-hand side is evaluated capped at N, so
/// a wrapped product can never make a small threshold appear (research note Sec 0 item 5).
constexpr bool exceeds_scaled_pow(std::uint64_t n, std::uint64_t coeff, std::uint64_t base,
                                  unsigned e) {
  const std::uint64_t rhs = capped_mul(coeff, capped_pow(base, e, n), n);
  return n > rhs;
}

/// Outer construction applies (Sec 2, l.286; Alg 8 l.4): N > L_N^{outer_exponent}.
constexpr bool outer_nontrivial(std::uint64_t n, const PaperParams& p) {
  return exceeds_scaled_pow(n, 1, L(n), p.outer_exponent);
}

/// Local construction applies (Sec 2, l.286; Alg 4 l.2): N > local_coeff * L_N^{local_exponent}.
constexpr bool local_nontrivial(std::uint64_t n, const PaperParams& p) {
  return exceeds_scaled_pow(n, p.local_coeff, L(n), p.local_exponent);
}

/// Mode A: N > L_N^64. False for every N < 2^64 (L_N >= 2, so L_N^64 >= 2^64): the outer
/// construction never runs on a representable input and Alg 8 takes its base case.
constexpr bool outer_nontrivial(std::uint64_t n) {
  return outer_nontrivial(n, kModeAParams);
}

/// Mode A: N > 2 * L_N^4; first true at N = 209953 (L_N = 18, 2 * 18^4 = 209952).
constexpr bool local_nontrivial(std::uint64_t n) {
  return local_nontrivial(n, kModeAParams);
}

}  // namespace sssp::paper::model
