#pragma once
// Saturating length arithmetic for the classification tests and sentinels of Sec 5 (open
// question Q-38): N^3 sigma (Sec 5.2 l.1026-1028), K^2 r (l.1032), a K^j (l.1035-1036), and
// later N^c sigma and N^16 sigma (Sec 5.4). These products are used only to classify lengths
// and as sentinels, never as distances, so a result that does not fit the type is replaced by
// WeightTraits<W>::infinity(). Comparing a finite length x <= product is then exact: if the true
// product overflows, it exceeds every finite length of the type.
//
// Only additions and doublings are performed (Sec 2 l.297: multiplication by 2^j takes j
// doublings; multiplication by an integer uses its binary expansion, Sec 5.2 l.1027-1028).
// For double the doublings are exact and the integer multiplication rounds each addition.

#include <cstdint>
#include <type_traits>

#include "sssp/graph.hpp"

namespace sssp::paper::sat {

/// a + b, or infinity if the sum does not fit (integers) or overflows (floating point).
template <class W>
W add(W a, W b) {
  if constexpr (std::is_integral_v<W>) {
    constexpr W kInf = WeightTraits<W>::infinity();
    if (a >= kInf - b) return kInf;
    return static_cast<W>(a + b);
  } else {
    return a + b;  // IEEE overflow gives +inf
  }
}

/// a * 2^k by k doublings.
template <class W>
W mul_pow2(W a, unsigned k) {
  for (unsigned i = 0; i < k && a != WeightTraits<W>::infinity(); ++i) a = add(a, a);
  return a;
}

/// a * n by the binary expansion of the integer n (doubling and adding, most significant bit
/// first).
template <class W>
W mul_int(W a, std::uint64_t n) {
  W r{0};
  for (int bit = 63; bit >= 0; --bit) {
    r = add(r, r);
    if (((n >> bit) & 1U) != 0) r = add(r, a);
  }
  return r;
}

/// a * n^e by e integer multiplications.
template <class W>
W mul_int_pow(W a, std::uint64_t n, unsigned e) {
  for (unsigned i = 0; i < e; ++i) a = mul_int(a, n);
  return a;
}

}  // namespace sssp::paper::sat
