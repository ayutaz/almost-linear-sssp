// Reproducible graph generators (see include/sssp/generators.hpp for the contract).

// Floating-point contraction (FMA) changes rounding and therefore the generated graphs.
// GCC in ISO mode (-std=c++20) already uses -ffp-contract=off; clang needs this pragma.
#if defined(__clang__)
#pragma clang fp contract(off)
#endif

#include "sssp/generators.hpp"

#include <cmath>

namespace sssp::gen {

namespace {
constexpr double kLn2 = 0.693147180559945309417232121458176568;
constexpr double kLn2Hi = 0x1.62e42feep-1;        // fdlibm ln2_hi
constexpr double kLn2Lo = 0x1.a39ef35793c76p-33;  // fdlibm ln2_lo
constexpr double kSqrtHalf = 0.707106781186547524400844362104849039;
}  // namespace

double det_log(double x) noexcept {
  // x = m * 2^e with m in [sqrt(1/2), sqrt(2)); frexp is exact.
  int e = 0;
  double m = std::frexp(x, &e);  // m in [0.5, 1)
  if (m < kSqrtHalf) {
    m *= 2.0;
    --e;
  }
  // log m = 2 atanh(s), s = (m-1)/(m+1), |s| <= 0.1716; the series converges geometrically
  // with ratio s^2 <= 0.0295, so 14 terms reach well below 2^-53.
  const double s = (m - 1.0) / (m + 1.0);
  const double s2 = s * s;
  double term = s;
  double sum = 0.0;
  for (int k = 0; k < 14; ++k) {
    sum += term / static_cast<double>(2 * k + 1);
    term *= s2;
  }
  return static_cast<double>(e) * kLn2 + 2.0 * sum;
}

double det_exp(double x) noexcept {
  // x = k ln2 + r with |r| <= ln2/2; exp(r) by Taylor series; ldexp is exact.
  const double kd = std::floor(x / kLn2 + 0.5);
  // Cody–Waite reduction: ln2 = kLn2Hi + kLn2Lo where kLn2Hi has 32 trailing zero bits, so
  // kd * kLn2Hi is exact for |kd| < 2^20 and the reduction loses no precision.
  const double r = (x - kd * kLn2Hi) - kd * kLn2Lo;
  double term = 1.0;
  double sum = 1.0;
  for (int i = 1; i < 24; ++i) {
    term *= r / static_cast<double>(i);
    sum += term;
  }
  return std::ldexp(sum, static_cast<int>(kd));
}

}  // namespace sssp::gen
