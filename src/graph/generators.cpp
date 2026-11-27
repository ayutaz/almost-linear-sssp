// Reproducible graph generators (see include/sssp/generators.hpp for the contract).

// Floating-point contraction (FMA) changes rounding and therefore the generated graphs.
// GCC in ISO mode (-std=c++20) already uses -ffp-contract=off; clang needs this pragma.
#if defined(__clang__)
#pragma clang fp contract(off)
#endif

#include "sssp/generators.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

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

// ------------------------------------------------------------------ weights

namespace {

constexpr double kMaxExactInt = 9007199254740992.0;  // 2^53

void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("invalid weight spec: ") + what);
}

/// Rounds a non-negative real draw to the weight type (nearest integer, ties up, for integers).
template <class W>
W from_real(double x) {
  if constexpr (std::is_integral_v<W>) {
    return static_cast<W>(std::floor(x + 0.5));
  } else {
    return static_cast<W>(x);
  }
}

std::uint64_t int_lo(const WeightSpec& s) {
  return static_cast<std::uint64_t>(std::ceil(s.lo));
}
std::uint64_t int_hi(const WeightSpec& s) {
  return static_cast<std::uint64_t>(std::floor(s.hi));
}

}  // namespace

void validate_weight_spec(const WeightSpec& s) {
  require(std::isfinite(s.lo) && std::isfinite(s.hi), "lo and hi must be finite");
  require(s.lo >= 0.0, "lo must be >= 0");
  require(s.hi <= kMaxExactInt, "hi must be <= 2^53");
  switch (s.dist) {
    case WeightDist::kUniformInt:
      require(s.lo <= s.hi && int_lo(s) <= int_hi(s), "[lo, hi] must contain an integer");
      break;
    case WeightDist::kUniformReal:
      require(s.lo <= s.hi, "lo must be <= hi");
      break;
    case WeightDist::kLogUniform:
      require(s.lo > 0.0 && s.lo <= s.hi, "log_uniform needs 0 < lo <= hi");
      break;
    case WeightDist::kExponential:
      require(std::isfinite(s.mean) && s.mean > 0.0 && s.mean <= 0x1p47,
              "exponential needs 0 < mean <= 2^47");
      break;
    case WeightDist::kConstant:
      break;
    case WeightDist::kZeroMixed:
      require(s.zero_probability >= 0.0 && s.zero_probability <= 1.0,
              "zero_probability must be in [0, 1]");
      require(s.lo <= s.hi && int_lo(s) <= int_hi(s), "[lo, hi] must contain an integer");
      break;
    case WeightDist::kWideRange:
      require(s.wide_log2 >= 0 && s.wide_log2 <= 62, "wide_log2 must be in [0, 62]");
      break;
  }
}

template <class W>
W sample_weight(const WeightSpec& s, Xoshiro256ss& rng) {
  switch (s.dist) {
    case WeightDist::kUniformInt:
      return static_cast<W>(rng.uniform_range(int_lo(s), int_hi(s)));
    case WeightDist::kUniformReal:
      return from_real<W>(s.lo + (s.hi - s.lo) * rng.uniform01());
    case WeightDist::kLogUniform: {
      const double a = det_log(s.lo);
      const double b = det_log(s.hi);
      const double x = det_exp(a + (b - a) * rng.uniform01());
      return from_real<W>(std::clamp(x, s.lo, s.hi));
    }
    case WeightDist::kExponential:
      // 1 - U is in (0, 1], so the logarithm is finite.
      return from_real<W>(-s.mean * det_log(1.0 - rng.uniform01()));
    case WeightDist::kConstant:
      return from_real<W>(s.lo);
    case WeightDist::kZeroMixed:
      if (rng.bernoulli(s.zero_probability)) return W{0};
      if constexpr (std::is_integral_v<W>) {
        return static_cast<W>(rng.uniform_range(int_lo(s), int_hi(s)));
      } else {
        return s.lo + (s.hi - s.lo) * rng.uniform01();
      }
    case WeightDist::kWideRange: {
      const auto k = static_cast<std::uint64_t>(s.wide_log2);
      const auto e = rng.uniform_range(0, k);
      if (e == k) return static_cast<W>(std::uint64_t{1} << k);
      if constexpr (std::is_integral_v<W>) {
        const std::uint64_t base = std::uint64_t{1} << e;
        return static_cast<W>(rng.uniform_range(base, 2 * base - 1));
      } else {
        return std::ldexp(1.0 + rng.uniform01(), static_cast<int>(e));
      }
    }
  }
  return W{0};  // unreachable
}

template std::uint64_t sample_weight<std::uint64_t>(const WeightSpec&, Xoshiro256ss&);
template double sample_weight<double>(const WeightSpec&, Xoshiro256ss&);

namespace {
constexpr std::array<std::pair<WeightDist, std::string_view>, 7> kWeightNames{{
    {WeightDist::kUniformInt, "uniform_int"},
    {WeightDist::kUniformReal, "uniform_real"},
    {WeightDist::kLogUniform, "log_uniform"},
    {WeightDist::kExponential, "exponential"},
    {WeightDist::kConstant, "constant"},
    {WeightDist::kZeroMixed, "zero_mixed"},
    {WeightDist::kWideRange, "wide_range"},
}};
}  // namespace

std::string_view to_string(WeightDist d) noexcept {
  for (const auto& [k, name] : kWeightNames) {
    if (k == d) return name;
  }
  return "unknown";
}

WeightDist parse_weight_dist(std::string_view name) {
  for (const auto& [k, n] : kWeightNames) {
    if (n == name) return k;
  }
  throw std::invalid_argument("unknown weight distribution: " + std::string(name));
}

}  // namespace sssp::gen
