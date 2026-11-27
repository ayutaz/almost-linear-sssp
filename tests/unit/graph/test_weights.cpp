#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "sssp/generators.hpp"

using sssp::gen::sample_weight;
using sssp::gen::WeightDist;
using sssp::gen::WeightSpec;
using sssp::gen::Xoshiro256ss;

namespace {

template <class W>
std::vector<W> draw(const WeightSpec& spec, std::uint64_t seed, int count) {
  sssp::gen::validate_weight_spec(spec);
  Xoshiro256ss rng(seed);
  std::vector<W> out;
  for (int i = 0; i < count; ++i) out.push_back(sample_weight<W>(spec, rng));
  return out;
}

template <class W>
bool is_integral_value(W w) {
  if constexpr (std::is_integral_v<W>) {
    return true;
  } else {
    return std::floor(w) == w;
  }
}

}  // namespace

TEMPLATE_TEST_CASE("uniform_int weights are integers within [lo, hi] and hit both ends",
                   "[graph][weights]", std::uint64_t, double) {
  WeightSpec spec{.dist = WeightDist::kUniformInt, .lo = 2, .hi = 6};
  const auto ws = draw<TestType>(spec, 1, 2000);
  std::set<TestType> seen(ws.begin(), ws.end());
  REQUIRE(seen == std::set<TestType>{2, 3, 4, 5, 6});
  for (auto w : ws) REQUIRE(is_integral_value(w));
}

TEST_CASE("uniform_real weights are in [lo, hi) with mean near the midpoint", "[graph][weights]") {
  WeightSpec spec{.dist = WeightDist::kUniformReal, .lo = 1.0, .hi = 3.0};
  const auto ws = draw<double>(spec, 2, 20000);
  double sum = 0;
  for (double w : ws) {
    REQUIRE(w >= 1.0);
    REQUIRE(w < 3.0);
    sum += w;
  }
  REQUIRE(sum / static_cast<double>(ws.size()) == Catch::Approx(2.0).margin(0.02));
  REQUIRE(std::set<double>(ws.begin(), ws.end()).size() > 19000);  // genuinely real-valued
}

TEMPLATE_TEST_CASE("log_uniform weights span orders of magnitude evenly in log scale",
                   "[graph][weights]", std::uint64_t, double) {
  WeightSpec spec{.dist = WeightDist::kLogUniform, .lo = 1.0, .hi = 1e6};
  const auto ws = draw<TestType>(spec, 3, 20000);
  int below_1000 = 0;
  for (auto w : ws) {
    REQUIRE(static_cast<double>(w) >= 1.0);
    REQUIRE(static_cast<double>(w) <= 1e6);
    if (static_cast<double>(w) < 1000.0) ++below_1000;
  }
  // log-uniform: P(w < 10^3) = 1/2
  REQUIRE(below_1000 == Catch::Approx(10000).margin(400));
}

TEST_CASE("exponential weights are non-negative with the requested mean", "[graph][weights]") {
  WeightSpec spec{.dist = WeightDist::kExponential, .mean = 5.0};
  const auto ws = draw<double>(spec, 4, 40000);
  double sum = 0;
  for (double w : ws) {
    REQUIRE(w >= 0.0);
    REQUIRE(std::isfinite(w));
    sum += w;
  }
  REQUIRE(sum / static_cast<double>(ws.size()) == Catch::Approx(5.0).margin(0.1));
}

TEMPLATE_TEST_CASE("constant weights all equal lo", "[graph][weights]", std::uint64_t, double) {
  WeightSpec spec{.dist = WeightDist::kConstant, .lo = 7};
  for (auto w : draw<TestType>(spec, 5, 100)) REQUIRE(w == TestType{7});
  WeightSpec zero{.dist = WeightDist::kConstant, .lo = 0};
  for (auto w : draw<TestType>(zero, 5, 100)) REQUIRE(w == TestType{0});
}

TEMPLATE_TEST_CASE("zero_mixed weights contain zeros at the requested rate", "[graph][weights]",
                   std::uint64_t, double) {
  WeightSpec spec{.dist = WeightDist::kZeroMixed, .lo = 1, .hi = 10, .zero_probability = 0.3};
  const auto ws = draw<TestType>(spec, 6, 20000);
  int zeros = 0;
  for (auto w : ws) {
    if (w == TestType{0}) {
      ++zeros;
    } else {
      REQUIRE(static_cast<double>(w) >= 1.0);
      REQUIRE(static_cast<double>(w) <= 10.0);
    }
  }
  REQUIRE(zeros == Catch::Approx(6000).margin(300));
}

TEMPLATE_TEST_CASE("wide_range weights cover [1, 2^K] across all binary magnitudes",
                   "[graph][weights]", std::uint64_t, double) {
  WeightSpec spec{.dist = WeightDist::kWideRange, .wide_log2 = 50};
  const auto ws = draw<TestType>(spec, 7, 20000);
  std::set<int> exponents;
  const double top = std::ldexp(1.0, 50);
  for (auto w : ws) {
    const double d = static_cast<double>(w);
    REQUIRE(d >= 1.0);
    REQUIRE(d <= top);
    exponents.insert(std::ilogb(d));
  }
  REQUIRE(exponents.size() >= 50);  // binades 2^0 .. 2^49 all present
  const auto [mn, mx] = std::minmax_element(ws.begin(), ws.end());
  REQUIRE(static_cast<double>(*mx) / static_cast<double>(*mn) > std::ldexp(1.0, 45));
}

TEST_CASE("weight specs with invalid parameters are rejected", "[graph][weights][validation]") {
  using sssp::gen::validate_weight_spec;
  REQUIRE_THROWS_AS(validate_weight_spec({.dist = WeightDist::kUniformInt, .lo = 5, .hi = 2}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(validate_weight_spec({.dist = WeightDist::kUniformReal, .lo = -1, .hi = 2}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(validate_weight_spec({.dist = WeightDist::kLogUniform, .lo = 0, .hi = 2}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(validate_weight_spec({.dist = WeightDist::kExponential, .mean = -1}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(validate_weight_spec({.dist = WeightDist::kZeroMixed, .zero_probability = 2}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(validate_weight_spec({.dist = WeightDist::kWideRange, .wide_log2 = 63}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(
      validate_weight_spec({.dist = WeightDist::kUniformReal, .lo = 0, .hi = std::nan("")}),
      std::invalid_argument);
}

TEST_CASE("weight distribution names round-trip", "[graph][weights]") {
  for (auto d : {WeightDist::kUniformInt, WeightDist::kUniformReal, WeightDist::kLogUniform,
                 WeightDist::kExponential, WeightDist::kConstant, WeightDist::kZeroMixed,
                 WeightDist::kWideRange}) {
    REQUIRE(sssp::gen::parse_weight_dist(sssp::gen::to_string(d)) == d);
  }
  REQUIRE(sssp::gen::parse_weight_dist("uniform_int") == WeightDist::kUniformInt);
  REQUIRE_THROWS_AS(sssp::gen::parse_weight_dist("gaussian"), std::invalid_argument);
}
