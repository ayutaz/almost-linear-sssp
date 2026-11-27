// Regression: CsrGraph<Counted<double>> accepted negative / NaN weights because
// validate_edges only checked signed and floating-point weight types (reported by the
// verification review, 2026-11-27).
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <stdexcept>

#include "sssp/counted.hpp"
#include "sssp/graph.hpp"

using C = sssp::Counted<double>;

TEST_CASE("Counted<double> graphs reject negative and NaN weights",
          "[graph][counted][regression]") {
  REQUIRE_THROWS_AS(sssp::CsrGraph<C>::from_edges(2, {{0, 1, C(-1.0)}}), std::invalid_argument);
  REQUIRE_THROWS_AS(sssp::CsrGraph<C>::from_edges(2, {{0, 1, C(std::nan(""))}}),
                    std::invalid_argument);
  REQUIRE_NOTHROW(sssp::CsrGraph<C>::from_edges(2, {{0, 1, C(0.0)}}));
}
