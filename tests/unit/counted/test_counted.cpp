// sssp::Counted<W> counts comparisons and additions on weights (comparison-addition model).

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <vector>

#include "sssp/counted.hpp"
#include "sssp/graph.hpp"

using sssp::Counted;
using sssp::counted_stats;
using sssp::CountedScope;
using sssp::reset_counted_stats;

using U64 = std::uint64_t;

// Only the operations of the comparison-addition model exist: code that compiles with
// Counted<W> provably performs no subtraction, multiplication, division or rounding on weights.
template <class T>
concept HasMinus = requires(T a) { a - a; };
template <class T>
concept HasTimes = requires(T a) { a * a; };
template <class T>
concept HasDivide = requires(T a) { a / a; };
static_assert(!HasMinus<Counted<double>> && !HasMinus<Counted<U64>>);
static_assert(!HasTimes<Counted<double>> && !HasTimes<Counted<U64>>);
static_assert(!HasDivide<Counted<double>> && !HasDivide<Counted<U64>>);
static_assert(sizeof(Counted<double>) == sizeof(double));
static_assert(std::is_trivially_copyable_v<Counted<U64>>);
static_assert(!std::is_convertible_v<Counted<double>, double>);  // no silent escape

TEMPLATE_TEST_CASE("Counted counts each addition and comparison", "[counted]", U64, double) {
  using C = Counted<TestType>;
  reset_counted_stats();
  const C a{TestType{2}};
  const C b{TestType{3}};
  REQUIRE(counted_stats().additions == 0);
  REQUIRE(counted_stats().comparisons == 0);

  const C s = a + b;
  REQUIRE(s.value() == TestType{5});
  REQUIRE(counted_stats().additions == 1);

  C acc = a;
  acc += b;
  REQUIRE(acc.value() == TestType{5});
  REQUIRE(counted_stats().additions == 2);

  REQUIRE(a < b);
  REQUIRE(a <= b);
  REQUIRE(b > a);
  REQUIRE(b >= a);
  REQUIRE(a == a);
  REQUIRE(a != b);
  REQUIRE(counted_stats().comparisons == 6);
  REQUIRE(counted_stats().additions == 2);

  // Mixed with a raw weight (implicit conversion of the raw value) still counts.
  REQUIRE(a < TestType{9});
  REQUIRE(TestType{1} < a);
  REQUIRE(counted_stats().comparisons == 8);
}

TEST_CASE("Counted copies, constructions and value() are free", "[counted]") {
  reset_counted_stats();
  Counted<double> a{1.5};
  Counted<double> b = a;
  a = b;
  REQUIRE(a.value() == 1.5);
  REQUIRE(Counted<double>{}.value() == 0.0);
  REQUIRE(counted_stats().additions == 0);
  REQUIRE(counted_stats().comparisons == 0);
}

TEST_CASE("CountedScope measures the delta of a region", "[counted]") {
  reset_counted_stats();
  const Counted<U64> x{1};
  (void)(x + x);
  const CountedScope scope;
  (void)(x + x);
  (void)(x + x);
  (void)(x < x);
  REQUIRE(scope.delta().additions == 2);
  REQUIRE(scope.delta().comparisons == 1);
  REQUIRE(counted_stats().additions == 3);
}

TEST_CASE("Counted statistics are per thread", "[counted]") {
  reset_counted_stats();
  const Counted<U64> x{1};
  (void)(x + x);
  std::uint64_t other_additions = 0;
  std::thread t([&other_additions] {
    const Counted<U64> y{2};
    (void)(y + y);
    (void)(y + y);
    (void)(y + y);
    other_additions = counted_stats().additions;
  });
  t.join();
  REQUIRE(other_additions == 3);
  REQUIRE(counted_stats().additions == 1);
}

TEMPLATE_TEST_CASE("WeightTraits support Counted weights", "[counted]", U64, double) {
  using C = Counted<TestType>;
  using T = sssp::WeightTraits<C>;
  using Base = sssp::WeightTraits<TestType>;
  STATIC_REQUIRE(T::is_exact == Base::is_exact);
  REQUIRE(T::infinity().value() == Base::infinity());
  REQUIRE(T::zero().value() == TestType{0});
  reset_counted_stats();
  REQUIRE(T::checked_add(C{TestType{2}}, C{TestType{3}}).value() == TestType{5});
  REQUIRE(counted_stats().additions == 1);  // one model addition; the overflow guard is free
  REQUIRE(counted_stats().comparisons == 0);
  if constexpr (Base::is_exact) {
    REQUIRE_THROWS_AS(T::checked_add(C{Base::infinity() - 1}, C{TestType{1}}), std::overflow_error);
  }
}

namespace {

// A tiny generic Bellman-Ford written only against the comparison-addition interface.
template <class W>
std::vector<W> generic_relax(const sssp::CsrGraph<W>& g, sssp::Vertex s) {
  using T = sssp::WeightTraits<W>;
  std::vector<W> d(g.num_vertices(), T::infinity());
  d[s] = T::zero();
  for (std::size_t round = 0; round + 1 < g.num_vertices(); ++round) {
    for (sssp::Vertex u = 0; u < g.num_vertices(); ++u) {
      if (d[u] == T::infinity()) continue;
      for (const auto& a : g.arcs(u)) {
        const W cand = d[u] + a.weight;
        if (cand < d[a.to]) d[a.to] = cand;
      }
    }
  }
  return d;
}

}  // namespace

TEST_CASE("Counted weights instrument a generic SSSP routine", "[counted]") {
  using C = Counted<U64>;
  // Path 0 -1- 1 -1- 2: 2 rounds; per round vertices are scanned in order 0, 1, 2.
  const auto g = sssp::CsrGraph<C>::from_edges(3, {{0, 1, C{1}}, {1, 2, C{1}}});
  reset_counted_stats();
  const auto d = generic_relax(g, 0);
  REQUIRE(d[0].value() == 0);
  REQUIRE(d[1].value() == 1);
  REQUIRE(d[2].value() == 2);
  // Round 1: u=0 finite (1 arc), u=1 finite after relaxation (2 arcs), u=2 finite (1 arc).
  // Round 2: the same 4 arcs. Each arc costs 1 addition + 1 comparison; each vertex costs one
  // infinity test (comparison): 2 rounds * 3 vertices.
  REQUIRE(counted_stats().additions == 8);
  REQUIRE(counted_stats().comparisons == 8 + 6);
}
