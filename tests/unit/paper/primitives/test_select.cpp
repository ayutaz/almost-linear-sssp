// Sec 3.1, Lemma 3.1: BFPRT selection (groups of five, median of medians, three-way partition)
// returns the rank-r key in O(p) comparisons. Reference: std::nth_element.
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/primitives/search.hpp"
#include "sssp/result.hpp"

namespace {

template <class T, class Cmp = std::less<>>
T nth_ref(std::vector<T> a, std::size_t r, Cmp cmp = {}) {
  std::nth_element(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(r), a.end(), cmp);
  return a[r];
}

std::vector<std::int64_t> shaped(int shape, std::size_t p, std::mt19937_64& rng) {
  std::vector<std::int64_t> a(p);
  for (std::size_t i = 0; i < p; ++i) {
    const auto si = static_cast<std::int64_t>(i);
    const auto sp = static_cast<std::int64_t>(p);
    switch (shape) {
      case 0:
        a[i] = static_cast<std::int64_t>(rng() % 1000000);
        break;  // random
      case 1:
        a[i] = static_cast<std::int64_t>(rng() % 3);
        break;  // many duplicates
      case 2:
        a[i] = si;
        break;  // sorted
      case 3:
        a[i] = sp - si;
        break;  // reversed
      case 4:
        a[i] = 7;
        break;  // all equal
      default:
        a[i] = std::min(si, sp - si);
        break;  // organ pipe
    }
  }
  return a;
}

}  // namespace

TEST_CASE("BFPRT equals nth_element for every rank on small inputs", "[paper][select]") {
  for (std::uint64_t seed = 0; seed < 30; ++seed) {
    std::mt19937_64 rng(seed);
    for (int shape = 0; shape < 6; ++shape) {
      const std::size_t p = 1 + (seed * 11 + static_cast<std::uint64_t>(shape)) % 80;
      const auto a = shaped(shape, p, rng);
      for (std::size_t r = 0; r < p; ++r) {
        REQUIRE(sssp::paper::select_bfprt(a, r) == nth_ref(a, r));
      }
    }
  }
}

TEST_CASE("BFPRT equals nth_element on large inputs of every shape", "[paper][select]") {
  std::mt19937_64 rng(99);
  for (int shape = 0; shape < 6; ++shape) {
    for (const std::size_t p : {1000u, 4321u, 50000u}) {
      const auto a = shaped(shape, p, rng);
      for (const std::size_t r : {std::size_t{0}, p / 3, (p - 1) / 2, p / 2, p - 1}) {
        REQUIRE(sssp::paper::select_bfprt(a, r) == nth_ref(a, r));
      }
    }
  }
}

TEST_CASE("BFPRT honours a comparator, including ties between distinct records",
          "[paper][select][comparator]") {
  std::mt19937_64 rng(5);
  std::vector<double> d(777);
  for (auto& x : d) x = static_cast<double>(rng() % 100) * 0.5;
  for (const std::size_t r : {std::size_t{0}, std::size_t{100}, std::size_t{776}}) {
    REQUIRE(sssp::paper::select_bfprt(d, r, std::greater<>{}) == nth_ref(d, r, std::greater<>{}));
  }
  // Records equivalent under the comparator but distinguishable: the result must be an input
  // record equivalent to the reference's rank-r record.
  using P = std::pair<int, int>;
  std::vector<P> recs;
  for (int i = 0; i < 500; ++i) recs.emplace_back(static_cast<int>(rng() % 9), i);
  const auto by_first = [](const P& a, const P& b) { return a.first < b.first; };
  for (std::size_t r = 0; r < recs.size(); r += 37) {
    const P got = sssp::paper::select_bfprt(recs, r, by_first);
    REQUIRE(got.first == nth_ref(recs, r, by_first).first);
    REQUIRE(std::find(recs.begin(), recs.end(), got) != recs.end());
  }
}

TEST_CASE("BFPRT boundary cases", "[paper][select]") {
  REQUIRE(sssp::paper::select_bfprt(std::vector<int>{42}, 0) == 42);
  REQUIRE_THROWS_AS(sssp::paper::select_bfprt(std::vector<int>{}, 0), std::out_of_range);
  REQUIRE_THROWS_AS(sssp::paper::select_bfprt(std::vector<int>{1, 2}, 2), std::out_of_range);
  const std::vector<int> a{5, 1, 4};
  REQUIRE(sssp::paper::select_bfprt(a, 0) == 1);
  REQUIRE(sssp::paper::select_bfprt(a, 1) == 4);
  REQUIRE(sssp::paper::select_bfprt(a, 2) == 5);
  REQUIRE(a == std::vector<int>{5, 1, 4});  // the input is not modified
}

#if SSSP_DIAGNOSTICS  // counters are compiled out when SSSP_DIAGNOSTICS=0
TEST_CASE("BFPRT uses a linear number of comparisons", "[paper][select][lemma3.1]") {
  // Bound derived from the algorithm, not measured: one level on n keys costs at most 2n
  // comparisons to sort the groups of five (<= 10 per group) and at most 2n to partition; the
  // two recursive sizes sum to at most 9n/10 + O(1) (Lemma 3.1 proof), so the total is about
  // 4n / (1 - 9/10) = 40n plus additive terms. 80 p leaves room for those terms.
  std::mt19937_64 rng(1234);
  for (int shape = 0; shape < 6; ++shape) {
    for (const std::size_t p : {10u, 100u, 1000u, 10000u, 200000u}) {
      const auto a = shaped(shape, p, rng);
      for (const std::size_t r : {std::size_t{0}, p / 2, p - 1}) {
        sssp::OpCounters ops;
        (void)sssp::paper::select_bfprt(a, r, std::less<>{}, &ops);
        REQUIRE(ops.comparisons > 0);
        REQUIRE(ops.comparisons <= 80 * p);
      }
    }
  }
}
#endif
