// Sec 3.1, Lemma 3.1: stable counting sort on an integer coordinate in [0, k), and stable
// multi-coordinate radix sort (passes from the last coordinate to the first).
// Reference: std::stable_sort.
#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include "paper/primitives/search.hpp"

namespace {

struct Rec {
  std::size_t key;
  std::size_t tag;  // input position, exposes stability
  bool operator==(const Rec&) const = default;
};

struct Rec3 {
  std::array<std::size_t, 3> c;
  std::size_t tag;
  bool operator==(const Rec3&) const = default;
};

}  // namespace

TEST_CASE("counting sort equals std::stable_sort on the key", "[paper][sort][counting]") {
  for (std::uint64_t seed = 0; seed < 60; ++seed) {
    std::mt19937_64 rng(seed);
    const std::size_t p = seed % 7 == 0 ? 0 : 1 + (seed * 37) % 500;
    const std::size_t k = 1 + (seed * 13) % 40;  // includes k = 1 (all keys equal)
    std::uniform_int_distribution<std::size_t> pick(0, k - 1);
    std::vector<Rec> recs;
    for (std::size_t i = 0; i < p; ++i) recs.push_back({pick(rng), i});
    auto expect = recs;
    std::stable_sort(expect.begin(), expect.end(),
                     [](const Rec& a, const Rec& b) { return a.key < b.key; });
    const auto got = sssp::paper::counting_sort_stable(recs, k, [](const Rec& r) { return r.key; });
    REQUIRE(got == expect);
  }
}

TEST_CASE("counting sort keeps input order for all-equal keys", "[paper][sort][counting]") {
  std::vector<Rec> recs;
  for (std::size_t i = 0; i < 100; ++i) recs.push_back({5, i});
  REQUIRE(sssp::paper::counting_sort_stable(recs, 6, [](const Rec& r) { return r.key; }) == recs);
}

TEST_CASE("counting sort rejects a key outside [0, k)", "[paper][sort][counting]") {
  const std::vector<Rec> recs{{0, 0}, {3, 1}};
  REQUIRE_THROWS_AS(sssp::paper::counting_sort_stable(recs, 3, [](const Rec& r) { return r.key; }),
                    std::out_of_range);
  const std::vector<Rec> none;
  REQUIRE(sssp::paper::counting_sort_stable(none, 0, [](const Rec& r) { return r.key; }).empty());
}

TEST_CASE("radix sort is lexicographic and stable on several coordinates", "[paper][sort][radix]") {
  for (std::uint64_t seed = 0; seed < 40; ++seed) {
    std::mt19937_64 rng(500 + seed);
    const std::size_t p = seed % 9 == 0 ? 0 : 1 + (seed * 53) % 400;
    const std::array<std::size_t, 3> ranges{1 + seed % 5, 1 + (seed * 3) % 7, 1 + seed % 2};
    std::vector<Rec3> recs;
    for (std::size_t i = 0; i < p; ++i) {
      Rec3 r{{}, i};
      for (std::size_t j = 0; j < 3; ++j) {
        r.c[j] = std::uniform_int_distribution<std::size_t>(0, ranges[j] - 1)(rng);
      }
      recs.push_back(r);
    }
    auto expect = recs;
    std::stable_sort(expect.begin(), expect.end(),
                     [](const Rec3& a, const Rec3& b) { return a.c < b.c; });
    const auto got =
        sssp::paper::radix_sort_stable(recs, std::span<const std::size_t>(ranges),
                                       [](const Rec3& r, std::size_t j) { return r.c[j]; });
    REQUIRE(got == expect);
  }
}

TEST_CASE("radix sort with no coordinates keeps the input order", "[paper][sort][radix]") {
  std::vector<Rec3> recs{{{2, 0, 0}, 0}, {{1, 0, 0}, 1}};
  const auto got = sssp::paper::radix_sort_stable(
      recs, std::span<const std::size_t>{}, [](const Rec3& r, std::size_t j) { return r.c[j]; });
  REQUIRE(got == recs);
}

TEST_CASE("radix sort groups undirected vertex pairs (min, max)", "[paper][sort][radix]") {
  // The use in Alg 8 l.5 / Alg 4 l.4: bring parallel edges of a vertex pair together.
  struct E {
    std::size_t u, v, id;
  };
  const std::vector<E> edges{{3, 1, 0}, {0, 2, 1}, {1, 3, 2}, {2, 0, 3}, {1, 1, 4}, {0, 1, 5}};
  const std::array<std::size_t, 2> ranges{4, 4};
  const auto got = sssp::paper::radix_sort_stable(
      edges, std::span<const std::size_t>(ranges),
      [](const E& e, std::size_t j) { return j == 0 ? std::min(e.u, e.v) : std::max(e.u, e.v); });
  std::vector<std::size_t> ids;
  for (const auto& e : got) ids.push_back(e.id);
  REQUIRE(ids == std::vector<std::size_t>{5, 1, 3, 4, 0, 2});
}
