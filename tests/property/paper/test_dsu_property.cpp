// Sec 3.3, Lemma 3.6 (proof): the fixed-band charging argument bounds the total number of
// parent-pointer traversals of all finds by 3 per find plus k*b1 + 2k, where
// b3 = 1 + L_k, b2 = ceil(log2(b3+1)), b1 = ceil(log2(b2+1)) and L_k = ceil(log2(k+2)) (Sec 2).
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <vector>

#include "paper/primitives/dsu.hpp"

using sssp::paper::DisjointSets;

// The bound is checked on counted operations, which exist only with SSSP_DIAGNOSTICS=1.
#if SSSP_DIAGNOSTICS
namespace {

std::uint64_t ceil_log2(std::uint64_t x) {
  std::uint64_t r = 0;
  while ((std::uint64_t{1} << r) < x) ++r;
  return r;
}

std::uint64_t band_bound(std::uint64_t k, std::uint64_t finds) {
  const std::uint64_t lk = ceil_log2(k + 2);
  const std::uint64_t b3 = 1 + lk;
  const std::uint64_t b2 = ceil_log2(b3 + 1);
  const std::uint64_t b1 = ceil_log2(b2 + 1);
  return 3 * finds + k * b1 + 2 * k;
}

}  // namespace

TEST_CASE("dsu find traversals obey the Lemma 3.6 band bound on random op sequences",
          "[paper][dsu][property]") {
  for (const std::uint64_t seed : {11u, 12u, 13u, 14u, 15u, 16u, 17u, 18u}) {
    for (const std::size_t k : {std::size_t{1}, std::size_t{3}, std::size_t{100}, std::size_t{4096},
                                std::size_t{50000}}) {
      std::mt19937_64 rng(seed ^ (k << 20));
      std::uniform_int_distribution<std::uint32_t> item(0, static_cast<std::uint32_t>(k - 1));
      DisjointSets ds(k);
      for (std::size_t op = 0; op < 3 * k; ++op) {
        if (rng() % 3 == 0) {
          (void)ds.find(item(rng));
        } else {
          (void)ds.unite(item(rng), item(rng));
        }
      }
      INFO("seed " << seed << " k " << k);
      REQUIRE(ds.find_steps() <= band_bound(k, ds.find_calls()));
    }
  }
}

TEST_CASE("dsu find traversals obey the Lemma 3.6 band bound on deep binomial trees",
          "[paper][dsu][property]") {
  // Linking equal-rank roots builds binomial trees of depth = rank, the deepest chains the
  // rank rule allows; then finds hit the deepest items first (worst order for compression).
  for (const std::uint64_t seed : {21u, 22u, 23u}) {
    const std::size_t k = std::size_t{1} << 14;
    DisjointSets ds(k);
    for (std::size_t step = 1; step < k; step <<= 1) {
      for (std::size_t i = 0; i + step < k; i += 2 * step) {
        REQUIRE(ds.unite(static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i + step)));
      }
    }
    std::vector<std::uint32_t> order(k);
    for (std::size_t i = 0; i < k; ++i) order[i] = static_cast<std::uint32_t>(i);
    std::mt19937_64 rng(seed);
    std::shuffle(order.begin(), order.end(), rng);
    for (const auto x : order) (void)ds.find(x);
    for (const auto x : order) (void)ds.find(x);
    INFO("seed " << seed);
    REQUIRE(ds.num_sets() == 1);
    REQUIRE(ds.find_steps() <= band_bound(k, ds.find_calls()));
  }
}
#endif
