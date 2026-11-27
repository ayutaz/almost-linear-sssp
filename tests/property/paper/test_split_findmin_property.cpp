// Property tests: counted work of the split-findmin structures stays within the bounds proved
// in Sec 3.5 (Lemma 3.10 for the dyadic representation).
//
// The bounds are instantiated with the constants that the proofs give:
//  Lemma 3.10 (k records, K = 2^D >= k the padded size, D = ceil(log2 k)):
//   - every rescan moves a record to a strict descendant of its node, so rescans <= k * D;
//   - each node is exposed at most once and the root never is: exposed nodes <= 2K - 2;
//   - a live interval has at most two maximal nodes per depth, so relabelling the two new lists
//     visits at most 2 * 2(D+1) nodes per split, and the cut path has at most D steps.
// Fixed seeds; several seeds and several split orders per size.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <vector>

#include "paper/primitives/split_findmin.hpp"

#if SSSP_DIAGNOSTICS

namespace {

using sssp::paper::sfm_ceil_log2;
using sssp::paper::SplitFindminCounters;

enum class Order { kRandom, kAscending, kDescending, kMiddleOut };

std::vector<std::uint32_t> split_order(std::uint32_t k, Order order, std::mt19937_64& rng) {
  std::vector<std::uint32_t> pos;
  for (std::uint32_t p = 1; p < k; ++p) pos.push_back(p);
  switch (order) {
    case Order::kRandom:
      std::shuffle(pos.begin(), pos.end(), rng);
      break;
    case Order::kAscending:
      break;
    case Order::kDescending:
      std::reverse(pos.begin(), pos.end());
      break;
    case Order::kMiddleOut:
      std::sort(pos.begin(), pos.end(), [k](std::uint32_t a, std::uint32_t b) {
        const auto da = a > k / 2 ? a - k / 2 : k / 2 - a;
        const auto db = b > k / 2 ? b - k / 2 : k / 2 - b;
        return da != db ? da < db : a < b;
      });
      break;
  }
  return pos;
}

/// Performs every split (plus some rejected ones and interleaved decreases and queries) and
/// calls per_split(before, after) around each successful split.
template <class S, class F>
SplitFindminCounters split_everything(S& s, std::uint32_t k, Order order, std::uint64_t seed,
                                      F&& per_split) {
  using Key = typename S::key_type;
  std::mt19937_64 rng(seed);
  for (const std::uint32_t pos : split_order(k, order, rng)) {
    const auto before = s.counters();
    REQUIRE(s.split(pos));
    per_split(before, s.counters());
    REQUIRE_FALSE(s.split(pos));
    const std::uint32_t i = std::uniform_int_distribution<std::uint32_t>(0, k - 1)(rng);
    const Key cur = s.key(i);
    const Key next = static_cast<Key>(std::uniform_int_distribution<std::uint64_t>(0, 1000)(rng));
    if (next < cur) s.decrease(i, next);
    (void)s.query(std::uniform_int_distribution<std::uint32_t>(0, k - 1)(rng));
  }
  return s.counters();
}

const std::vector<std::uint32_t> kSizes{1, 2, 3, 5, 17, 64, 100, 1000, 4096, 5000};
const std::vector<Order> kOrders{Order::kRandom, Order::kAscending, Order::kDescending,
                                 Order::kMiddleOut};

}  // namespace

TEST_CASE("split-findmin property: dyadic work stays within the Lemma 3.10 bounds",
          "[split_findmin][property]") {
  for (const std::uint32_t k : kSizes) {
    const std::uint64_t depth = sfm_ceil_log2(k);
    const std::uint64_t padded = std::uint64_t{1} << depth;
    for (const Order order : kOrders) {
      for (std::uint64_t seed = 1; seed <= 3; ++seed) {
        INFO("k=" << k << " order=" << static_cast<int>(order) << " seed=" << seed);
        std::vector<std::uint64_t> keys(k);
        std::mt19937_64 krng(seed * 7919 + k);
        for (auto& x : keys) x = std::uniform_int_distribution<std::uint64_t>(0, 5000)(krng);
        sssp::paper::DyadicSplitFindmin<std::uint64_t> s{std::span<const std::uint64_t>(keys)};
        const auto c = split_everything(
            s, k, order, seed, [&](const SplitFindminCounters& a, const SplitFindminCounters& b) {
              REQUIRE(b.list_node_visits - a.list_node_visits <= 4 * (depth + 1));
              REQUIRE(b.cut_path_steps - a.cut_path_steps <= depth);
              REQUIRE(b.rescans - a.rescans <= k);
            });
        REQUIRE(c.splits == k - 1);
        REQUIRE(c.rejected_splits == k - 1);
        REQUIRE(c.init_scans == k);
        REQUIRE(c.rescans <= std::uint64_t{k} * depth);
        REQUIRE(c.exposed_nodes <= 2 * padded - 2);
        // After all k-1 splits every record sits alone in a node, so every record of a
        // multi-record array was rescanned at least once.
        if (k > 1) REQUIRE(c.rescans >= k);
      }
    }
  }
}

namespace {

// Rescan bounds of the blocked structures, assembled from the Lemma 3.10 bound of every dyadic
// component exactly as in the proof of Lemma 3.11:
//  one level:  upper dyadic on k' = ceil(k/B) proxies + a dyadic structure per block;
//  two levels: upper dyadic on k' proxies + a one-level structure (block size B') per block.
std::uint64_t dyadic_rescan_bound(std::uint64_t k) {
  return k * sfm_ceil_log2(k);
}

std::uint64_t dyadic_exposed_bound(std::uint64_t k) {
  return 2 * (std::uint64_t{1} << sfm_ceil_log2(k)) - 2;
}

template <class F>
std::uint64_t blocked_sum(std::uint64_t k, std::uint64_t block, F&& lower_bound,
                          std::uint64_t (*upper_bound)(std::uint64_t)) {
  const std::uint64_t blocks = (k + block - 1) / block;
  std::uint64_t total = upper_bound(blocks);
  for (std::uint64_t b = 0; b < blocks; ++b) {
    total += lower_bound(std::min(block, k - b * block));
  }
  return total;
}

std::uint64_t one_level_rescan_bound(std::uint64_t k, std::uint64_t block) {
  return blocked_sum(k, block, dyadic_rescan_bound, dyadic_rescan_bound);
}

std::uint64_t two_level_rescan_bound(std::uint64_t k) {
  const std::uint64_t outer = sssp::paper::sfm_block_size(k);
  const std::uint64_t inner = sssp::paper::sfm_block_size(outer);
  return blocked_sum(
      k, outer, [inner](std::uint64_t len) { return one_level_rescan_bound(len, inner); },
      dyadic_rescan_bound);
}

std::uint64_t one_level_exposed_bound(std::uint64_t k, std::uint64_t block) {
  return blocked_sum(k, block, dyadic_exposed_bound, dyadic_exposed_bound);
}

std::uint64_t two_level_exposed_bound(std::uint64_t k) {
  const std::uint64_t outer = sssp::paper::sfm_block_size(k);
  const std::uint64_t inner = sssp::paper::sfm_block_size(outer);
  return blocked_sum(
      k, outer, [inner](std::uint64_t len) { return one_level_exposed_bound(len, inner); },
      dyadic_exposed_bound);
}

template <class S>
void check_blocked(std::uint64_t (*rescan_bound)(std::uint64_t),
                   std::uint64_t (*exposed_bound)(std::uint64_t),
                   std::uint64_t (*closed_form)(std::uint64_t)) {
  for (const std::uint32_t k : kSizes) {
    for (const Order order : kOrders) {
      for (std::uint64_t seed = 1; seed <= 3; ++seed) {
        INFO("k=" << k << " order=" << static_cast<int>(order) << " seed=" << seed);
        std::vector<double> keys(k);
        std::mt19937_64 krng(seed * 104729 + k);
        for (auto& x : keys) {
          x = static_cast<double>(std::uniform_int_distribution<int>(0, 5000)(krng));
        }
        S s{std::span<const double>(keys)};
        const auto c = split_everything(
            s, k, order, seed, [](const SplitFindminCounters& a, const SplitFindminCounters& b) {
              // Proof of Lemma 3.11: a bounded number of calls into the upper and lower
              // structures per operation (at most one lower split and two upper splits).
              REQUIRE(b.lower_splits - a.lower_splits <= 1);
              REQUIRE(b.upper_splits - a.upper_splits <= 2);
            });
        REQUIRE(c.splits == k - 1);
        REQUIRE(c.rejected_splits == k - 1);
        REQUIRE(c.lower_splits + c.upper_splits <= 3 * c.splits);
        REQUIRE(c.rescans <= rescan_bound(k));
        REQUIRE(c.exposed_nodes <= exposed_bound(k));
        REQUIRE(c.rescans <= closed_form(k));
      }
    }
  }
}

}  // namespace

TEST_CASE("split-findmin property: one blocking step stays within the Lemma 3.11 bounds",
          "[split_findmin][property]") {
  // Closed form of the first blocking step: upper k'ceil(log2 k') <= k + B, lower blocks
  // sum len*ceil(log2 len) <= k*ceil(log2 B), i.e. O(k(1 + log log(k+4))).
  check_blocked<sssp::paper::OneLevelBlockedSplitFindmin<double>>(
      [](std::uint64_t k) { return one_level_rescan_bound(k, sssp::paper::sfm_block_size(k)); },
      [](std::uint64_t k) { return one_level_exposed_bound(k, sssp::paper::sfm_block_size(k)); },
      [](std::uint64_t k) {
        const std::uint64_t b = sssp::paper::sfm_block_size(k);
        return k + b + k * sfm_ceil_log2(b);
      });
}

TEST_CASE("split-findmin property: two blocking steps stay within the Lemma 3.11 bounds",
          "[split_findmin][property]") {
  // Closed form of the second blocking step (B' = max(2, ceil(log2(B+1)))): outer upper
  // <= k + B; per block of length l the inner upper is <= l + B' and the inner lowers are
  // <= l*ceil(log2 B'); in total k(2 + ceil(log2 B')) + B + k'B', i.e. O(k(1 + log log B)).
  check_blocked<sssp::paper::SplitFindmin<double>>(
      two_level_rescan_bound, two_level_exposed_bound, [](std::uint64_t k) {
        const std::uint64_t b = sssp::paper::sfm_block_size(k);
        const std::uint64_t b2 = sssp::paper::sfm_block_size(b);
        const std::uint64_t blocks = (k + b - 1) / b;
        return k * (2 + sfm_ceil_log2(b2)) + b + blocks * b2;
      });
}

TEST_CASE("split-findmin property: blocked work per record shrinks with the blocking depth",
          "[split_findmin][property]") {
  // On k = 2^16 records with every split performed, the dyadic bound is 16 rescans per record;
  // the two-level closed form gives 2 + ceil(log2 B') + small terms (B = 17, B' = 5).
  const std::uint32_t k = 1U << 16;
  sssp::paper::SplitFindmin<std::uint64_t> s(std::size_t{k});
  std::mt19937_64 rng(2026);
  for (const std::uint32_t pos : split_order(k, Order::kRandom, rng)) REQUIRE(s.split(pos));
  const auto c = s.counters();
  REQUIRE(c.rescans <= two_level_rescan_bound(k));
  REQUIRE(c.rescans <= std::uint64_t{k} * 6);
  REQUIRE(c.rescans < dyadic_rescan_bound(k) / 2);
}

#else

TEST_CASE("split-findmin property: counters are compiled out without SSSP_DIAGNOSTICS",
          "[split_findmin][property]") {
  sssp::paper::DyadicSplitFindmin<std::uint64_t> s(std::size_t{16});
  REQUIRE(s.split(5));
  REQUIRE(s.counters().rescans == 0);
  REQUIRE(s.counters().splits == 0);
}

#endif
