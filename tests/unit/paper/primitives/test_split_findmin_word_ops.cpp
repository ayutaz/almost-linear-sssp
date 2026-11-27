// Sec 3.5 Lemma 3.11 with the block lookups of Sec 2 l.305-306 (Q-51 decision (b), resolves
// D-18): the blocked split-findmin finds a record's block and offset through a scan-built
// layout instead of integer division. These tests pin the observable behaviour over every
// block size from 1 to beyond k, including short final blocks, against an independent naive
// model, so that replacing the division by stored lookups cannot change any result.
#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/primitives/split_findmin.hpp"

namespace {

using sssp::paper::IntervalBounds;
using sssp::paper::SfmMin;
using sssp::paper::model::uniform_below;
using sssp::paper::model::WordRng;

/// Linear-scan oracle over explicit interval starts.
template <class Key>
struct Naive {
  std::vector<Key> keys;
  std::vector<std::uint8_t> start;

  explicit Naive(std::vector<Key> k) : keys(std::move(k)), start(keys.size(), 0) { start[0] = 1; }

  [[nodiscard]] IntervalBounds bounds_of(std::uint32_t i) const {
    std::uint32_t lo = i;
    while (start[lo] == 0) --lo;
    std::uint32_t hi = i + 1;
    while (hi < keys.size() && start[hi] == 0) ++hi;
    return {lo, hi};
  }
  bool split(std::uint32_t pos) {
    if (pos == 0 || pos >= keys.size() || start[pos] != 0) return false;
    start[pos] = 1;
    return true;
  }
  [[nodiscard]] SfmMin<Key> query(std::uint32_t i) const {
    const IntervalBounds b = bounds_of(i);
    SfmMin<Key> m{keys[b.lo], b.lo};
    for (std::uint32_t j = b.lo + 1; j < b.hi; ++j) {
      if (keys[j] < m.key) m = {keys[j], j};
    }
    return m;
  }
};

template <class S, class Key>
void require_same(const S& s, const Naive<Key>& n) {
  for (std::uint32_t i = 0; i < n.keys.size(); ++i) {
    REQUIRE(s.key(i) == n.keys[i]);
    REQUIRE(s.query(i) == n.query(i));
    REQUIRE(s.bounds(s.interval_of(i)) == n.bounds_of(i));
  }
}

/// Random splits and decreases on k records with block size b; checks after every operation.
template <class S>
void run_random(std::uint32_t k, std::uint32_t b, std::uint64_t seed) {
  using Key = typename S::key_type;
  WordRng rng(seed);
  std::vector<Key> keys(k);
  for (auto& x : keys) x = static_cast<Key>(uniform_below(rng, 4 * std::uint64_t{k} + 1));
  S s(std::span<const Key>(keys), b);
  Naive<Key> n(keys);
  REQUIRE(s.block_size() == b);
  REQUIRE(s.num_blocks() == (k + b - 1) / b);
  for (std::uint32_t j = 0; j < s.num_blocks(); ++j) {
    const std::uint32_t lo = j * b;
    REQUIRE(s.lower(j).size() == std::min(b, k - lo));
  }
  require_same(s, n);
  for (std::uint32_t step = 0; step < 2 * k; ++step) {
    if (uniform_below(rng, 2) == 0) {
      const auto pos = static_cast<std::uint32_t>(uniform_below(rng, std::uint64_t{k} + 1));
      REQUIRE(s.split(pos) == n.split(pos));
    } else {
      const auto i = static_cast<std::uint32_t>(uniform_below(rng, k));
      const auto key =
          static_cast<Key>(uniform_below(rng, static_cast<std::uint64_t>(n.keys[i]) + 1));
      s.decrease(i, key);
      n.keys[i] = key;
    }
    require_same(s, n);
  }
}

}  // namespace

TEMPLATE_TEST_CASE("blocked split-findmin: every block size gives the naive results",
                   "[split_findmin][word]",
                   (sssp::paper::OneLevelBlockedSplitFindmin<std::uint64_t>),
                   (sssp::paper::OneLevelBlockedSplitFindmin<double>),
                   (sssp::paper::SplitFindmin<std::uint64_t>),
                   (sssp::paper::SplitFindmin<double>)) {
  for (std::uint32_t k = 1; k <= 40; ++k) {
    for (std::uint32_t b = 1; b <= k + 2; ++b) {
      run_random<TestType>(k, b, 1000 * std::uint64_t{k} + b);
    }
  }
}

TEMPLATE_TEST_CASE("blocked split-findmin: default block sizes on larger arrays",
                   "[split_findmin][word]", (sssp::paper::SplitFindmin<std::uint64_t>),
                   (sssp::paper::OneLevelBlockedSplitFindmin<std::uint64_t>)) {
  for (const std::uint32_t k : {63u, 64u, 65u, 127u, 128u, 129u, 300u}) {
    run_random<TestType>(k, sssp::paper::sfm_block_size(k), k);
  }
}
