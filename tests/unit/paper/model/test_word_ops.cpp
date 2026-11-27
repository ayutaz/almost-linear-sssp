// Sec 2 word operations (arXiv 2611.34911v1, txt l.301-308), decision (b) of Q-51:
//  - "For fixed-size blocks of any array, a scan stores the block index and offset for every
//    member. No integer division is needed for these lookups" (l.305-306): BlockLayout must
//    agree with i / B and i % B on every member, and is built in one scan.
//  - "An integer logarithm on a range of k indices can be tabulated in O(k) time" (l.303-304):
//    the table is built in one step per index and agrees with the word logarithms.
//  - Uniform draws (l.307-308) that read ceil(log2 h) from the table request exactly the same
//    bits as the word version, so every random result is unchanged.
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/model/ilog.hpp"
#include "paper/model/random.hpp"
#include "paper/model/word_ops.hpp"
#include "sssp/result.hpp"

using sssp::paper::model::BlockLayout;
using sssp::paper::model::IlogTable;
using sssp::paper::model::random_permutation;
using sssp::paper::model::random_shuffle;
using sssp::paper::model::uniform_below;
using sssp::paper::model::WordRng;

namespace {

// Records every bit request and forwards it to a real generator.
struct RecordingBits {
  WordRng rng;
  std::vector<unsigned> requested;
  std::uint64_t next_bits(unsigned c) {
    requested.push_back(c);
    return rng.next_bits(c);
  }
};

// Fisher-Yates with the word-model draw of Sec 2 l.307-308 and ceil(log2 h) computed per draw
// by the word logarithm: the behaviour random_shuffle had before the table was introduced.
template <class T>
void reference_shuffle(std::vector<T>& a, WordRng& rng) {
  for (std::size_t i = a.size(); i > 1; --i) {
    const auto j = static_cast<std::size_t>(uniform_below(rng, i));
    std::swap(a[i - 1], a[j]);
  }
}

}  // namespace

TEST_CASE("block layout agrees with integer division on every member", "[paper][model][word]") {
  for (std::uint32_t k = 0; k <= 130; ++k) {
    for (std::uint32_t b = 1; b <= k + 3; ++b) {
      const BlockLayout lay(k, b);
      REQUIRE(lay.size() == k);
      REQUIRE(lay.block_size() == b);
      REQUIRE(lay.num_blocks() == (k + b - 1) / b);
      for (std::uint32_t i = 0; i < k; ++i) {
        REQUIRE(lay.block_of(i) == i / b);
        REQUIRE(lay.offset_of(i) == i % b);
      }
      for (std::uint32_t j = 0; j < lay.num_blocks(); ++j) {
        REQUIRE(lay.block_start(j) == j * b);
        REQUIRE(lay.block_end(j) == std::min(k, (j + 1) * b));
      }
      REQUIRE(lay.block_start(lay.num_blocks()) == k);
    }
  }
}

TEST_CASE("block layout agrees with integer division on large arrays", "[paper][model][word]") {
  for (const std::uint32_t k : {1000u, 4096u, 65535u, 65536u, 100003u}) {
    for (const std::uint32_t b : {2u, 7u, 17u, 32u, 1000u}) {
      const BlockLayout lay(k, b);
      for (std::uint32_t i = 0; i < k; ++i) {
        REQUIRE(lay.block_of(i) == i / b);
        REQUIRE(lay.offset_of(i) == i % b);
      }
    }
  }
}

TEST_CASE("block layout rejects a zero block size", "[paper][model][word]") {
  REQUIRE_THROWS_AS(BlockLayout(5, 0), std::invalid_argument);
  const BlockLayout empty;
  REQUIRE(empty.size() == 0);
  REQUIRE(empty.num_blocks() == 0);
}

TEST_CASE("ilog table agrees with the word logarithms on every index of its range",
          "[paper][model][word]") {
  constexpr std::uint64_t kRange = std::uint64_t{1} << 18;
  const IlogTable lg(kRange);
  for (std::uint64_t i = 1; i <= kRange; ++i) {
    REQUIRE(lg.floor_log2(i) == sssp::paper::model::floor_log2(i));
    REQUIRE(lg.ceil_log2(i) == sssp::paper::model::ceil_log2(i));
  }
}

#if SSSP_DIAGNOSTICS
TEST_CASE("block layout and ilog table are built with one step per member",
          "[paper][model][word]") {
  // Linear construction (Sec 2 l.303-306): the step counters equal the array length.
  for (const std::uint32_t k : {0u, 1u, 2u, 3u, 100u, 4097u}) {
    for (const std::uint32_t b : {1u, 3u, 64u}) {
      REQUIRE(BlockLayout(k, b).build_steps() == k);
    }
    REQUIRE(IlogTable(k).build_steps() == k);
  }
}
#else
TEST_CASE("word-op build counters are compiled out without SSSP_DIAGNOSTICS",
          "[paper][model][word]") {
  REQUIRE(BlockLayout(100, 3).build_steps() == 0);
  REQUIRE(IlogTable(100).build_steps() == 0);
}
#endif

TEST_CASE("uniform draws through the ilog table request the same bits as the word version",
          "[paper][model][word]") {
  constexpr std::uint64_t kRange = 3000;
  const IlogTable lg(kRange);
  for (std::uint64_t h = 1; h <= kRange; ++h) {
    RecordingBits viaTable{WordRng(h), {}};
    RecordingBits viaWord{WordRng(h), {}};
    for (int r = 0; r < 4; ++r) {
      REQUIRE(uniform_below(viaTable, h, lg) == uniform_below(viaWord, h));
    }
    REQUIRE(viaTable.requested == viaWord.requested);
  }
}

TEST_CASE("uniform draws through the ilog table check their range", "[paper][model][word]") {
  const IlogTable lg(10);
  RecordingBits src{WordRng(1), {}};
  REQUIRE_THROWS_AS(uniform_below(src, 0, lg), std::invalid_argument);
  REQUIRE_THROWS_AS(uniform_below(src, 11, lg), std::out_of_range);
  REQUIRE(uniform_below(src, 1, lg) == 0);  // h = 1 draws nothing
  REQUIRE(src.requested.empty());
}

TEST_CASE("random_shuffle and random_permutation are unchanged by the tabulated logarithm",
          "[paper][model][word]") {
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    for (const std::size_t n : {std::size_t{0}, std::size_t{1}, std::size_t{2}, std::size_t{5},
                                std::size_t{64}, std::size_t{65}, std::size_t{1000}}) {
      std::vector<std::uint32_t> expected(n);
      for (std::size_t i = 0; i < n; ++i) expected[i] = static_cast<std::uint32_t>(i);
      WordRng ref(seed);
      reference_shuffle(expected, ref);

      WordRng rng(seed);
      REQUIRE(random_permutation(n, rng) == expected);

      std::vector<std::uint32_t> a(n);
      for (std::size_t i = 0; i < n; ++i) a[i] = static_cast<std::uint32_t>(i);
      WordRng rng2(seed);
      random_shuffle(std::span<std::uint32_t>(a), rng2);
      REQUIRE(a == expected);
    }
  }
}
