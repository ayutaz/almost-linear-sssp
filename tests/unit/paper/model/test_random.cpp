// Sec 2, randomness (arXiv 2611.34911v1, txt l.307-310): independent uniformly random words; a
// uniform integer in 0..h-1 takes ceil(log2 h) unbiased bits and rejects outcomes >= h; h = 1
// returns zero; each attempt succeeds with probability greater than 1/2.
#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/model/random.hpp"

using sssp::paper::model::random_permutation;
using sssp::paper::model::random_shuffle;
using sssp::paper::model::uniform_below;
using sssp::paper::model::WordRng;

namespace {
constexpr std::uint64_t kMax = std::numeric_limits<std::uint64_t>::max();

// A bit source that replays a fixed script and records every request.
struct ScriptedBits {
  ScriptedBits() = default;
  explicit ScriptedBits(std::vector<std::uint64_t> s) : script(std::move(s)) {}
  std::vector<std::uint64_t> script;
  std::size_t next = 0;
  std::vector<unsigned> requested;
  std::uint64_t next_bits(unsigned c) {
    requested.push_back(c);
    if (next >= script.size()) throw std::logic_error("script exhausted");
    return script[next++];
  }
};

// Counts attempts while forwarding to a real generator.
struct CountingBits {
  WordRng rng;
  std::uint64_t calls = 0;
  std::uint64_t next_bits(unsigned c) {
    ++calls;
    return rng.next_bits(c);
  }
};

// Upper critical value of chi-square with df degrees of freedom at tail probability 1e-6
// (Wilson-Hilferty approximation, z = 4.753). A deterministic seed either passes or fails
// forever, so the tail is chosen small enough that a correct sampler cannot fail by chance.
double chi_square_critical(unsigned df) {
  const double d = df;
  const double z = 4.753;
  const double t = 1.0 - 2.0 / (9.0 * d) + z * std::sqrt(2.0 / (9.0 * d));
  return d * t * t * t;
}

double chi_square(const std::vector<std::uint64_t>& counts, double expected) {
  double s = 0;
  for (const auto c : counts) {
    const double diff = static_cast<double>(c) - expected;
    s += diff * diff / expected;
  }
  return s;
}
}  // namespace

TEST_CASE("SplitMix64 matches the published reference sequence", "[paper][model][random]") {
  // Reference values computed from the published SplitMix64 algorithm (Steele, Lea, Flood 2014;
  // Vigna's splitmix64.c), independently of this implementation.
  WordRng a(0);
  REQUIRE(a.next_word() == 0xe220a8397b1dcdafull);
  REQUIRE(a.next_word() == 0x6e789e6aa1b965f4ull);
  REQUIRE(a.next_word() == 0x06c45d188009454full);
  WordRng b(1234567);
  REQUIRE(b.next_word() == 0x599ed017fb08fc85ull);
  REQUIRE(b.next_word() == 0x2c73f08458540fa5ull);
  REQUIRE(b.next_word() == 0x883ebce5a3f27c77ull);
}

TEST_CASE("generator is deterministic under a fixed seed", "[paper][model][random]") {
  WordRng a(42);
  WordRng b(42);
  WordRng c(43);
  bool differs = false;
  for (int i = 0; i < 1000; ++i) {
    const auto x = a.next_word();
    REQUIRE(x == b.next_word());
    differs = differs || x != c.next_word();
  }
  REQUIRE(differs);
}

TEST_CASE("next_bits returns the c most significant bits of the next word",
          "[paper][model][random]") {
  WordRng twin(7);
  WordRng rng(7);
  REQUIRE(rng.next_bits(0) == 0);  // zero bits: nothing is drawn, the twin stays in step
  for (unsigned c = 1; c <= 64; ++c) {
    const auto w = twin.next_word();
    const auto x = rng.next_bits(c);
    REQUIRE(x == (c == 64 ? w : (w >> (64 - c))));
    if (c < 64) REQUIRE(x < (std::uint64_t{1} << c));
  }
  REQUIRE_THROWS_AS(rng.next_bits(65), std::invalid_argument);
}

TEST_CASE("uniform_below(1) is zero and draws no bits", "[paper][model][random]") {
  ScriptedBits src;
  REQUIRE(uniform_below(src, 1) == 0);
  REQUIRE(src.requested.empty());
  REQUIRE_THROWS_AS(uniform_below(src, 0), std::invalid_argument);
}

TEST_CASE("uniform_below rejects outcomes at least h instead of reducing modulo h",
          "[paper][model][random]") {
  {
    ScriptedBits src{{5, 7, 3}};  // h = 5 needs 3 bits; 5 and 7 are rejected (5 % 5 would be 0)
    REQUIRE(uniform_below(src, 5) == 3);
    REQUIRE(src.requested == std::vector<unsigned>{3, 3, 3});
  }
  {
    ScriptedBits src{{6, 7, 0}};
    REQUIRE(uniform_below(src, 6) == 0);
    REQUIRE(src.requested.size() == 3);
  }
  {
    ScriptedBits src{{15}};  // h = 16 = 2^4: 4 bits, never rejected
    REQUIRE(uniform_below(src, 16) == 15);
    REQUIRE(src.requested == std::vector<unsigned>{4});
  }
  {
    ScriptedBits src{{1}};  // h = 2: one bit
    REQUIRE(uniform_below(src, 2) == 1);
    REQUIRE(src.requested == std::vector<unsigned>{1});
  }
  {
    ScriptedBits src{{kMax, kMax - 1}};  // h = 2^64 - 1: 64 bits, only 2^64 - 1 is rejected
    REQUIRE(uniform_below(src, kMax) == kMax - 1);
    REQUIRE(src.requested == std::vector<unsigned>{64, 64});
  }
}

TEST_CASE("uniform_below stays in range for arbitrary h", "[paper][model][random]") {
  WordRng rng(2024);
  for (int i = 0; i < 20000; ++i) {
    const std::uint64_t h = (rng.next_word() >> (rng.next_word() & 63)) | 1;
    REQUIRE(uniform_below(rng, h) < h);
  }
}

TEST_CASE("each rejection attempt succeeds with probability above one half",
          "[paper][model][random]") {
  // Worst case h = 2^k + 1: the success probability is h / 2^(k+1), just above 1/2, so the mean
  // number of attempts is below 2. 10000 draws give a standard error of about 0.014.
  for (const unsigned k : {1u, 5u, 20u, 40u, 62u}) {
    CountingBits src{WordRng(k)};
    const std::uint64_t h = (std::uint64_t{1} << k) + 1;
    constexpr int kDraws = 10000;
    for (int i = 0; i < kDraws; ++i) REQUIRE(uniform_below(src, h) < h);
    const double mean = static_cast<double>(src.calls) / kDraws;
    REQUIRE(mean < 2.1);
  }
}

TEST_CASE("chi-square on small ranges is consistent with uniformity", "[paper][model][random]") {
  for (const std::uint64_t seed : {1ull, 2ull, 3ull, 4ull, 5ull}) {
    WordRng rng(seed);
    for (const std::uint64_t h : {2ull, 3ull, 5ull, 6ull, 7ull, 10ull, 17ull, 100ull}) {
      const std::uint64_t per_bucket = 500;
      std::vector<std::uint64_t> counts(h, 0);
      for (std::uint64_t i = 0; i < per_bucket * h; ++i) ++counts[uniform_below(rng, h)];
      const double stat = chi_square(counts, static_cast<double>(per_bucket));
      REQUIRE(stat < chi_square_critical(static_cast<unsigned>(h - 1)));
    }
  }
}

TEST_CASE("random_permutation returns a permutation and is deterministic",
          "[paper][model][random]") {
  WordRng rng(99);
  REQUIRE(random_permutation(0, rng).empty());
  REQUIRE(random_permutation(1, rng) == std::vector<std::uint32_t>{0});
  for (const std::size_t n : {2u, 3u, 10u, 1000u}) {
    WordRng a(n);
    WordRng b(n);
    auto p = random_permutation(n, a);
    REQUIRE(p == random_permutation(n, b));
    std::sort(p.begin(), p.end());
    for (std::size_t i = 0; i < n; ++i) REQUIRE(p[i] == i);
  }
}

TEST_CASE("random_permutation is uniform over all orders of four elements",
          "[paper][model][random]") {
  for (const std::uint64_t seed : {11ull, 12ull, 13ull}) {
    WordRng rng(seed);
    std::map<std::vector<std::uint32_t>, std::uint64_t> seen;
    const std::uint64_t per_order = 400;
    for (std::uint64_t i = 0; i < 24 * per_order; ++i) ++seen[random_permutation(4, rng)];
    REQUIRE(seen.size() == 24);
    std::vector<std::uint64_t> counts;
    for (const auto& [perm, c] : seen) counts.push_back(c);
    REQUIRE(chi_square(counts, static_cast<double>(per_order)) < chi_square_critical(23));
  }
}

TEST_CASE("random_shuffle permutes a span in place", "[paper][model][random]") {
  WordRng rng(5);
  std::array<int, 6> a{4, 4, 1, 9, 2, 7};
  auto sorted_before = a;
  std::sort(sorted_before.begin(), sorted_before.end());
  random_shuffle(std::span<int>(a), rng);
  auto sorted_after = a;
  std::sort(sorted_after.begin(), sorted_after.end());
  REQUIRE(sorted_after == sorted_before);
  std::array<int, 0> empty{};
  random_shuffle(std::span<int>(empty), rng);  // no-op
}
