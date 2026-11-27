#pragma once
// Deterministic PRNG for tests. SplitMix64 (Steele, Lea, Flood 2014) is used instead of
// <random> engines + distributions because std::uniform_int_distribution is not specified
// bit-for-bit across standard libraries; a fixed seed must give the same graph everywhere.

#include <cstdint>
#include <limits>

namespace sssp::test {

class SplitMix64 {
 public:
  explicit constexpr SplitMix64(std::uint64_t seed) noexcept : state_(seed) {}

  constexpr std::uint64_t next() noexcept {
    std::uint64_t z = (state_ += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }

  /// Uniform integer in [0, bound) by rejection sampling (unbiased). bound must be > 0.
  constexpr std::uint64_t below(std::uint64_t bound) noexcept {
    const std::uint64_t limit = std::numeric_limits<std::uint64_t>::max() -
                                std::numeric_limits<std::uint64_t>::max() % bound;
    std::uint64_t x = next();
    while (x >= limit) x = next();
    return x % bound;
  }

  /// Uniform integer in [lo, hi] (inclusive).
  constexpr std::uint64_t between(std::uint64_t lo, std::uint64_t hi) noexcept {
    if (hi - lo == std::numeric_limits<std::uint64_t>::max()) return next();
    return lo + below(hi - lo + 1);
  }

  /// Uniform double in [0, 1) with 53 random bits.
  constexpr double unit() noexcept {
    return static_cast<double>(next() >> 11) * 0x1.0p-53;  // NOLINT(readability-magic-numbers)
  }

 private:
  std::uint64_t state_;
};

}  // namespace sssp::test
