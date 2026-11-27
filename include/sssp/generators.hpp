#pragma once
// Reproducible graph generators.
//
// Reproducibility contract: every generator is a pure function of (parameters, seed) and produces
// bit-identical output on every platform and standard library. To guarantee that we
//   * use our own PRNGs (splitmix64 for seeding, xoshiro256** for streams) instead of <random>
//     engines + std::*_distribution, whose algorithms are implementation-defined;
//   * derive integers with exact integer arithmetic (rejection sampling, no modulo bias);
//   * derive reals only from IEEE-754 basic operations (+ - * / sqrt, all correctly rounded),
//     with log/exp implemented in this library (det_log / det_exp) because libm's
//     transcendental functions are not correctly rounded and differ between platforms;
//   * never iterate over unordered containers when producing output.
// FMA contraction is disabled in the implementation file (it would change rounding).

#include <array>
#include <cstddef>
#include <cstdint>

namespace sssp::gen {

// ------------------------------------------------------------------ PRNG

/// splitmix64 (Steele, Lea, Flood 2014; reference code by S. Vigna). Used to expand a 64-bit
/// seed into xoshiro state and to derive independent sub-stream seeds.
class SplitMix64 {
 public:
  explicit constexpr SplitMix64(std::uint64_t seed) noexcept : state_(seed) {}
  constexpr std::uint64_t next() noexcept {
    std::uint64_t z = (state_ += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }

 private:
  std::uint64_t state_;
};

/// xoshiro256** 1.0 (Blackman, Vigna 2018).
class Xoshiro256ss {
 public:
  using State = std::array<std::uint64_t, 4>;

  /// Seeds the 256-bit state with four splitmix64 outputs (the reference recommendation).
  explicit constexpr Xoshiro256ss(std::uint64_t seed) noexcept : s_{} {
    SplitMix64 sm(seed);
    for (auto& w : s_) w = sm.next();
  }
  /// Raw state (must not be all zero); used for known-answer tests.
  static constexpr Xoshiro256ss from_state(State s) noexcept {
    Xoshiro256ss r(0);
    r.s_ = s;
    return r;
  }

  constexpr std::uint64_t next() noexcept {
    const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }

  /// Uniform integer in [0, bound), bound >= 1. Unbiased: rejects the low (2^64 mod bound)
  /// values so that every residue class is equally likely.
  constexpr std::uint64_t uniform_below(std::uint64_t bound) noexcept {
    const std::uint64_t threshold = (0 - bound) % bound;  // == 2^64 mod bound
    for (;;) {
      const std::uint64_t r = next();
      if (r >= threshold) return r % bound;
    }
  }

  /// Uniform integer in [lo, hi] (inclusive), lo <= hi.
  constexpr std::uint64_t uniform_range(std::uint64_t lo, std::uint64_t hi) noexcept {
    const std::uint64_t span = hi - lo;
    if (span == ~std::uint64_t{0}) return next();
    return lo + uniform_below(span + 1);
  }

  /// Uniform double in [0, 1) with 53 random bits (exact conversion, no rounding).
  constexpr double uniform01() noexcept { return static_cast<double>(next() >> 11) * 0x1.0p-53; }

  /// Bernoulli(p).
  constexpr bool bernoulli(double p) noexcept { return uniform01() < p; }

 private:
  static constexpr std::uint64_t rotl(std::uint64_t x, int k) noexcept {
    return (x << k) | (x >> (64 - k));
  }
  State s_;
};

/// Derives the seed of an independent sub-stream (e.g. topology vs. weights) from a user seed.
constexpr std::uint64_t substream_seed(std::uint64_t seed, std::uint64_t stream) noexcept {
  SplitMix64 sm(seed ^ (0xD1B54A32D192ED03ULL * (stream + 1)));
  return sm.next();
}

// ------------------------------------------------------------------ deterministic math

/// Natural logarithm built only from correctly rounded IEEE operations (x > 0, finite).
/// Accurate to a few ulp; identical on every IEEE-754 platform.
double det_log(double x) noexcept;
/// Exponential built only from correctly rounded IEEE operations (|x| <= ~709).
double det_exp(double x) noexcept;

}  // namespace sssp::gen
