#pragma once
// Randomness in the model of Sec 2 (arXiv 2611.34911v1, txt l.307-310): independent uniformly
// random words; a uniform integer in 0..h-1 takes ceil(log2 h) unbiased bits and rejects
// outcomes >= h (each attempt succeeds with probability > 1/2); h = 1 returns zero.
// Only Sec 6 consumes randomness (sampling, uniform permutations); callers pass the generator
// explicitly and fix the seed, so runs are reproducible on every platform.
//
// std::uniform_int_distribution / std::shuffle are not used: their output is
// implementation-defined, and reduction modulo h is biased.

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/model/ilog.hpp"

namespace sssp::paper::model {

/// SplitMix64 (Steele, Lea, Flood 2014): 64-bit words from a 64-bit state; fully specified by
/// integer arithmetic, hence identical on every platform. There is no default seed.
class WordRng {
 public:
  explicit constexpr WordRng(std::uint64_t seed) noexcept : state_(seed) {}

  /// One uniformly random 64-bit word.
  constexpr std::uint64_t next_word() noexcept {
    std::uint64_t z = (state_ += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
  }

  /// c unbiased bits (0 <= c <= 64): the c most significant bits of the next word. c = 0 draws
  /// nothing and returns 0. Throws std::invalid_argument for c > 64.
  constexpr std::uint64_t next_bits(unsigned c) {
    if (c > 64) throw std::invalid_argument("next_bits: at most 64 bits per word");
    if (c == 0) return 0;
    const std::uint64_t w = next_word();
    return c == 64 ? w : (w >> (64 - c));
  }

 private:
  std::uint64_t state_;
};

/// Anything that yields c unbiased bits on request.
template <class S>
concept BitSource = requires(S& s, unsigned c) {
  { s.next_bits(c) } -> std::convertible_to<std::uint64_t>;
};

/// Uniform integer in 0..h-1 (Sec 2, l.307-308): h = 1 returns 0 without drawing; otherwise draw
/// ceil(log2 h) bits and reject outcomes >= h. Throws std::invalid_argument for h = 0.
template <BitSource S>
std::uint64_t uniform_below(S& src, std::uint64_t h) {
  if (h == 0) throw std::invalid_argument("uniform_below: empty range");
  if (h == 1) return 0;
  const unsigned c = ceil_log2(h);
  for (;;) {
    const std::uint64_t x = static_cast<std::uint64_t>(src.next_bits(c));
    if (x < h) return x;
  }
}

/// Uniformly random order of a (Fisher-Yates: position i receives a uniform element of 0..i).
template <class T, BitSource S>
void random_shuffle(std::span<T> a, S& src) {
  for (std::size_t i = a.size(); i > 1; --i) {
    const auto j = static_cast<std::size_t>(uniform_below(src, i));
    using std::swap;
    swap(a[i - 1], a[j]);
  }
}

/// Uniformly random permutation of 0..n-1 (vertex / edge index range). Throws std::length_error
/// if n exceeds the 32-bit index range.
template <BitSource S>
std::vector<std::uint32_t> random_permutation(std::size_t n, S& src) {
  if (n > (std::size_t{1} << 32)) throw std::length_error("random_permutation: n > 2^32");
  std::vector<std::uint32_t> p(n);
  for (std::size_t i = 0; i < n; ++i) p[i] = static_cast<std::uint32_t>(i);
  random_shuffle(std::span<std::uint32_t>(p), src);
  return p;
}

}  // namespace sssp::paper::model
