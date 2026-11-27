#pragma once
// Test support for the Sec 5.4 assembly: a connected graph family whose light clusters are
// joined by edges of weight about 2^54, so that some blocks are unaligned (sigma_U > N^10 s_U
// for N <= 64, Sec 5.4 l.1206-1208). Written independently of src/paper.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "sssp/graph.hpp"
#include "support/prng.hpp"

namespace sssp::test::local {

/// Clustered graph: clusters of 3..8 vertices with light edges (1..3, plus a few chords),
/// joined in a random tree by heavy edges 2^U[54, 56], plus a few extra heavy edges. Path sums
/// stay below 2^63 for n <= 64.
inline std::vector<Edge<std::uint64_t>> clustered(std::size_t n, std::uint64_t seed) {
  SplitMix64 rng(seed);
  std::vector<Edge<std::uint64_t>> e;
  std::vector<std::size_t> start{0};
  while (start.back() < n)
    start.push_back(std::min<std::size_t>(n, start.back() + 3 + rng.below(6)));
  const std::size_t k = start.size() - 1;
  for (std::size_t c = 0; c < k; ++c) {
    for (std::size_t v = start[c] + 1; v < start[c + 1]; ++v) {
      e.push_back({static_cast<Vertex>(start[c] + rng.below(v - start[c])), static_cast<Vertex>(v),
                   rng.between(1, 3)});
    }
    if (start[c + 1] - start[c] >= 3) {
      e.push_back({static_cast<Vertex>(start[c]), static_cast<Vertex>(start[c + 1] - 1),
                   rng.between(1, 3)});
    }
  }
  const auto heavy = [&] { return std::uint64_t{1} << rng.between(54, 56); };
  const auto pick = [&](std::size_t c) {
    return static_cast<Vertex>(start[c] + rng.below(start[c + 1] - start[c]));
  };
  for (std::size_t c = 1; c < k; ++c) e.push_back({pick(rng.below(c)), pick(c), heavy()});
  for (std::size_t i = 0; i < k / 2; ++i)
    e.push_back({pick(rng.below(k)), pick(rng.below(k)), heavy()});
  return e;
}

}  // namespace sssp::test::local
