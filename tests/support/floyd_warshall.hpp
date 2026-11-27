#pragma once
// All-pairs Floyd–Warshall reference for tiny graphs. It shares no code with the library's
// SSSP implementations and works from the edge list, not the CSR, so it is an independent
// oracle. Unreachable pairs are WeightTraits<W>::infinity(); integer overflow throws.
//
// Note: for floating-point weights FW sums sub-paths in a different association than
// relaxation-based SSSP, so results may differ by rounding (see docs/design/numerics.md).

#include <cstddef>
#include <span>
#include <vector>

#include "sssp/graph.hpp"

namespace sssp::test {

template <class W>
std::vector<std::vector<W>> floyd_warshall(std::size_t n, std::span<const Edge<W>> edges) {
  using T = WeightTraits<W>;
  std::vector<std::vector<W>> d(n, std::vector<W>(n, T::infinity()));
  for (std::size_t i = 0; i < n; ++i) d[i][i] = T::zero();
  for (const auto& e : edges) {
    if (e.weight < d[e.u][e.v]) {
      d[e.u][e.v] = e.weight;
      d[e.v][e.u] = e.weight;
    }
  }
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      if (d[i][k] == T::infinity()) continue;
      for (std::size_t j = 0; j < n; ++j) {
        if (d[k][j] == T::infinity()) continue;
        const W via = T::checked_add(d[i][k], d[k][j]);
        if (via < d[i][j]) d[i][j] = via;
      }
    }
  }
  return d;
}

}  // namespace sssp::test
