#pragma once
// Exhaustive enumeration of small weighted simple graphs: every subset of the C(n,2) vertex
// pairs, times every assignment of a weight from `weights` to each chosen pair.
// The number of graphs is (|weights| + 1)^(n(n-1)/2).

#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "sssp/graph.hpp"

namespace sssp::test {

[[nodiscard]] constexpr std::uint64_t count_small_graphs(std::size_t n, std::size_t num_weights) {
  std::uint64_t total = 1;
  const std::size_t pairs = n == 0 ? 0 : n * (n - 1) / 2;
  for (std::size_t i = 0; i < pairs; ++i) total *= num_weights + 1;
  return total;
}

/// Calls f(const std::vector<Edge<W>>&) for every weighted simple graph on n vertices.
template <class W, class F>
void for_each_small_graph(std::size_t n, std::span<const W> weights, F&& f) {
  std::vector<std::pair<Vertex, Vertex>> pairs;
  for (Vertex u = 0; u < n; ++u) {
    for (Vertex v = u + 1; v < n; ++v) pairs.emplace_back(u, v);
  }
  const std::size_t k = weights.size();
  std::vector<Edge<W>> edges;
  // Mixed radix over pairs: digit 0 = absent, digit d > 0 = present with weights[d-1].
  std::vector<std::size_t> digit(pairs.size(), 0);
  while (true) {
    edges.clear();
    for (std::size_t i = 0; i < pairs.size(); ++i) {
      if (digit[i] != 0) edges.push_back({pairs[i].first, pairs[i].second, weights[digit[i] - 1]});
    }
    f(std::as_const(edges));
    std::size_t i = 0;
    while (i < digit.size() && digit[i] == k) digit[i++] = 0;
    if (i == digit.size()) break;
    ++digit[i];
  }
}

}  // namespace sssp::test
