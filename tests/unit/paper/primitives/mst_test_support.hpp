#pragma once
// Test helpers for the Sec 3.2 MST primitives: seeded random multigraphs and naive checks that
// do not use the code under test.

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>
#include <type_traits>
#include <vector>

#include "sssp/graph.hpp"

namespace mst_test {

/// Weight r for integer types; r * 0.1 for double (not dyadic: only comparisons are made).
template <class W>
W make_weight(std::uint64_t r) {
  if constexpr (std::is_floating_point_v<W>) {
    return static_cast<W>(r) * W{0.1};
  } else {
    return static_cast<W>(r);
  }
}

/// Random multigraph on n vertices: if `connected`, a random spanning tree first, then
/// `extra` random edges (parallel edges and loops allowed). Weights in [0, max_weight];
/// a small max_weight produces many ties. The edge list is shuffled so ids are not ordered.
template <class W>
std::vector<sssp::Edge<W>> random_edges(std::size_t n, std::size_t extra, std::uint64_t seed,
                                        std::uint64_t max_weight, bool connected) {
  std::mt19937_64 rng(seed);
  std::vector<sssp::Edge<W>> edges;
  auto weight = [&] { return make_weight<W>(rng() % (max_weight + 1)); };
  if (connected) {
    for (std::size_t v = 1; v < n; ++v) {
      const auto u = static_cast<sssp::Vertex>(rng() % v);
      edges.push_back({u, static_cast<sssp::Vertex>(v), weight()});
    }
  }
  if (n > 0) {
    for (std::size_t i = 0; i < extra; ++i) {
      const auto u = static_cast<sssp::Vertex>(rng() % n);
      const auto v = static_cast<sssp::Vertex>(rng() % n);
      edges.push_back({u, v, weight()});
    }
  }
  std::shuffle(edges.begin(), edges.end(), rng);
  return edges;
}

/// Naive union-find with full relabelling (O(n) per union), independent of the library.
class NaiveLabels {
 public:
  explicit NaiveLabels(std::size_t n) : label_(n) { std::iota(label_.begin(), label_.end(), 0); }
  [[nodiscard]] std::size_t label(std::size_t v) const { return label_[v]; }
  bool unite(std::size_t a, std::size_t b) {
    const std::size_t la = label_[a];
    const std::size_t lb = label_[b];
    if (la == lb) return false;
    for (auto& l : label_) {
      if (l == lb) l = la;
    }
    return true;
  }

 private:
  std::vector<std::size_t> label_;
};

template <class W>
std::size_t count_components(std::size_t n, const std::vector<sssp::Edge<W>>& edges) {
  NaiveLabels d(n);
  std::size_t c = n;
  for (const auto& e : edges) {
    if (d.unite(e.u, e.v)) --c;
  }
  return c;
}

}  // namespace mst_test
