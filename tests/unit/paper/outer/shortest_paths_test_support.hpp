#pragma once
// Test support for the outer construction (Alg 8 l.8-17 of arXiv 2611.34911v1, Sec 8.3;
// P3-OUTER-PIPELINE): connected test inputs built here (independent of the library's
// generators), the hand examples L0-L4 of docs/research/04_outer_scales_pipeline.md Sec 3.5,
// and the comparison with the independent oracle (Bellman-Ford on the same graph).

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "paper/outer/shortest_paths.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/graphs.hpp"
#include "support/prng.hpp"

namespace sssp::test::pipeline {

/// Mode B style parameters (D-7): a constant separator b, everything else the paper's.
inline paper::outer::OuterConstructionParams constant_separator(std::uint64_t b) {
  paper::outer::OuterConstructionParams p;
  p.separator = paper::outer::OuterSeparator{b, 1};
  return p;
}

/// A connected loopless multigraph: a random recursive tree plus `extra` random nonloop edges
/// (parallel edges possible), weights uniform in [lo, hi] (lo >= 1).
template <class W>
std::vector<Edge<W>> connected_edges(std::size_t n, std::size_t extra, std::uint64_t lo,
                                     std::uint64_t hi, std::uint64_t seed) {
  std::vector<Edge<W>> e = random_tree_edges<W>(n, lo, hi, seed);
  if (n < 2) return e;
  SplitMix64 rng(seed ^ 0x5bd1e995u);
  for (std::size_t i = 0; i < extra; ++i) {
    const auto u = static_cast<Vertex>(rng.below(n));
    auto v = static_cast<Vertex>(rng.below(n - 1));
    if (v >= u) ++v;
    e.push_back({u, v, weight_of<W>(rng.between(lo, hi))});
  }
  return e;
}

/// A named connected input.
template <class W>
struct Input {
  std::string name;
  std::size_t n = 0;
  std::vector<Edge<W>> edges;
};

/// Level L0-L4 of note 04 Sec 3.5 with its separator size.
template <class W>
struct Level {
  std::string name;
  std::size_t n = 0;
  std::vector<Edge<W>> edges;
  std::uint64_t b = 1;
  std::vector<W> dist_from_0;  ///< by hand
};

template <class W>
std::vector<Level<W>> levels() {
  std::vector<Level<W>> l;
  // L0: one vertex.
  l.push_back({"L0", 1, {}, 1, {W{0}}});
  // L1: one edge, both ends pivots (b = 1).
  l.push_back({"L1", 2, {{0, 1, W{7}}}, 1, {W{0}, W{7}}});
  // L2: the path 0 - 1 - 2 - 3 (10, 3, 3), b = 2: pivots {0, 2}, nonpivots 1 and 3.
  l.push_back(
      {"L2", 4, {{0, 1, W{10}}, {1, 2, W{3}}, {2, 3, W{3}}}, 2, {W{0}, W{10}, W{13}, W{16}}});
  // L3: L2 plus the chord 1 - 3 (5): d(0, 3) = 10 + 5.
  l.push_back({"L3",
               4,
               {{0, 1, W{10}}, {1, 2, W{3}}, {2, 3, W{3}}, {1, 3, W{5}}},
               2,
               {W{0}, W{10}, W{13}, W{15}}});
  // L4: the tree 0-1 (4), 1-2 (4), 0-3 (100), 3-4 (8), 4-5 (1), b = 3.
  l.push_back({"L4",
               6,
               {{0, 1, W{4}}, {1, 2, W{4}}, {0, 3, W{100}}, {3, 4, W{8}}, {4, 5, W{1}}},
               3,
               {W{0}, W{4}, W{8}, W{100}, W{108}, W{109}}});
  return l;
}

/// dist equals Bellman-Ford on g from s, exactly.
template <class W>
void require_oracle(const CsrGraph<W>& g, Vertex s, const std::vector<W>& dist) {
  const auto oracle = baseline::bellman_ford(g, s);
  REQUIRE(dist.size() == oracle.dist.size());
  for (Vertex v = 0; v < dist.size(); ++v) {
    INFO("source " << s << " vertex " << v);
    REQUIRE(dist[v] == oracle.dist[v]);
  }
}

/// Number of edges of class `c`.
template <class W>
std::size_t count_class(const paper::owner::OuterEdgeClassification<W>& cls,
                        paper::owner::OuterEdgeClass c) {
  std::size_t k = 0;
  for (const auto x : cls.edge_class) k += x == c ? 1 : 0;
  return k;
}

/// Number of blocks of kind `k` with at least two states.
template <class W>
std::size_t count_blocks(const paper::outer::OuterScales<W>& sc, paper::outer::OuterBlockKind k) {
  std::size_t c = 0;
  for (const auto& b : sc.blocks) c += (b.kind == k && b.size >= 2) ? 1 : 0;
  return c;
}

}  // namespace sssp::test::pipeline
