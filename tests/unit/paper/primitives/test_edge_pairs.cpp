// Tests for the lightest edge per vertex pair (Alg 4 l.4, Alg 5 l.4, Alg 8 l.5; Sec 2
// l.286-288): one (weight, id)-least survivor per unordered pair, loops dropped, and a pointer
// from every original edge to its pair's survivor. The reference is a std::map keyed by the
// unordered pair, written here independently of the library.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <map>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/primitives/edge_pairs.hpp"
#include "sssp/graph.hpp"

using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::Vertex;
using sssp::paper::lightest_edge_pairs;
using sssp::paper::lightest_edge_per_pair;

namespace {

template <class W>
std::vector<Edge<W>> random_multigraph(std::size_t n, std::size_t m, std::uint64_t max_weight,
                                       std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::vector<Edge<W>> edges;
  for (std::size_t i = 0; i < m; ++i) {
    const auto u = static_cast<Vertex>(rng() % n);
    const auto v = static_cast<Vertex>(rng() % n);
    edges.push_back({u, v, static_cast<W>(rng() % (max_weight + 1))});
  }
  return edges;
}

// Reference: survivor per unordered pair by a std::map, ties by the least index.
template <class W>
std::pair<std::vector<EdgeId>, std::vector<EdgeId>> map_reference(
    const std::vector<Edge<W>>& edges) {
  std::map<std::pair<Vertex, Vertex>, EdgeId> best;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u == e.v) continue;
    const auto key = std::make_pair(std::min(e.u, e.v), std::max(e.u, e.v));
    const auto it = best.find(key);
    if (it == best.end() || e.weight < edges[it->second].weight) {
      best[key] = static_cast<EdgeId>(i);  // equal weight: the earlier (smaller) index stays
    }
  }
  std::vector<EdgeId> kept;
  for (const auto& [key, id] : best) kept.push_back(id);
  std::sort(kept.begin(), kept.end());
  std::vector<EdgeId> rep(edges.size(), kNoEdge);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u != e.v) rep[i] = best.at({std::min(e.u, e.v), std::max(e.u, e.v)});
  }
  return {kept, rep};
}

}  // namespace

TEMPLATE_TEST_CASE("lightest_edge_pairs on a hand example", "[paper][edge_pairs]", std::uint64_t,
                   double) {
  using T = TestType;
  // Pair {0,1}: ids 0 (5), 1 (3), 2 (3) -> id 1 (equal weight, smaller id). Loop 3 dropped.
  // Pair {1,2}: ids 4, 5 (both 4) -> id 4. Pair {0,2}: id 6 alone.
  const std::vector<Edge<T>> edges{{0, 1, T(5)}, {1, 0, T(3)}, {0, 1, T(3)}, {2, 2, T(1)},
                                   {1, 2, T(4)}, {2, 1, T(4)}, {0, 2, T(9)}};
  const auto r = lightest_edge_pairs<T>(3, edges);
  REQUIRE(r.kept == std::vector<EdgeId>{1, 4, 6});
  REQUIRE(r.representative == std::vector<EdgeId>{1, 1, 1, kNoEdge, 4, 4, 6});
  REQUIRE(lightest_edge_per_pair<T>(3, edges) == r.kept);
}

TEMPLATE_TEST_CASE("lightest_edge_pairs boundary inputs", "[paper][edge_pairs]", std::uint64_t,
                   double) {
  using T = TestType;
  const auto empty = lightest_edge_pairs<T>(0, std::vector<Edge<T>>{});
  REQUIRE(empty.kept.empty());
  REQUIRE(empty.representative.empty());
  // Only loops: nothing kept, every edge points nowhere.
  const auto loops = lightest_edge_pairs<T>(2, std::vector<Edge<T>>{{0, 0, T(1)}, {1, 1, T(0)}});
  REQUIRE(loops.kept.empty());
  REQUIRE(loops.representative == std::vector<EdgeId>{kNoEdge, kNoEdge});
  // Isolated vertices and a single edge.
  const auto one = lightest_edge_pairs<T>(5, std::vector<Edge<T>>{{4, 2, T(0)}});
  REQUIRE(one.kept == std::vector<EdgeId>{0});
  REQUIRE(one.representative == std::vector<EdgeId>{0});
  // Endpoint out of range.
  REQUIRE_THROWS_AS(lightest_edge_pairs<T>(1, std::vector<Edge<T>>{{0, 1, T(1)}}),
                    std::out_of_range);
}

TEMPLATE_TEST_CASE("lightest_edge_pairs equals a std::map reference on random multigraphs",
                   "[paper][edge_pairs]", std::uint64_t, double) {
  using T = TestType;
  for (std::uint64_t seed = 1; seed <= 300; ++seed) {
    const std::size_t n = 1 + seed % 23;
    const std::size_t m = (seed * 7) % (5 * n + 1);
    const std::uint64_t max_weight = seed % 5;  // many ties, including all-zero weights
    const auto edges = random_multigraph<T>(n, m, max_weight, seed);
    const auto r = lightest_edge_pairs<T>(n, edges);
    const auto [kept, rep] = map_reference(edges);
    INFO("seed " << seed << " n " << n << " m " << m);
    REQUIRE(r.kept == kept);
    REQUIRE(r.representative == rep);
    // Every representative is itself kept and joins the same unordered pair.
    for (std::size_t i = 0; i < edges.size(); ++i) {
      if (rep[i] == kNoEdge) continue;
      REQUIRE(r.representative[r.representative[i]] == r.representative[i]);
    }
  }
}

TEST_CASE("lightest_edge_pairs on many vertices and few edges", "[paper][edge_pairs]") {
  // O(m + N): the counting passes are indexed by vertex, so a sparse input on a large vertex
  // range must still be handled (and quickly).
  const std::size_t n = 1'000'000;
  const std::vector<Edge<std::uint64_t>> edges{
      {999'999, 0, 7}, {0, 999'999, 7}, {5, 5, 1}, {123'456, 654'321, 2}, {654'321, 123'456, 1}};
  const auto r = lightest_edge_pairs<std::uint64_t>(n, edges);
  REQUIRE(r.kept == std::vector<EdgeId>{0, 4});
  REQUIRE(r.representative == std::vector<EdgeId>{0, 0, kNoEdge, 4, 4});
}
