// Differential test of Algorithm 4 (LocalHierarchy, Sec 5.4, Lemma 5.10) with the paper's
// Mode A constants on an input whose auxiliary tree distances leave 64 bits although every
// distance from the source fits (regression 20261127-owner-decomposition-overflow at the
// paper's size).
//
// Graph: the path 0 - 1 - ... - (k - 1) with k = 209953 (the first size of the radial branch,
// L = 18, b = 18^4); the first four edges have length 2^62, the others length 1. Source 2.
// The MST is the path itself, rooted at vertex 0 (Q-31): D_T(0, 4) = 2^64, and every later root
// distance exceeds 64 bits.
//
// Exact reference: on a path the distance is the sum of the edges between the endpoints,
//   d(2, v) = |v - 2| * 2^62 for v <= 4,   d(2, v) = 2^63 + (v - 4) for v >= 4,
// at most 2^63 + 209948. The same values are printed by tests/python/reference/exact_sssp.py
// on this edge list (checked when the test was written; the closed form is used here because
// the C++ test cannot run the Python reference).
//
// Checked: the radial branch runs with no FallbackRecord; the hierarchy satisfies Def 4.1
// (project validator); Lemma 5.10 path realization: the least state-graph distance from the
// native entry of the source over the copies of v equals d(2, v) at every vertex, and no copy
// is below it. The state-graph search is a Dijkstra of the test support (shares no code with
// src/paper); it is an oracle on the output, not a step of the algorithm.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "property/paper/local_hierarchy_large_support.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace {

using W = std::uint64_t;
using sssp::Vertex;
namespace ll = sssp::test::local_large;

constexpr W kQ = W{1} << 62;
constexpr Vertex kSource = 2;

std::vector<sssp::Edge<W>> heavy_prefix_path(std::size_t k) {
  std::vector<sssp::Edge<W>> edges;
  edges.reserve(k - 1);
  for (Vertex v = 0; v + 1 < k; ++v) edges.push_back({v, v + 1, v < 4 ? kQ : W{1}});
  return edges;
}

/// d(2, v) on heavy_prefix_path (closed form, see the header).
W exact_distance(Vertex v) {
  if (v <= 4) return (v < kSource ? kSource - v : v - kSource) * kQ;
  return 2 * kQ + (v - 4);
}

}  // namespace

TEST_CASE("Alg 4 Mode A at k = 209953: four edges of 2^62 at the root of the MST",
          "[slow][local-large][overflow]") {
  const std::size_t k = ll::kFirstRadial;
  const auto edges = heavy_prefix_path(k);
  const auto g = sssp::CsrGraph<W>::from_edges(k, edges);
  REQUIRE(sssp::paper::local_construction_applies(k, sssp::paper::kPaperLocalParams));

  sssp::Diagnostics diag;
  sssp::paper::LocalHierarchyResult<W> lh;
  REQUIRE_NOTHROW(lh = sssp::paper::local_hierarchy(g, &diag));
  CHECK(diag.fallbacks.empty());
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters.at("local_hierarchy_radial") == 1);
  CHECK(diag.counters.count("local_hierarchy_base_case") == 0);
#endif
  REQUIRE(lh.native_entry.size() == k);

  const auto& h = lh.hierarchy;
  const auto valid = sssp::paper::validate_ordered_hierarchy(h, g);
  INFO(valid.summary());
  CHECK(valid.ok());

  // Lemma 5.10: every tentative length of this search is a path length from the source, at
  // most 2^63 + 209948 + 2^62 < 2^64, so the unchecked sums of the support Dijkstra are exact.
  const auto ds = ll::dijkstra(ll::state_digraph(h), {{lh.native_entry[kSource], W{0}}});
  std::vector<W> best(k, ll::infinity<W>());
  std::size_t below = 0;
  for (std::size_t x = 0; x < ds.size(); ++x) {
    const Vertex v = h.projection[x];
    if (ds[x] < exact_distance(v)) ++below;
    if (ds[x] < best[v]) best[v] = ds[x];
  }
  CHECK(below == 0);
  std::size_t mismatches = 0;
  Vertex first_bad = sssp::kNoVertex;
  for (Vertex v = 0; v < k; ++v) {
    if (best[v] != exact_distance(v)) {
      if (mismatches == 0) first_bad = v;
      ++mismatches;
    }
  }
  INFO("first mismatching vertex " << first_bad);
  CHECK(mismatches == 0);
}
