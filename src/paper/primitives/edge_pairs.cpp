// Lightest edge per vertex pair (Alg 4 l.4, Alg 5 l.4, Alg 8 l.5). See edge_pairs.hpp.
#include "paper/primitives/edge_pairs.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace sssp::paper {

template <class W>
LightestEdgePairs lightest_edge_pairs(std::size_t n, std::span<const Edge<W>> edges) {
  validate_edges<W>(n, edges);
  LightestEdgePairs out;
  out.representative.assign(edges.size(), kNoEdge);

  // Non-loop edges in increasing index order (loops are dropped).
  std::vector<EdgeId> ids;
  ids.reserve(edges.size());
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (edges[i].u != edges[i].v) ids.push_back(static_cast<EdgeId>(i));
  }
  auto lo = [&](EdgeId id) -> std::size_t { return std::min(edges[id].u, edges[id].v); };
  auto hi = [&](EdgeId id) -> std::size_t { return std::max(edges[id].u, edges[id].v); };

  // Lemma 3.1: two stable counting-sort passes, first by the larger endpoint, then by the
  // smaller, group equal unordered pairs. Stability keeps the indices increasing inside a group.
  std::vector<EdgeId> buf(ids.size());
  std::vector<std::size_t> count(n + 1, 0);
  auto stable_pass = [&](auto key) {
    std::fill(count.begin(), count.end(), 0);
    for (const EdgeId id : ids) ++count[key(id) + 1];
    for (std::size_t v = 0; v < n; ++v) count[v + 1] += count[v];
    for (const EdgeId id : ids) buf[count[key(id)]++] = id;
    ids.swap(buf);
  };
  stable_pass(hi);
  stable_pass(lo);

  // One scan per group: indices increase inside a group, so the first edge of strictly least
  // length is the (length, index)-least one.
  std::vector<std::uint8_t> keep(edges.size(), 0);
  for (std::size_t i = 0; i < ids.size();) {
    EdgeId best = ids[i];
    std::size_t j = i + 1;
    for (; j < ids.size() && lo(ids[j]) == lo(best) && hi(ids[j]) == hi(best); ++j) {
      if (edges[ids[j]].weight < edges[best].weight) best = ids[j];
    }
    for (std::size_t k = i; k < j; ++k) out.representative[ids[k]] = best;
    keep[best] = 1;
    i = j;
  }
  for (std::size_t i = 0; i < keep.size(); ++i) {
    if (keep[i] != 0) out.kept.push_back(static_cast<EdgeId>(i));
  }
  return out;
}

template <class W>
std::vector<EdgeId> lightest_edge_per_pair(std::size_t n, std::span<const Edge<W>> edges) {
  return lightest_edge_pairs<W>(n, edges).kept;
}

template LightestEdgePairs lightest_edge_pairs<std::uint64_t>(std::size_t,
                                                              std::span<const Edge<std::uint64_t>>);
template LightestEdgePairs lightest_edge_pairs<double>(std::size_t, std::span<const Edge<double>>);
template std::vector<EdgeId> lightest_edge_per_pair<std::uint64_t>(
    std::size_t, std::span<const Edge<std::uint64_t>>);
template std::vector<EdgeId> lightest_edge_per_pair<double>(std::size_t,
                                                            std::span<const Edge<double>>);

}  // namespace sssp::paper
