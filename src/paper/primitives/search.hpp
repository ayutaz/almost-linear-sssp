#pragma once
// Sec 3.1 "Searches and Selection" (Lemma 3.1), Mode A (paper-faithful).
//
// Depth-first search, breadth-first search, connected components and topological sorting cost
// O(k + l) on k vertices and l arcs; a stable counting-sort pass costs O(p + k); BFPRT
// selection costs O(p) comparisons. Searches and sorts use only integer (index) operations;
// selection touches keys only through the comparator, so the comparison-addition model of Sec 2
// is respected.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// One directed arc tail -> head.
struct DiArc {
  Vertex tail;
  Vertex head;
};

/// Directed adjacency lists. Arcs of a tail keep their input order; each arc remembers its input
/// index. Construction is one stable counting-sort pass on the tail (Lemma 3.1), O(k + l).
class Digraph {
 public:
  Digraph() = default;

  /// Throws std::out_of_range for an endpoint >= k, std::length_error for oversized input.
  static Digraph from_arcs(std::size_t k, std::span<const DiArc> arcs);

  [[nodiscard]] std::size_t num_vertices() const noexcept {
    return offsets_.empty() ? 0 : offsets_.size() - 1;
  }
  [[nodiscard]] std::size_t num_arcs() const noexcept { return heads_.size(); }
  [[nodiscard]] std::span<const Vertex> heads(Vertex v) const noexcept {
    return {heads_.data() + offsets_[v], heads_.data() + offsets_[v + 1]};
  }
  [[nodiscard]] std::span<const EdgeId> arc_ids(Vertex v) const noexcept {
    return {ids_.data() + offsets_[v], ids_.data() + offsets_[v + 1]};
  }
  /// Number of arcs entering v (counted once at construction, used as the initial "current
  /// indegree" of the topological sort).
  [[nodiscard]] std::size_t indegree(Vertex v) const noexcept { return indegree_[v]; }

 private:
  std::vector<std::size_t> offsets_;
  std::vector<Vertex> heads_;
  std::vector<EdgeId> ids_;
  std::vector<std::size_t> indegree_;
};

/// Result of a single-root search. `order` lists the vertices reachable from the root in
/// discovery order. `parent`/`parent_edge` (size k) describe the first-discovery tree;
/// the root and unreached vertices have kNoVertex / kNoEdge. For a CsrGraph `parent_edge` is the
/// undirected edge id, for a Digraph it is the input arc index.
struct Traversal {
  std::vector<Vertex> order;
  std::vector<Vertex> parent;
  std::vector<EdgeId> parent_edge;
};

/// Sec 3.1: "marks a vertex at discovery and scans its adjacency list once". DFS yields the
/// preorder of the recursive search that visits arcs in adjacency order; it uses an explicit
/// stack. Counters (optional): vertex_scans per discovered vertex, edge_scans per arc examined.
/// Throws std::out_of_range if root >= k.
Traversal dfs_order(const Digraph& g, Vertex root, OpCounters* ops = nullptr);
Traversal bfs_order(const Digraph& g, Vertex root, OpCounters* ops = nullptr);
template <class W>
Traversal dfs_order(const CsrGraph<W>& g, Vertex root, OpCounters* ops = nullptr);
template <class W>
Traversal bfs_order(const CsrGraph<W>& g, Vertex root, OpCounters* ops = nullptr);

/// Connected components of an undirected graph by repeated BFS. Labels are 0..count-1, numbered
/// in order of each component's smallest vertex.
struct Components {
  std::vector<Vertex> label;
  std::size_t count = 0;
};
template <class W>
Components connected_components(const CsrGraph<W>& g, OpCounters* ops = nullptr);

/// Sec 3.1 / Lemma 3.1: repeatedly remove a vertex of zero current indegree (FIFO, seeded in
/// index order), decrementing the indegrees of its successors. Every arc of the result points
/// forward. Throws std::invalid_argument if the digraph has a cycle (including a self-loop).
std::vector<Vertex> topological_sort(const Digraph& g, OpCounters* ops = nullptr);

// ------------------------------------------------------------------------- stable integer sorts

/// Sec 3.1 / Lemma 3.1: stable counting sort of p records on an integer coordinate in [0, k):
/// count each coordinate, take prefix sums as starting positions, and write each record to its
/// coordinate's next position in input order. O(p + k). `key(rec)` is evaluated once per record
/// and must be an integer index (never a real length: Sec 2 forbids using lengths as addresses).
/// Throws std::out_of_range if some key is >= k.
template <std::ranges::random_access_range R, class KeyFn>
std::vector<std::ranges::range_value_t<R>> counting_sort_stable(const R& recs, std::size_t k,
                                                                KeyFn key) {
  const std::size_t p = std::ranges::size(recs);
  std::vector<std::size_t> start(k + 1, 0);
  std::vector<std::size_t> keys(p);
  for (std::size_t i = 0; i < p; ++i) {
    keys[i] = static_cast<std::size_t>(key(recs[i]));
    if (keys[i] >= k) {
      throw std::out_of_range("counting_sort_stable: key " + std::to_string(keys[i]) +
                              " outside [0, " + std::to_string(k) + ")");
    }
    ++start[keys[i] + 1];
  }
  for (std::size_t c = 0; c < k; ++c) start[c + 1] += start[c];
  std::vector<std::size_t> source(p);  // source[output position] = input position
  for (std::size_t i = 0; i < p; ++i) source[start[keys[i]]++] = i;
  std::vector<std::ranges::range_value_t<R>> out;
  out.reserve(p);
  for (std::size_t j = 0; j < p; ++j) out.push_back(recs[source[j]]);
  return out;
}

/// Sec 3.1 / Lemma 3.1: lexicographic stable sort on c = ranges.size() integer coordinates,
/// coordinate j being `key(rec, j)` in [0, ranges[j]). Applies the stable counting pass from the
/// last coordinate to the first; O(c (p + max k)). Stability of each pass preserves the order
/// of the less significant coordinates. Records equal on every coordinate keep input order.
template <std::ranges::random_access_range R, class KeyFn>
std::vector<std::ranges::range_value_t<R>> radix_sort_stable(const R& recs,
                                                             std::span<const std::size_t> ranges,
                                                             KeyFn key) {
  using T = std::ranges::range_value_t<R>;
  std::vector<T> cur(std::ranges::begin(recs), std::ranges::end(recs));
  for (std::size_t j = ranges.size(); j-- > 0;) {
    cur = counting_sort_stable(cur, ranges[j], [&](const T& r) { return key(r, j); });
  }
  return cur;
}

// ----------------------------------------------------------------------------------- selection

/// Inputs of at most this many keys are solved by direct comparison (Sec 3.1: "Use direct
/// comparison for bounded p"). The paper only requires some fixed bound; the value is our choice.
inline constexpr std::size_t kSelectDirectCutoff = 10;

namespace detail {

template <class T, class Cmp>
void insertion_sort_counted(std::span<T> a, Cmp& cmp, OpCounters* ops) {
  for (std::size_t i = 1; i < a.size(); ++i) {
    for (std::size_t j = i; j > 0; --j) {
      if (ops) SSSP_COUNT(++ops->comparisons);
      if (!cmp(a[j], a[j - 1])) break;
      std::swap(a[j], a[j - 1]);
    }
  }
}

template <class T, class Cmp>
T select_bfprt_rec(std::vector<T> a, std::size_t r, Cmp& cmp, OpCounters* ops) {
  // The second recursive call of Sec 3.1 is a tail call; it is performed by this loop.
  for (;;) {
    if (a.size() <= kSelectDirectCutoff) {
      insertion_sort_counted(std::span<T>(a), cmp, ops);
      return a[r];
    }
    // Sec 3.1: groups of five and at most one smaller group; take each group's median.
    std::vector<T> medians;
    medians.reserve((a.size() + 4) / 5);
    for (std::size_t g = 0; g < a.size(); g += 5) {
      const std::size_t len = std::min<std::size_t>(5, a.size() - g);
      insertion_sort_counted(std::span<T>(a).subspan(g, len), cmp, ops);
      medians.push_back(a[g + (len - 1) / 2]);
    }
    // Sec 3.1: the median of the group medians, recursively.
    const std::size_t mid = (medians.size() - 1) / 2;
    const T pivot = select_bfprt_rec(std::move(medians), mid, cmp, ops);
    // Sec 3.1: partition into the keys below, equal to, and above the pivot.
    std::vector<T> below;
    std::vector<T> above;
    std::size_t equal = 0;
    for (auto& x : a) {
      if (ops) SSSP_COUNT(++ops->comparisons);
      if (cmp(x, pivot)) {
        below.push_back(std::move(x));
        continue;
      }
      if (ops) SSSP_COUNT(++ops->comparisons);
      if (cmp(pivot, x)) {
        above.push_back(std::move(x));
      } else {
        ++equal;
      }
    }
    // Return the pivot if its equal block contains rank r; otherwise continue in the strict side
    // with the adjusted rank. Lemma 3.1 proof: the strict side has at most 7p/10 + O(1) keys,
    // including when keys tie.
    if (r < below.size()) {
      a = std::move(below);
    } else if (r < below.size() + equal) {
      return pivot;
    } else {
      r -= below.size() + equal;
      a = std::move(above);
    }
  }
}

}  // namespace detail

/// Sec 3.1 / Lemma 3.1: the selection method of Blum, Floyd, Pratt, Rivest and Tarjan. Returns
/// (a copy of) a key of rank r (0-based) among the p keys under the strict weak order `cmp`, i.e.
/// a key equivalent to what std::nth_element would place at position r. Worst-case O(p)
/// comparisons and elementary operations; only `cmp` is applied to keys, so it is valid for real
/// lengths in the comparison-addition model (Sec 2). The input is not modified.
/// Counter (optional): ops->comparisons. Throws std::out_of_range if r >= p.
template <std::ranges::random_access_range R, class Cmp = std::ranges::less>
std::ranges::range_value_t<R> select_bfprt(const R& keys, std::size_t r, Cmp cmp = {},
                                           OpCounters* ops = nullptr) {
  const std::size_t p = std::ranges::size(keys);
  if (r >= p) {
    throw std::out_of_range("select_bfprt: rank " + std::to_string(r) + " with " +
                            std::to_string(p) + " keys");
  }
  using T = std::ranges::range_value_t<R>;
  return detail::select_bfprt_rec(std::vector<T>(std::ranges::begin(keys), std::ranges::end(keys)),
                                  r, cmp, ops);
}

}  // namespace sssp::paper
