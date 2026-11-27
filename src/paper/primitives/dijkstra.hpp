#pragma once
// Sec 3.1 / Lemma 3.2: Dijkstra's search for non-negative arc lengths, Mode A (paper-faithful).
//
// "Start a source at zero and all other tentative values at infinity. Extract a vertex of
// smallest finite unsettled value, settle it, and offer that value plus each outgoing arc length
// to the arc's head. Stop when the queue is empty, or when the smallest finite value exceeds a
// specified radius." (Sec 3.1, paraphrased.) Consequently a vertex at distance exactly R is
// settled; the search stops only on a value strictly greater than R (research note A13).
//
// Two queue disciplines are provided, both from Lemma 3.2:
//   * an indexed binary heap, O(1 + log(k + 1)) per insertion, decrease or extraction;
//   * a scan for each minimum, O(k^2 + l) for a full search (used by the outer base case,
//     Alg 8 l.5).
// "Unvisited entries in a reused array carry a search identifier, so each new search need only
// initialize the entries it touches" (Sec 3.1): SearchArena implements this.
//
// Lengths are only copied, compared and added (comparison-addition model, Sec 2).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// Indexed binary min-heap over vertex ids. Lemma 3.2 proof: each operation restores heap order
/// along one root-leaf path, recording each moved vertex's array address (`pos_`).
/// `pos_` is never cleared: membership of v is certified by heap_[pos_[v]].v == v, so reset()
/// costs O(1) beyond growing the array.
template <class W>
class IndexedBinaryHeap {
 public:
  struct Entry {
    W key;
    Vertex v;
  };

  /// Empties the heap and makes room for vertex ids < n.
  void reset(std::size_t n) {
    heap_.clear();
    if (pos_.size() < n) pos_.resize(n, 0);
  }
  /// Optional comparison counter (ops->comparisons).
  void set_counters(OpCounters* ops) noexcept { ops_ = ops; }

  [[nodiscard]] bool empty() const noexcept { return heap_.empty(); }
  [[nodiscard]] std::size_t size() const noexcept { return heap_.size(); }
  [[nodiscard]] bool contains(Vertex v) const noexcept {
    if (v >= pos_.size()) return false;
    const std::size_t p = pos_[v];
    return p < heap_.size() && heap_[p].v == v;
  }
  [[nodiscard]] const Entry& top() const noexcept { return heap_.front(); }

  /// Precondition: !contains(v), v < n of the last reset().
  void push(Vertex v, W key) {
    heap_.push_back({key, v});
    sift_up(heap_.size() - 1);
  }
  /// Precondition: contains(v) and key is not greater than v's current key.
  void decrease(Vertex v, W key) {
    const std::size_t i = pos_[v];
    heap_[i].key = key;
    sift_up(i);
  }
  /// Precondition: !empty().
  Entry pop_min() {
    const Entry min = heap_.front();
    const Entry last = heap_.back();
    heap_.pop_back();
    if (!heap_.empty()) {
      heap_.front() = last;
      sift_down(0);
    }
    return min;
  }

 private:
  bool less(const W& a, const W& b) const noexcept {
    if (ops_) SSSP_COUNT(++ops_->comparisons);
    return a < b;
  }
  void place(std::size_t i, const Entry& e) noexcept {
    heap_[i] = e;
    pos_[e.v] = i;
  }
  void sift_up(std::size_t i) noexcept {
    const Entry e = heap_[i];
    while (i > 0) {
      const std::size_t p = (i - 1) / 2;
      if (!less(e.key, heap_[p].key)) break;
      place(i, heap_[p]);
      i = p;
    }
    place(i, e);
  }
  void sift_down(std::size_t i) noexcept {
    const Entry e = heap_[i];
    const std::size_t n = heap_.size();
    for (;;) {
      std::size_t c = 2 * i + 1;
      if (c >= n) break;
      if (c + 1 < n && less(heap_[c + 1].key, heap_[c].key)) ++c;
      if (!less(heap_[c].key, e.key)) break;
      place(i, heap_[c]);
      i = c;
    }
    place(i, e);
  }

  std::vector<Entry> heap_;
  std::vector<std::size_t> pos_;
  OpCounters* ops_ = nullptr;
};

/// Per-vertex search state reused across searches (Sec 3.1, search identifiers). An entry
/// belongs to the current search iff its stamp equals the current identifier; begin_search()
/// only advances the identifier, so a search costs O(entries it touches). When the identifier
/// would wrap, all stamps are cleared once (counted by full_resets()). `Stamp` is a template
/// parameter only so that tests can force wrap-around with a narrow type.
template <class W, class Stamp = std::uint32_t>
class SearchArena {
  static_assert(std::is_unsigned_v<Stamp>, "search identifiers are unsigned words");

 public:
  SearchArena() = default;
  explicit SearchArena(std::size_t n) { reserve(n); }

  /// Grows the arrays to hold vertex ids < n (new entries are untouched).
  void reserve(std::size_t n) {
    if (n > stamp_.size()) {
      stamp_.resize(n, 0);
      entry_.resize(n);
    }
  }
  [[nodiscard]] std::size_t capacity() const noexcept { return stamp_.size(); }

  /// Starts a new search over vertex ids < n.
  void begin_search(std::size_t n) {
    reserve(n);
    if (cur_ == std::numeric_limits<Stamp>::max()) {
      std::fill(stamp_.begin(), stamp_.end(), Stamp{0});
      cur_ = 0;
      ++full_resets_;
    }
    ++cur_;
    touched_.clear();
    settled_.clear();
    scratch_.clear();
    heap_.reset(n);
  }

  // ---- queries on the current search
  [[nodiscard]] bool touched(Vertex v) const noexcept {
    return cur_ != 0 && v < stamp_.size() && stamp_[v] == cur_;
  }
  [[nodiscard]] bool settled(Vertex v) const noexcept { return touched(v) && entry_[v].settled; }
  /// Settled distance, tentative value of a touched unsettled vertex, or infinity.
  [[nodiscard]] W dist(Vertex v) const noexcept {
    return touched(v) ? entry_[v].dist : WeightTraits<W>::infinity();
  }
  [[nodiscard]] Vertex parent(Vertex v) const noexcept {
    return touched(v) ? entry_[v].parent : kNoVertex;
  }
  [[nodiscard]] EdgeId parent_edge(Vertex v) const noexcept {
    return touched(v) ? entry_[v].parent_edge : kNoEdge;
  }
  /// Vertices touched by the current search, in order of first touch.
  [[nodiscard]] std::span<const Vertex> touched_vertices() const noexcept { return touched_; }
  /// Vertices settled by the current search, in extraction order (nondecreasing distance).
  [[nodiscard]] std::span<const Vertex> settled_order() const noexcept { return settled_; }
  [[nodiscard]] Stamp search_id() const noexcept { return cur_; }
  [[nodiscard]] std::uint64_t full_resets() const noexcept { return full_resets_; }

  // ---- mutation, for searches built on the arena
  /// First touch of v in this search (precondition: !touched(v)).
  void discover(Vertex v, W d, Vertex p, EdgeId e) {
    stamp_[v] = cur_;
    entry_[v] = Entry{d, p, e, false};
    touched_.push_back(v);
  }
  /// Lower the tentative value of a touched, unsettled vertex.
  void improve(Vertex v, W d, Vertex p, EdgeId e) noexcept {
    entry_[v].dist = d;
    entry_[v].parent = p;
    entry_[v].parent_edge = e;
  }
  void settle(Vertex v) {
    entry_[v].settled = true;
    settled_.push_back(v);
  }
  [[nodiscard]] IndexedBinaryHeap<W>& heap() noexcept { return heap_; }
  /// Reusable vertex list (the scan version keeps its frontier here).
  [[nodiscard]] std::vector<Vertex>& scratch() noexcept { return scratch_; }

 private:
  struct Entry {
    W dist{};
    Vertex parent = kNoVertex;
    EdgeId parent_edge = kNoEdge;
    bool settled = false;
  };
  std::vector<Stamp> stamp_;
  std::vector<Entry> entry_;
  std::vector<Vertex> touched_;
  std::vector<Vertex> settled_;
  std::vector<Vertex> scratch_;
  IndexedBinaryHeap<W> heap_;
  Stamp cur_ = 0;
  std::uint64_t full_resets_ = 0;
};

namespace detail {

template <class W>
void check_search_args(const CsrGraph<W>& g, Vertex s, W radius) {
  if (s >= g.num_vertices()) {
    throw std::out_of_range("dijkstra: source " + std::to_string(s) + " out of range");
  }
  if constexpr (std::is_floating_point_v<W>) {
    if (std::isnan(radius)) throw std::invalid_argument("dijkstra: radius is NaN");
  }
}

inline void count(OpCounters* ops, [[maybe_unused]] std::uint64_t OpCounters::* field) noexcept {
  if (ops) SSSP_COUNT(++(ops->*field));
}

// Sec 3.1: "offer that value plus each outgoing arc length to the arc's head". Calls
// on_discover for a head touched for the first time and on_improve for a lowered tentative value.
template <class W, class Stamp, class OnDiscover, class OnImprove>
void offer_arcs(const CsrGraph<W>& g, Vertex v, W d, SearchArena<W, Stamp>& arena, OpCounters* ops,
                OnDiscover on_discover, OnImprove on_improve) {
  count(ops, &OpCounters::vertex_scans);
  for (const auto& a : g.arcs(v)) {
    count(ops, &OpCounters::edge_scans);
    count(ops, &OpCounters::additions);
    const W nd = WeightTraits<W>::checked_add(d, a.weight);
    count(ops, &OpCounters::comparisons);
    if (!arena.touched(a.to)) {  // tentative value infinity: every finite offer improves it
      arena.discover(a.to, nd, v, a.edge);
      count(ops, &OpCounters::relaxations);
      on_discover(a.to, nd);
    } else if (nd < arena.dist(a.to)) {  // never true for a settled head (non-negative lengths)
      arena.improve(a.to, nd, v, a.edge);
      count(ops, &OpCounters::relaxations);
      on_improve(a.to, nd);
    }
  }
}

}  // namespace detail

/// Lemma 3.2, binary-heap version, on the reusable arena. Settles exactly the vertices reachable
/// from s at distance <= radius (radius = infinity: all reachable vertices). Returns the number
/// of settled vertices; distances, parents and the settled order are read from the arena.
/// Throws std::out_of_range for an invalid source, std::invalid_argument for a NaN radius and,
/// std::overflow_error if an offered sum would reach infinity (integer lengths) or round to
/// infinity (floating-point lengths).
template <class W, class Stamp>
std::size_t dijkstra_search(const CsrGraph<W>& g, Vertex s, W radius, SearchArena<W, Stamp>& arena,
                            OpCounters* ops = nullptr) {
  detail::check_search_args(g, s, radius);
  arena.begin_search(g.num_vertices());
  auto& heap = arena.heap();
  heap.set_counters(ops);
  arena.discover(s, W{0}, kNoVertex, kNoEdge);
  heap.push(s, W{0});
  detail::count(ops, &OpCounters::queue_pushes);
  while (!heap.empty()) {
    // Extract a vertex of smallest unsettled value.
    const auto [d, v] = heap.pop_min();
    detail::count(ops, &OpCounters::queue_pops);
    // Stop when the smallest finite value exceeds the radius (d == radius is settled).
    detail::count(ops, &OpCounters::comparisons);
    if (radius < d) break;
    arena.settle(v);
    detail::offer_arcs(
        g, v, d, arena, ops,
        [&](Vertex x, W nd) {
          heap.push(x, nd);
          detail::count(ops, &OpCounters::queue_pushes);
        },
        [&](Vertex x, W nd) {
          heap.decrease(x, nd);
          detail::count(ops, &OpCounters::decrease_keys);
        });
  }
  heap.set_counters(nullptr);
  return arena.settled_order().size();
}

/// Lemma 3.2, scan version: each minimum is found by scanning the discovered unsettled vertices
/// (at most k), so a full search costs O(k^2 + l). Same contract as dijkstra_search.
template <class W, class Stamp>
std::size_t dijkstra_scan_search(const CsrGraph<W>& g, Vertex s, W radius,
                                 SearchArena<W, Stamp>& arena, OpCounters* ops = nullptr) {
  detail::check_search_args(g, s, radius);
  arena.begin_search(g.num_vertices());
  auto& frontier = arena.scratch();
  arena.discover(s, W{0}, kNoVertex, kNoEdge);
  frontier.push_back(s);
  detail::count(ops, &OpCounters::queue_pushes);
  while (!frontier.empty()) {
    // A scan for the minimum: |frontier| - 1 comparisons.
    std::size_t best = 0;
    for (std::size_t i = 1; i < frontier.size(); ++i) {
      detail::count(ops, &OpCounters::comparisons);
      if (arena.dist(frontier[i]) < arena.dist(frontier[best])) best = i;
    }
    const Vertex v = frontier[best];
    const W d = arena.dist(v);
    detail::count(ops, &OpCounters::comparisons);
    if (radius < d) break;  // the smallest finite value exceeds the radius
    frontier[best] = frontier.back();
    frontier.pop_back();
    detail::count(ops, &OpCounters::queue_pops);
    arena.settle(v);
    detail::offer_arcs(
        g, v, d, arena, ops,
        [&](Vertex x, W) {
          frontier.push_back(x);
          detail::count(ops, &OpCounters::queue_pushes);
        },
        [&](Vertex, W) { detail::count(ops, &OpCounters::decrease_keys); });
  }
  return arena.settled_order().size();
}

namespace detail {

template <class W, class Stamp>
SsspResult<W> export_result(std::size_t n, const SearchArena<W, Stamp>& arena, SsspResult<W> r) {
  r.dist.assign(n, WeightTraits<W>::infinity());
  r.parent.assign(n, kNoVertex);
  r.parent_edge.assign(n, kNoEdge);
  for (const Vertex v : arena.settled_order()) {
    r.dist[v] = arena.dist(v);
    r.parent[v] = arena.parent(v);
    r.parent_edge[v] = arena.parent_edge(v);
  }
  return r;
}

}  // namespace detail

/// Convenience form of dijkstra_search with a fresh arena: vertices not settled (unreachable or
/// beyond the radius) get infinity and kNoVertex. Operation counts go to diag.ops.
template <class W>
SsspResult<W> dijkstra(const CsrGraph<W>& g, Vertex s, W radius = WeightTraits<W>::infinity()) {
  SearchArena<W> arena;
  SsspResult<W> r;
  dijkstra_search(g, s, radius, arena, &r.diag.ops);
  return detail::export_result(g.num_vertices(), arena, std::move(r));
}

/// Convenience form of dijkstra_scan_search (Lemma 3.2 scan version, O(k^2 + l)).
template <class W>
SsspResult<W> dijkstra_scan(const CsrGraph<W>& g, Vertex s,
                            W radius = WeightTraits<W>::infinity()) {
  SearchArena<W> arena;
  SsspResult<W> r;
  dijkstra_scan_search(g, s, radius, arena, &r.diag.ops);
  return detail::export_result(g.num_vertices(), arena, std::move(r));
}

extern template class IndexedBinaryHeap<std::uint64_t>;
extern template class IndexedBinaryHeap<double>;
extern template class SearchArena<std::uint64_t>;
extern template class SearchArena<double>;
extern template SsspResult<std::uint64_t> dijkstra<std::uint64_t>(const CsrGraph<std::uint64_t>&,
                                                                  Vertex, std::uint64_t);
extern template SsspResult<double> dijkstra<double>(const CsrGraph<double>&, Vertex, double);
extern template SsspResult<std::uint64_t> dijkstra_scan<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, Vertex, std::uint64_t);
extern template SsspResult<double> dijkstra_scan<double>(const CsrGraph<double>&, Vertex, double);

}  // namespace sssp::paper
