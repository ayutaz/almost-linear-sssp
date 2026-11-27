#pragma once
// Test-only reference values for Algorithm 3 (Evaluate) on a seeded ordered hierarchy
// (Def 4.2). Nothing here is part of the paper's algorithm and nothing here is shared with an
// Evaluate implementation: the oracles below are written from the definitions only.
//
//   (a) target distance  delta(v) = min_i (c_i + d(v_i, v))     (Def 4.2)
//       by a multi-source binary-heap Dijkstra on the undirected physical graph;
//   (b) protected states (Def 4.2, "tight witness"): a breadth-first search that starts at
//       the states whose least initial offer equals delta(pi(x)) and follows only the state
//       arcs y -> x with delta(pi(y)) + l(y -> x) = delta(pi(x)) (research notes 02 sec 2.2);
//   (c) the multi-source Dijkstra on the DIRECTED state graph D seeded with the least initial
//       offer of every state, i.e. dist_D(x) = min_y (q(y) + d_D(y, x)).
//
// How the three relate (used by the self-tests and by the differential tests):
//   * Under Def 4.1 (H1: every arc projects to a physical walk no longer than itself) and
//     feasibility of all initial offers, every state path gives q(y) + len >= delta(pi(x)), so
//     dist_D(x) >= delta(pi(x)) for every state, and x is protected iff dist_D(x) ==
//     delta(pi(x)). The BFS in (b) and this Dijkstra test are two independent computations of
//     the same set; test_eval_oracle.cpp checks that they agree on valid instances.
//   * Lemma 4.6 promises delta(pi(x)) only at PROTECTED states. An unprotected state may emit a
//     value strictly below dist_D(x) (supplemental offers use physical gap walks, not state
//     arcs) or may never emit at all (it can retire unemitted). Hence (c) is the reference for
//     protected states only; for every state the reference is the lower bound delta(pi(x)).
//     (research notes 02 A9 / open question Q-29.)
//
// Weights: exact for unsigned integers (overflow throws std::overflow_error via
// WeightTraits::checked_add). For double the equality tests in (b) are exact IEEE comparisons;
// use integer-valued doubles below 2^53 so that every sum is exact.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/graph.hpp"

namespace sssp::test::paper {

using sssp::paper::PhysicalSource;
using sssp::paper::SeededOrderedHierarchy;
using sssp::paper::StateArc;
using sssp::paper::StateGraph;
using sssp::paper::StateId;
using sssp::paper::StateOffer;

/// Every reference value for one seeded ordered hierarchy.
template <class W>
struct EvalReference {
  std::vector<W> target;           ///< delta(v) per physical vertex (infinity if unreachable)
  std::vector<W> state_target;     ///< delta(pi(x)) per state
  std::vector<W> best_offer;       ///< least initial offer per state (infinity if none)
  std::vector<bool> is_protected;  ///< (b) protected states (tight-arc BFS)
  std::vector<W> state_distance;   ///< (c) dist_D(x) from the initial offers
  /// States with a finite dist_D in nondecreasing (dist_D, state id) order: a valid
  /// emission order for a plain state-graph Dijkstra.
  std::vector<StateId> dijkstra_order;

  [[nodiscard]] std::size_t num_protected() const noexcept {
    std::size_t k = 0;
    for (const bool p : is_protected) k += p ? 1 : 0;
    return k;
  }
  /// The value every protected state must emit; for an unprotected state the lower bound.
  [[nodiscard]] W expected_emission(StateId x) const { return state_target.at(x); }
};

// ------------------------------------------------------------------ (a) target distances

/// delta(v) = min_i (c_i + d(v_i, v)) on the undirected physical graph `g`.
/// Throws std::out_of_range for a source vertex outside g, std::invalid_argument for an
/// invalid offset (negative, NaN, infinite).
template <class W>
std::vector<W> target_distances(const CsrGraph<W>& g, std::span<const PhysicalSource<W>> sources) {
  using T = WeightTraits<W>;
  const std::size_t n = g.num_vertices();
  std::vector<W> dist(n, T::infinity());
  using Item = std::pair<W, Vertex>;
  std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
  for (const auto& s : sources) {
    if (s.vertex >= n) throw std::out_of_range("physical source vertex out of range");
    if (!T::is_valid_weight(s.offset) || s.offset == T::infinity()) {
      throw std::invalid_argument("physical source offset must be finite and nonnegative");
    }
    if (s.offset < dist[s.vertex]) {
      dist[s.vertex] = s.offset;
      heap.emplace(s.offset, s.vertex);
    }
  }
  while (!heap.empty()) {
    const auto [d, u] = heap.top();
    heap.pop();
    if (d != dist[u]) continue;
    for (const auto& a : g.arcs(u)) {
      const W nd = T::checked_add(d, a.weight);
      if (nd < dist[a.to]) {
        dist[a.to] = nd;
        heap.emplace(nd, a.to);
      }
    }
  }
  return dist;
}

/// Least initial offer per state (infinity if none). Throws std::out_of_range for an offer at a
/// non-existent state and std::invalid_argument for a non-finite or negative offer value
/// (Def 4.2: an offer is a finite value; +infinity is "no offer" and is not stored).
template <class W>
std::vector<W> best_initial_offers(std::size_t num_states, std::span<const StateOffer<W>> offers) {
  using T = WeightTraits<W>;
  std::vector<W> best(num_states, T::infinity());
  for (const auto& o : offers) {
    if (o.state >= num_states) throw std::out_of_range("initial offer at a non-existent state");
    if (!T::is_valid_weight(o.value) || o.value == T::infinity()) {
      throw std::invalid_argument("initial offer must be a finite nonnegative value");
    }
    if (o.value < best[o.state]) best[o.state] = o.value;
  }
  return best;
}

// ------------------------------------------------------------------ (b) protected states

/// Protected states (Def 4.2) given delta(pi(x)) per state and the least initial offer per
/// state: BFS from {x : best_offer[x] == state_target[x] < infinity} along the arcs y -> x
/// with state_target[y] + l == state_target[x]. Every prefix of a tight witness is tight
/// (Def 4.2, l.670), which is why following tight arcs only is complete under H1.
template <class W>
std::vector<bool> protected_states(const StateGraph<W>& graph, const std::vector<W>& state_target,
                                   const std::vector<W>& best_offer) {
  using T = WeightTraits<W>;
  const std::size_t n = graph.num_states;
  if (state_target.size() != n || best_offer.size() != n) {
    throw std::invalid_argument("protected_states: per-state vectors have the wrong size");
  }
  std::vector<std::vector<std::size_t>> out(n);
  for (std::size_t i = 0; i < graph.arcs.size(); ++i) {
    const auto& a = graph.arcs[i];
    if (a.tail >= n || a.head >= n) throw std::out_of_range("state arc endpoint out of range");
    out[a.tail].push_back(i);
  }
  std::vector<bool> prot(n, false);
  std::queue<StateId> q;
  for (StateId x = 0; x < n; ++x) {
    if (best_offer[x] != T::infinity() && best_offer[x] == state_target[x]) {
      prot[x] = true;
      q.push(x);
    }
  }
  while (!q.empty()) {
    const StateId y = q.front();
    q.pop();
    for (const std::size_t i : out[y]) {
      const auto& a = graph.arcs[i];
      if (prot[a.head] || state_target[a.head] == T::infinity()) continue;
      // state_target[y] is finite (y is protected). Saturate instead of throwing: a sum that
      // would overflow cannot equal a finite target.
      if (state_target[y] >= T::infinity() - a.length) continue;
      if (T::checked_add(state_target[y], a.length) == state_target[a.head]) {
        prot[a.head] = true;
        q.push(a.head);
      }
    }
  }
  return prot;
}

// ------------------------------------------------------------------ (c) state-graph Dijkstra

/// Multi-source Dijkstra on the directed state graph from the per-state initial values
/// `initial` (infinity = no offer). Returns dist_D and fills `order` (if non-null) with the
/// states of finite distance in nondecreasing (distance, id) order.
template <class W>
std::vector<W> state_graph_distances(const StateGraph<W>& graph, const std::vector<W>& initial,
                                     std::vector<StateId>* order = nullptr) {
  using T = WeightTraits<W>;
  const std::size_t n = graph.num_states;
  if (initial.size() != n) throw std::invalid_argument("state_graph_distances: wrong size");
  std::vector<std::vector<std::size_t>> out(n);
  for (std::size_t i = 0; i < graph.arcs.size(); ++i) {
    const auto& a = graph.arcs[i];
    if (a.tail >= n || a.head >= n) throw std::out_of_range("state arc endpoint out of range");
    out[a.tail].push_back(i);
  }
  std::vector<W> dist = initial;
  std::vector<bool> done(n, false);
  using Item = std::pair<W, StateId>;
  std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
  for (StateId x = 0; x < n; ++x) {
    if (dist[x] != T::infinity()) heap.emplace(dist[x], x);
  }
  if (order != nullptr) order->clear();
  while (!heap.empty()) {
    const auto [d, x] = heap.top();
    heap.pop();
    if (done[x] || d != dist[x]) continue;
    done[x] = true;
    if (order != nullptr) order->push_back(x);
    for (const std::size_t i : out[x]) {
      const auto& a = graph.arcs[i];
      const W nd = T::checked_add(d, a.length);
      if (nd < dist[a.head]) {
        dist[a.head] = nd;
        heap.emplace(nd, a.head);
      }
    }
  }
  return dist;
}

// ------------------------------------------------------------------ everything at once

/// All reference values for `s` on the physical graph `physical`. Throws std::out_of_range
/// for out-of-range sources / offer states / projections and std::invalid_argument for an
/// infeasible initial offer (value < delta(pi(x))): the definitions of Def 4.2 assume
/// feasibility, so the oracle refuses such inputs instead of producing a meaningless set.
template <class W>
EvalReference<W> compute_eval_reference(const SeededOrderedHierarchy<W>& s,
                                        const CsrGraph<W>& physical) {
  const auto& h = s.hierarchy;
  const std::size_t n = h.num_states();
  if (h.projection.size() != n) throw std::invalid_argument("projection size != num_states");
  EvalReference<W> r;
  r.target = target_distances<W>(physical, s.sources);
  r.state_target.resize(n);
  for (StateId x = 0; x < n; ++x) {
    const Vertex v = h.projection[x];
    if (v >= r.target.size()) throw std::out_of_range("projection outside the physical graph");
    r.state_target[x] = r.target[v];
  }
  r.best_offer = best_initial_offers<W>(n, s.offers);
  for (StateId x = 0; x < n; ++x) {
    if (r.best_offer[x] != WeightTraits<W>::infinity() && r.best_offer[x] < r.state_target[x]) {
      throw std::invalid_argument("infeasible initial offer at state " + std::to_string(x));
    }
  }
  r.is_protected = protected_states<W>(h.graph, r.state_target, r.best_offer);
  r.state_distance = state_graph_distances<W>(h.graph, r.best_offer, &r.dijkstra_order);
  return r;
}

/// Minimum emitted value per physical vertex from a list of (state, value) emissions, the
/// quantity the caller of Evaluate extracts (Sec 8 l.2192). infinity where nothing emitted.
template <class W>
std::vector<W> physical_minima(std::size_t num_physical,
                               const std::vector<sssp::Vertex>& projection,
                               std::span<const std::pair<StateId, W>> emissions) {
  std::vector<W> m(num_physical, WeightTraits<W>::infinity());
  for (const auto& [x, v] : emissions) {
    if (x >= projection.size()) throw std::out_of_range("emission at a non-existent state");
    const Vertex p = projection[x];
    if (p >= num_physical) throw std::out_of_range("projection outside the physical graph");
    if (v < m[p]) m[p] = v;
  }
  return m;
}

/// Physical vertices that have at least one protected state (the "requested targets" of
/// Def 4.2 that this instance actually serves). Lemma 4.6 then forces the physical minimum of
/// the emissions to equal delta(v) at exactly these vertices.
template <class W>
std::vector<bool> served_targets(std::size_t num_physical,
                                 const std::vector<sssp::Vertex>& projection,
                                 const EvalReference<W>& ref) {
  std::vector<bool> served(num_physical, false);
  for (StateId x = 0; x < projection.size(); ++x) {
    if (ref.is_protected.at(x)) served.at(projection[x]) = true;
  }
  return served;
}

}  // namespace sssp::test::paper
