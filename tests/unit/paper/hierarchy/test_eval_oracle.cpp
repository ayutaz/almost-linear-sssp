// Self-tests of the Evaluate oracle (tests/support/paper/eval_oracle.hpp): the target distance
// of Def 4.2 against hand values and Floyd-Warshall, the protected set (tight witnesses) on
// boundary examples worked out from the definition, and the state-graph Dijkstra against a
// directed Floyd-Warshall written here.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "support/floyd_warshall.hpp"
#include "support/graphs.hpp"
#include "support/paper/eval_oracle.hpp"
#include "support/prng.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::PhysicalSource;
using sssp::paper::SeededOrderedHierarchy;
using sssp::paper::StateGraph;
using sssp::paper::StateId;
using sssp::paper::StateOffer;
using namespace sssp::test::paper;

namespace {

using U = std::uint64_t;
constexpr U kInf = std::numeric_limits<U>::max();

// Directed all-pairs closure on the state graph (independent of the oracle's Dijkstra).
std::vector<std::vector<U>> directed_fw(const StateGraph<U>& g) {
  const std::size_t n = g.num_states;
  std::vector<std::vector<U>> d(n, std::vector<U>(n, kInf));
  for (std::size_t i = 0; i < n; ++i) d[i][i] = 0;
  for (const auto& a : g.arcs) {
    if (a.length < d[a.tail][a.head]) d[a.tail][a.head] = a.length;
  }
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        if (d[i][k] != kInf && d[k][j] != kInf && d[i][k] + d[k][j] < d[i][j])
          d[i][j] = d[i][k] + d[k][j];
  return d;
}

// Seeded hierarchy over a hand-made state graph; the forest is irrelevant to the oracle, so a
// leaf forest suffices.
SeededOrderedHierarchy<U> seeded(StateGraph<U> g, std::vector<Vertex> proj,
                                 std::vector<PhysicalSource<U>> src,
                                 std::vector<StateOffer<U>> offers) {
  SeededOrderedHierarchy<U> s;
  s.hierarchy = sssp::paper::make_leaf_forest<U>(std::move(g), std::move(proj));
  sssp::paper::collect_roots(s.hierarchy);
  s.sources = std::move(src);
  s.offers = std::move(offers);
  return s;
}

}  // namespace

TEST_CASE("eval oracle: target distance by hand on a path with one and several sources",
          "[eval-oracle]") {
  // 0 -1- 1 -4- 2     3 (isolated)
  const auto g = CsrGraph<U>::from_edges(4, {{0, 1, 1}, {1, 2, 4}});
  const std::vector<PhysicalSource<U>> one{{0, 0}};
  CHECK(target_distances<U>(g, one) == std::vector<U>{0, 1, 5, kInf});
  const std::vector<PhysicalSource<U>> two{{0, 3}, {2, 0}};
  CHECK(target_distances<U>(g, two) == std::vector<U>{3, 4, 0, kInf});
  // A source whose offset is dominated by another source is irrelevant.
  const std::vector<PhysicalSource<U>> dominated{{0, 0}, {1, 7}};
  CHECK(target_distances<U>(g, dominated) == std::vector<U>{0, 1, 5, kInf});
  // The same vertex twice: the smaller offset wins.
  const std::vector<PhysicalSource<U>> twice{{2, 9}, {2, 2}};
  CHECK(target_distances<U>(g, twice) == std::vector<U>{7, 6, 2, kInf});
  CHECK(target_distances<U>(g, std::vector<PhysicalSource<U>>{}) == std::vector<U>(4, kInf));
}

TEST_CASE("eval oracle: target distance equals min_i c_i + Floyd-Warshall", "[eval-oracle]") {
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    sssp::test::SplitMix64 rng(seed);
    const std::size_t n = 1 + rng.below(12);
    const auto edges = sssp::test::random_edges<U>(n, rng.below(3 * n + 1), 0, 9, seed);
    const auto g = CsrGraph<U>::from_edges(n, edges);
    const auto fw = sssp::test::floyd_warshall<U>(n, edges);
    std::vector<PhysicalSource<U>> src;
    for (std::size_t i = 0, k = 1 + rng.below(3); i < k; ++i) {
      src.push_back({static_cast<Vertex>(rng.below(n)), rng.below(8)});
    }
    const auto delta = target_distances<U>(g, src);
    for (Vertex v = 0; v < n; ++v) {
      U want = kInf;
      for (const auto& s : src) {
        if (fw[s.vertex][v] != kInf) want = std::min(want, s.offset + fw[s.vertex][v]);
      }
      INFO("seed " << seed << " vertex " << v);
      CHECK(delta[v] == want);
    }
  }
}

TEST_CASE("eval oracle: invalid sources and offers are rejected", "[eval-oracle]") {
  const auto g = CsrGraph<U>::from_edges(2, {{0, 1, 1}});
  CHECK_THROWS_AS(target_distances<U>(g, std::vector<PhysicalSource<U>>{{2, 0}}),
                  std::out_of_range);
  CHECK_THROWS_AS(target_distances<U>(g, std::vector<PhysicalSource<U>>{{0, kInf}}),
                  std::invalid_argument);
  CHECK_THROWS_AS(best_initial_offers<U>(2, std::vector<StateOffer<U>>{{2, 0}}), std::out_of_range);
  CHECK_THROWS_AS(best_initial_offers<U>(2, std::vector<StateOffer<U>>{{0, kInf}}),
                  std::invalid_argument);
  // Infeasible initial offer (value < delta(pi(x))): Def 4.2 forbids it.
  StateGraph<U> d;
  d.num_states = 2;
  d.arcs = {{0, 1, 1}};
  const auto s = seeded(d, {0, 1}, {{0, 0}}, {{0, 0}, {1, 0}});
  CHECK_THROWS_AS(compute_eval_reference<U>(s, g), std::invalid_argument);
}

TEST_CASE("eval oracle: best initial offer is the minimum per state", "[eval-oracle]") {
  const std::vector<StateOffer<U>> offers{{1, 5}, {1, 3}, {2, 0}, {1, 4}};
  CHECK(best_initial_offers<U>(4, offers) == std::vector<U>{kInf, 3, 0, kInf});
}

TEST_CASE("eval oracle: protected set boundary cases from Def 4.2", "[eval-oracle]") {
  // Physical path 0 -1- 1 -4- 2, source (0, 0): delta = [0, 1, 5].
  const auto g = CsrGraph<U>::from_edges(3, {{0, 1, 1}, {1, 2, 4}});

  SECTION("tight offer and tight chain are protected; a reverse arc is not needed") {
    // states a=0 (pi 0), b=1 (pi 1), c=2 (pi 2); arcs a->b (1), b->c (4).
    StateGraph<U> d;
    d.num_states = 3;
    d.arcs = {{0, 1, 1}, {1, 2, 4}};
    const auto r = compute_eval_reference<U>(seeded(d, {0, 1, 2}, {{0, 0}}, {{0, 0}}), g);
    CHECK(r.state_target == std::vector<U>{0, 1, 5});
    CHECK(r.is_protected == std::vector<bool>{true, true, true});
    CHECK(r.state_distance == std::vector<U>{0, 1, 5});
    CHECK(r.num_protected() == 3);
  }
  SECTION("an offer strictly above delta does not protect its state") {
    StateGraph<U> d;
    d.num_states = 3;
    d.arcs = {{0, 1, 1}, {1, 2, 4}};
    const auto r = compute_eval_reference<U>(seeded(d, {0, 1, 2}, {{0, 0}}, {{0, 1}}), g);
    CHECK(r.is_protected == std::vector<bool>{false, false, false});
    CHECK(r.state_distance == std::vector<U>{1, 2, 6});
  }
  SECTION("a state arc longer than the physical distance breaks the witness") {
    // b->c has length 5 > d(1, 2) = 4: c is reached at 6 > delta 5.
    StateGraph<U> d;
    d.num_states = 3;
    d.arcs = {{0, 1, 1}, {1, 2, 5}};
    const auto r = compute_eval_reference<U>(seeded(d, {0, 1, 2}, {{0, 0}}, {{0, 0}}), g);
    CHECK(r.is_protected == std::vector<bool>{true, true, false});
    CHECK(r.state_distance[2] == 6);
  }
  SECTION("a state on a tight path is protected although its own offer is not tight") {
    StateGraph<U> d;
    d.num_states = 3;
    d.arcs = {{0, 1, 1}, {1, 2, 4}};
    const auto r = compute_eval_reference<U>(seeded(d, {0, 1, 2}, {{0, 0}}, {{0, 0}, {1, 3}}), g);
    CHECK(r.is_protected == std::vector<bool>{true, true, true});
  }
  SECTION("a tight offer in the middle protects its successors only") {
    StateGraph<U> d;
    d.num_states = 3;
    d.arcs = {{0, 1, 1}, {1, 2, 4}};
    const auto r = compute_eval_reference<U>(seeded(d, {0, 1, 2}, {{0, 0}}, {{1, 1}}), g);
    CHECK(r.is_protected == std::vector<bool>{false, true, true});
    CHECK(r.state_distance == std::vector<U>{kInf, 1, 5});
  }
  SECTION("two copies of one vertex: only the one on a directed tight path is protected") {
    // states a (pi 0), b1 (pi 1), b2 (pi 1), c (pi 2); arcs a->b1 (1), b2->c (4), b1->b2 (3).
    StateGraph<U> d;
    d.num_states = 4;
    d.arcs = {{0, 1, 1}, {2, 3, 4}, {1, 2, 3}};
    const auto r = compute_eval_reference<U>(seeded(d, {0, 1, 1, 2}, {{0, 0}}, {{0, 0}}), g);
    CHECK(r.is_protected == std::vector<bool>{true, true, false, false});
    CHECK(r.state_distance == std::vector<U>{0, 1, 4, 8});
  }
  SECTION("zero-length arcs between copies are tight") {
    StateGraph<U> d;
    d.num_states = 4;
    d.arcs = {{0, 1, 1}, {1, 2, 0}, {2, 3, 4}};
    const auto r = compute_eval_reference<U>(seeded(d, {0, 1, 1, 2}, {{0, 0}}, {{0, 0}}), g);
    CHECK(r.is_protected == std::vector<bool>{true, true, true, true});
  }
  SECTION("unreachable images are never protected and get no offers") {
    const auto g2 = CsrGraph<U>::from_edges(4, {{0, 1, 1}, {1, 2, 4}});
    StateGraph<U> d;
    d.num_states = 2;
    d.arcs = {{0, 1, 0}};
    const auto r = compute_eval_reference<U>(seeded(d, {3, 3}, {{0, 0}}, {}), g2);
    CHECK(r.state_target == std::vector<U>{kInf, kInf});
    CHECK(r.is_protected == std::vector<bool>{false, false});
  }
  SECTION("multiple weighted sources: protection follows the nearest source") {
    // sources (0, 3), (2, 0): delta = [3, 4, 0]. Offer 3 at a, 0 at c; arcs a->b, c->b.
    StateGraph<U> d;
    d.num_states = 3;
    d.arcs = {{0, 1, 1}, {2, 1, 4}};
    const auto r =
        compute_eval_reference<U>(seeded(d, {0, 1, 2}, {{0, 3}, {2, 0}}, {{0, 3}, {2, 0}}), g);
    CHECK(r.state_target == std::vector<U>{3, 4, 0});
    CHECK(r.is_protected == std::vector<bool>{true, true, true});
    CHECK(r.state_distance == std::vector<U>{3, 4, 0});
  }
}

TEST_CASE("eval oracle: state-graph Dijkstra equals directed Floyd-Warshall from the offers",
          "[eval-oracle]") {
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    sssp::test::SplitMix64 rng(seed);
    StateGraph<U> d;
    d.num_states = 1 + rng.below(12);
    for (std::size_t i = 0, m = rng.below(4 * d.num_states); i < m; ++i) {
      d.arcs.push_back({static_cast<StateId>(rng.below(d.num_states)),
                        static_cast<StateId>(rng.below(d.num_states)), rng.below(7)});
    }
    std::vector<U> init(d.num_states, kInf);
    for (std::size_t i = 0, k = 1 + rng.below(3); i < k; ++i) {
      const auto x = rng.below(d.num_states);
      init[x] = std::min(init[x], rng.below(5));
    }
    std::vector<StateId> order;
    const auto dist = state_graph_distances<U>(d, init, &order);
    const auto fw = directed_fw(d);
    std::size_t finite = 0;
    for (StateId x = 0; x < d.num_states; ++x) {
      U want = kInf;
      for (StateId y = 0; y < d.num_states; ++y) {
        if (init[y] != kInf && fw[y][x] != kInf) want = std::min(want, init[y] + fw[y][x]);
      }
      INFO("seed " << seed << " state " << x);
      CHECK(dist[x] == want);
      finite += want != kInf ? 1 : 0;
    }
    REQUIRE(order.size() == finite);
    for (std::size_t i = 1; i < order.size(); ++i) {
      CHECK(dist[order[i - 1]] <= dist[order[i]]);
    }
  }
}

TEST_CASE("eval oracle: tight-arc BFS equals {x : dist_D(x) = delta(pi(x))} on Def 4.1 inputs",
          "[eval-oracle]") {
  // Reference hierarchy instances satisfy H1, so the two characterisations must agree.
  std::size_t total_protected = 0;
  std::size_t total_unprotected = 0;
  for (std::uint64_t seed = 1; seed <= 300; ++seed) {
    sssp::test::SplitMix64 rng(seed);
    StateGraph<U> d;
    d.num_states = 2 + rng.below(14);
    for (std::size_t i = 0, m = rng.below(3 * d.num_states + 1); i < m; ++i) {
      d.arcs.push_back({static_cast<StateId>(rng.below(d.num_states)),
                        static_cast<StateId>(rng.below(d.num_states)), 1 + rng.below(6)});
    }
    const std::size_t np = 1 + rng.below(d.num_states);
    std::vector<Vertex> proj(d.num_states);
    for (std::size_t x = 0; x < d.num_states; ++x) {
      proj[x] = static_cast<Vertex>(x < np ? x : rng.below(np));
    }
    auto edges_g = sssp::paper::physical_graph_of_arcs<U>(d, proj, np);
    std::vector<Edge<U>> edges(edges_g.edges().begin(), edges_g.edges().end());
    for (std::size_t i = 0, k = rng.below(4); i < k; ++i) {  // physical shortcuts
      edges.push_back({static_cast<Vertex>(rng.below(np)), static_cast<Vertex>(rng.below(np)),
                       1 + rng.below(3)});
    }
    const auto g = CsrGraph<U>::from_edges(np, edges);
    std::vector<PhysicalSource<U>> src{{proj[0], 0}};
    std::vector<StateOffer<U>> offers{{0, 0}};
    const auto delta = target_distances<U>(g, src);
    for (std::size_t i = 0, k = rng.below(4); i < k; ++i) {
      const auto x = static_cast<StateId>(rng.below(d.num_states));
      if (delta[proj[x]] != kInf) offers.push_back({x, delta[proj[x]] + rng.below(2)});
    }
    const auto r = compute_eval_reference<U>(seeded(d, proj, src, offers), g);
    for (StateId x = 0; x < d.num_states; ++x) {
      INFO("seed " << seed << " state " << x);
      CHECK(r.state_distance[x] >= r.state_target[x]);
      const bool tight = r.state_distance[x] != kInf && r.state_distance[x] == r.state_target[x];
      CHECK(r.is_protected[x] == tight);
      (r.is_protected[x] ? total_protected : total_unprotected) += 1;
    }
  }
  // Both outcomes occur, so the comparison is not vacuous.
  CHECK(total_protected > 300);
  CHECK(total_unprotected > 100);
}

TEST_CASE("eval oracle: physical minima and served targets", "[eval-oracle]") {
  const std::vector<Vertex> proj{0, 1, 1, 2};
  const std::vector<std::pair<StateId, U>> em{{1, 7}, {2, 4}, {0, 0}};
  CHECK(physical_minima<U>(3, proj, em) == std::vector<U>{0, 4, kInf});
  EvalReference<U> r;
  r.is_protected = {false, false, true, false};
  CHECK(served_targets<U>(3, proj, r) == std::vector<bool>{false, true, false});
}
