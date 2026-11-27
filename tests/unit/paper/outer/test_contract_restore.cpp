// Sec 8.3, Algorithm 8 l.2 (contract zero-edge components, discard loops, retain the source
// component, keep original edge identifiers) and Lemma 8.5 (restore all original distances and a
// tree of first-discovery tight original arcs; unreachable vertices get infinity; the output
// distances are sums of original weights down the tree), arXiv 2611.34911v1, txt l.2145-2185.
#include <array>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "paper/outer/contract_restore.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::outer::contract_zero_edges;
using sssp::paper::outer::restore_predecessor_tree;
using sssp::paper::outer::ZeroContraction;

namespace {

template <class W>
CsrGraph<W> make_graph(std::size_t n, const std::vector<std::array<std::uint64_t, 3>>& es) {
  std::vector<Edge<W>> edges;
  edges.reserve(es.size());
  for (const auto& e : es) {
    edges.push_back(
        Edge<W>{static_cast<Vertex>(e[0]), static_cast<Vertex>(e[1]), static_cast<W>(e[2])});
  }
  return CsrGraph<W>::from_edges(n, std::move(edges));
}

/// Test-local union-find over zero-weight edges (independent of the implementation).
struct NaiveDsu {
  std::vector<std::size_t> p;
  explicit NaiveDsu(std::size_t n) : p(n) { std::iota(p.begin(), p.end(), std::size_t{0}); }
  std::size_t find(std::size_t x) {
    while (p[x] != x) x = p[x];
    return x;
  }
  void unite(std::size_t a, std::size_t b) { p[find(a)] = find(b); }
};

/// Structural checks of Alg 8 l.2 that hold for every input.
template <class W>
void check_contraction_structure(const CsrGraph<W>& g, Vertex s, const ZeroContraction<W>& c) {
  const std::size_t n = g.num_vertices();
  REQUIRE(c.vertex_map.size() == n);
  REQUIRE(c.zero_parent.size() == n);
  REQUIRE(c.zero_parent_edge.size() == n);
  const std::size_t big_n = c.graph.num_vertices();
  REQUIRE(c.representative.size() == big_n);
  REQUIRE(c.original_edge.size() == c.graph.num_edges());
  REQUIRE(c.source < big_n);
  REQUIRE(c.vertex_map[s] == c.source);
  REQUIRE(c.representative[c.source] == s);

  // Reference: zero components and the source's connected component, by naive union-find.
  NaiveDsu zero(n);
  NaiveDsu any(n);
  for (const auto& e : g.edges()) {
    any.unite(e.u, e.v);
    if (e.weight == W{0}) zero.unite(e.u, e.v);
  }
  for (Vertex v = 0; v < n; ++v) {
    const bool in_source_component = any.find(v) == any.find(s);
    if (!in_source_component) {
      CHECK(c.vertex_map[v] == kNoVertex);
      CHECK(c.zero_parent[v] == kNoVertex);
      CHECK(c.zero_parent_edge[v] == kNoEdge);
      continue;
    }
    REQUIRE(c.vertex_map[v] < big_n);
    // Representatives belong to their own class.
    CHECK(c.vertex_map[c.representative[c.vertex_map[v]]] == c.vertex_map[v]);
    for (Vertex u = 0; u < v; ++u) {
      if (any.find(u) != any.find(s)) continue;
      CHECK((c.vertex_map[u] == c.vertex_map[v]) == (zero.find(u) == zero.find(v)));
    }
    // Zero spanning tree inside the class: parent edges have weight 0 and lead to the
    // representative.
    if (c.representative[c.vertex_map[v]] == v) {
      CHECK(c.zero_parent[v] == kNoVertex);
      CHECK(c.zero_parent_edge[v] == kNoEdge);
    } else {
      Vertex x = v;
      std::size_t steps = 0;
      while (x != c.representative[c.vertex_map[v]] && steps <= n) {
        const EdgeId pe = c.zero_parent_edge[x];
        REQUIRE(pe < g.num_edges());
        const auto& e = g.edge(pe);
        CHECK(e.weight == W{0});
        const Vertex p = c.zero_parent[x];
        CHECK(((e.u == x && e.v == p) || (e.v == x && e.u == p)));
        CHECK(c.vertex_map[p] == c.vertex_map[v]);
        x = p;
        ++steps;
      }
      CHECK(x == c.representative[c.vertex_map[v]]);
    }
  }
  // Every contracted vertex is the class of some original vertex.
  std::vector<bool> hit(big_n, false);
  for (Vertex v = 0; v < n; ++v) {
    if (c.vertex_map[v] != kNoVertex) hit[c.vertex_map[v]] = true;
  }
  for (std::size_t x = 0; x < big_n; ++x) CHECK(hit[x]);

  // Edges: exactly the original edges of the source component whose endpoints lie in different
  // classes, in original-id order, with the original weight (positive) and mapped endpoints.
  std::vector<EdgeId> expected;
  for (EdgeId id = 0; id < g.num_edges(); ++id) {
    const auto& e = g.edge(id);
    if (c.vertex_map[e.u] == kNoVertex) continue;
    if (c.vertex_map[e.u] != c.vertex_map[e.v]) expected.push_back(id);
  }
  REQUIRE(c.original_edge == expected);
  for (EdgeId i = 0; i < c.graph.num_edges(); ++i) {
    const auto& ce = c.graph.edge(i);
    const auto& oe = g.edge(c.original_edge[i]);
    CHECK(ce.u != ce.v);
    CHECK(ce.weight > W{0});
    CHECK(ce.weight == oe.weight);
    CHECK(ce.u == c.vertex_map[oe.u]);
    CHECK(ce.v == c.vertex_map[oe.v]);
  }
}

/// Contract, solve on the contracted graph with a baseline Dijkstra, restore, verify.
template <class W>
sssp::SsspResult<W> contract_solve_restore(const CsrGraph<W>& g, Vertex s) {
  const auto c = contract_zero_edges(g, s);
  check_contraction_structure(g, s, c);
  const auto r = sssp::baseline::dijkstra_binary_heap(c.graph, c.source);
  auto res = restore_predecessor_tree(g, s, c, r.dist);
  const auto bf = sssp::baseline::bellman_ford(g, s);
  const auto rep = sssp::verify_sssp(g, s, res, &bf.dist);
  INFO(rep.summary());
  REQUIRE(rep.ok);
  REQUIRE(res.dist == bf.dist);
  return res;
}

// n = 7, source 2. Zero triangle {0,1,2} with a zero self-loop, zero edge {3,4}, parallel edges
// between the two classes, positive self-loops, and an unreachable edge {5,6}.
template <class W>
CsrGraph<W> mixed_example() {
  return make_graph<W>(7, {{0, 1, 0},
                           {1, 2, 0},
                           {2, 0, 0},
                           {1, 3, 5},
                           {0, 3, 7},
                           {3, 3, 4},
                           {3, 4, 0},
                           {4, 1, 2},
                           {2, 2, 0},
                           {5, 6, 1},
                           {0, 0, 3}});
}

}  // namespace

TEMPLATE_TEST_CASE("contraction without zero edges is the identity up to renaming",
                   "[paper][outer][contract]", std::uint64_t, double) {
  using W = TestType;
  // Connected, positive weights, with a parallel pair {1,2}.
  const auto g = make_graph<W>(
      5, {{0, 1, 3}, {1, 2, 1}, {2, 1, 4}, {2, 3, 2}, {3, 4, 6}, {4, 0, 9}, {1, 3, 1}});
  const auto c = contract_zero_edges(g, 0);
  check_contraction_structure(g, 0, c);
  REQUIRE(c.graph.num_vertices() == 5);
  REQUIRE(c.graph.num_edges() == g.num_edges());
  std::vector<EdgeId> all(g.num_edges());
  std::iota(all.begin(), all.end(), EdgeId{0});
  REQUIRE(c.original_edge == all);
  for (Vertex v = 0; v < 5; ++v) {
    CHECK(c.representative[c.vertex_map[v]] == v);
    CHECK(c.zero_parent[v] == kNoVertex);
  }
  contract_solve_restore(g, 0);
  contract_solve_restore(g, 3);
}

TEMPLATE_TEST_CASE("all-zero graph contracts to a single vertex", "[paper][outer][contract]",
                   std::uint64_t, double) {
  using W = TestType;
  // K4 with zero weights, a zero self-loop and a zero parallel edge.
  const auto g = make_graph<W>(
      4, {{0, 1, 0}, {0, 2, 0}, {0, 3, 0}, {1, 2, 0}, {1, 3, 0}, {2, 3, 0}, {2, 2, 0}, {3, 1, 0}});
  const auto c = contract_zero_edges(g, 1);
  check_contraction_structure(g, 1, c);
  REQUIRE(c.graph.num_vertices() == 1);
  REQUIRE(c.graph.num_edges() == 0);
  const auto res = contract_solve_restore(g, 1);
  REQUIRE(res.dist == std::vector<W>(4, W{0}));
}

TEMPLATE_TEST_CASE("zero cycles are contracted and the restored tree has no cycle",
                   "[paper][outer][contract]", std::uint64_t, double) {
  using W = TestType;
  // Zero 5-cycle 0..4, pendant 5 at weight 3 from vertex 2, a second zero 3-cycle {5,6,7},
  // and a positive chord back to 0.
  const auto g = make_graph<W>(8, {{0, 1, 0},
                                   {1, 2, 0},
                                   {2, 3, 0},
                                   {3, 4, 0},
                                   {4, 0, 0},
                                   {2, 5, 3},
                                   {5, 6, 0},
                                   {6, 7, 0},
                                   {7, 5, 0},
                                   {7, 0, 8}});
  const auto c = contract_zero_edges(g, 0);
  check_contraction_structure(g, 0, c);
  REQUIRE(c.graph.num_vertices() == 2);
  REQUIRE(c.graph.num_edges() == 2);
  const auto res = contract_solve_restore(g, 0);
  REQUIRE(res.dist == std::vector<W>{0, 0, 0, 0, 0, 3, 3, 3});
  // The only tight arc into the second class is {2,5}.
  std::size_t into_second = 0;
  for (Vertex v = 5; v < 8; ++v) into_second += res.parent[v] < 5 ? 1 : 0;
  REQUIRE(into_second == 1);
  REQUIRE(res.parent[5] == 2);
  REQUIRE(res.parent_edge[5] == 5);
}

TEMPLATE_TEST_CASE("source inside a zero component", "[paper][outer][contract]", std::uint64_t,
                   double) {
  using W = TestType;
  // Zero path 1-2-3; the source 2 is an interior member; exits leave from 1 and from 3.
  const auto g = make_graph<W>(5, {{0, 1, 4}, {1, 2, 0}, {2, 3, 0}, {3, 4, 1}, {0, 4, 2}});
  const auto c = contract_zero_edges(g, 2);
  check_contraction_structure(g, 2, c);
  REQUIRE(c.graph.num_vertices() == 3);
  REQUIRE(c.representative[c.source] == 2);
  const auto res = contract_solve_restore(g, 2);
  REQUIRE(res.dist == std::vector<W>{3, 0, 0, 0, 1});
  REQUIRE(res.parent[2] == kNoVertex);
  REQUIRE(res.parent[0] == 4);
}

TEMPLATE_TEST_CASE("only the source component is retained", "[paper][outer][contract]",
                   std::uint64_t, double) {
  using W = TestType;
  const W inf = WeightTraits<W>::infinity();
  // Two components {0,1,2} and {3,4,5} plus isolated vertex 6; zero edges in both.
  const auto g = make_graph<W>(7, {{0, 1, 0}, {1, 2, 2}, {3, 4, 0}, {4, 5, 0}, {5, 3, 1}});
  const auto c = contract_zero_edges(g, 1);
  check_contraction_structure(g, 1, c);
  REQUIRE(c.graph.num_vertices() == 2);
  REQUIRE(c.graph.num_edges() == 1);
  const auto res = contract_solve_restore(g, 1);
  REQUIRE(res.dist == std::vector<W>{0, 0, 2, inf, inf, inf, inf});
  for (Vertex v = 3; v < 7; ++v) {
    CHECK(res.parent[v] == kNoVertex);
    CHECK(res.parent_edge[v] == kNoEdge);
  }

  // An isolated source retains a single vertex and no edge.
  const auto ci = contract_zero_edges(g, 6);
  check_contraction_structure(g, 6, ci);
  REQUIRE(ci.graph.num_vertices() == 1);
  REQUIRE(ci.graph.num_edges() == 0);
  const auto resi = contract_solve_restore(g, 6);
  REQUIRE(resi.dist == std::vector<W>{inf, inf, inf, inf, inf, inf, 0});
}

TEMPLATE_TEST_CASE("parallel edges are kept, every loop is discarded", "[paper][outer][contract]",
                   std::uint64_t, double) {
  using W = TestType;
  const auto g = mixed_example<W>();
  const auto c = contract_zero_edges(g, 2);
  check_contraction_structure(g, 2, c);
  REQUIRE(c.graph.num_vertices() == 2);
  // e3 {1,3,5}, e4 {0,3,7} and e7 {4,1,2} all join the classes {0,1,2} and {3,4}.
  REQUIRE(c.original_edge == std::vector<EdgeId>{3, 4, 7});
  const auto res = contract_solve_restore(g, 2);
  const W inf = WeightTraits<W>::infinity();
  REQUIRE(res.dist == std::vector<W>{0, 0, 0, 2, 2, inf, inf});
  // Neither e3 nor e4 is tight (0 + 5 != 2, 0 + 7 != 2): class {3,4} is entered via e7 at 4.
  REQUIRE(res.parent[4] == 1);
  REQUIRE(res.parent_edge[4] == 7);
  REQUIRE(res.parent[3] == 4);
  REQUIRE(res.parent_edge[3] == 6);
}

TEMPLATE_TEST_CASE("contraction counters report removed vertices and edges",
                   "[paper][outer][contract]", std::uint64_t, double) {
  using W = TestType;
  const auto g = mixed_example<W>();
  sssp::Diagnostics diag;
  const auto c = contract_zero_edges(g, 2, &diag);
  REQUIRE(c.graph.num_vertices() == 2);
  const auto& k = diag.counters;
  CHECK(k.at("contract_input_vertices") == 7);
  CHECK(k.at("contract_input_edges") == 11);
  CHECK(k.at("contract_retained_vertices") == 2);
  CHECK(k.at("contract_retained_edges") == 3);
  CHECK(k.at("contract_merged_vertices") == 3);     // {0,1,2} -> 1 and {3,4} -> 1
  CHECK(k.at("contract_unreached_vertices") == 2);  // 5, 6
  CHECK(k.at("contract_removed_vertices") == 5);
  CHECK(k.at("contract_dropped_loops") == 7);    // e0 e1 e2 e5 e6 e8 e10
  CHECK(k.at("contract_unreached_edges") == 1);  // e9
  CHECK(k.at("contract_removed_edges") == 8);
  CHECK(k.at("contract_zero_edges") == 5);  // e0 e1 e2 e6 e8
}

TEMPLATE_TEST_CASE("restore reports tree and unreachable counts", "[paper][outer][contract]",
                   std::uint64_t, double) {
  using W = TestType;
  const auto g = mixed_example<W>();
  const auto c = contract_zero_edges(g, 2);
  const auto res = restore_predecessor_tree(g, 2, c, std::vector<W>{0, 2});
  CHECK(res.diag.counters.at("restore_tree_vertices") == 5);
  CHECK(res.diag.counters.at("restore_unreachable_vertices") == 2);
  CHECK(res.diag.fallbacks.empty());
}

TEMPLATE_TEST_CASE("restore rejects inputs that are not exact contracted distances",
                   "[paper][outer][contract]", std::uint64_t, double) {
  using W = TestType;
  const auto g = mixed_example<W>();
  const auto c = contract_zero_edges(g, 2);
  // Wrong length.
  REQUIRE_THROWS_AS(restore_predecessor_tree(g, 2, c, std::vector<W>{0}), std::invalid_argument);
  // Source label not zero.
  REQUIRE_THROWS_AS(restore_predecessor_tree(g, 2, c, std::vector<W>{1, 2}), std::invalid_argument);
  // Source does not match the contraction.
  REQUIRE_THROWS_AS(restore_predecessor_tree(g, 0, c, std::vector<W>{0, 2}), std::invalid_argument);
  // Label 3 for class {3,4} is not attained by any tight arc (0+5, 0+7, 0+2): not exact.
  REQUIRE_THROWS_AS(restore_predecessor_tree(g, 2, c, std::vector<W>{0, 3}), std::invalid_argument);
  // Contraction of another graph.
  const auto other = make_graph<W>(3, {{0, 1, 1}});
  REQUIRE_THROWS_AS(restore_predecessor_tree(other, 2, c, std::vector<W>{0, 2}),
                    std::invalid_argument);
}

TEMPLATE_TEST_CASE("contraction rejects an invalid source", "[paper][outer][contract]",
                   std::uint64_t, double) {
  using W = TestType;
  const auto g = make_graph<W>(3, {{0, 1, 0}});
  REQUIRE_THROWS_AS(contract_zero_edges(g, 3), std::out_of_range);
  const CsrGraph<W> empty = CsrGraph<W>::from_edges(0, {});
  REQUIRE_THROWS_AS(contract_zero_edges(empty, 0), std::out_of_range);
}

TEST_CASE("integer weights at the infinity sentinel are never tight", "[paper][outer][contract]") {
  using W = std::uint64_t;
  const W big = WeightTraits<W>::infinity();
  // Edge {0,1} has weight 2^64-1 (a valid input weight, but no finite label can use it), the
  // path 0-2-1 has length 2^63 + (2^63 - 2) = 2^64 - 2.
  const W half = W{1} << 63;
  const auto g = CsrGraph<W>::from_edges(
      3, {Edge<W>{0, 1, big}, Edge<W>{0, 2, half}, Edge<W>{2, 1, half - 2}});
  // The baselines reject a weight equal to their infinity, so the labels are given directly.
  const auto c = contract_zero_edges(g, 0);
  check_contraction_structure(g, 0, c);
  std::vector<W> labels(3);
  labels[c.vertex_map[0]] = 0;
  labels[c.vertex_map[1]] = big - 1;
  labels[c.vertex_map[2]] = half;
  const auto res = restore_predecessor_tree(g, 0, c, labels);
  REQUIRE(res.dist == std::vector<W>{0, big - 1, half});
  REQUIRE(res.parent[1] == 2);
  REQUIRE(res.parent_edge[1] == 2);
}
