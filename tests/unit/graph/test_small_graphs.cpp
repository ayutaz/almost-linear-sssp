#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"

using sssp::EdgeList;
using sssp::Vertex;

namespace {

template <class W>
const sssp::gen::NamedGraph<W>& find(const std::vector<sssp::gen::NamedGraph<W>>& suite,
                                     const std::string& name) {
  for (const auto& g : suite) {
    if (g.name == name) return g;
  }
  FAIL("missing small graph: " << name);
  return suite.front();
}

template <class W>
std::size_t components(const EdgeList<W>& g) {
  std::vector<std::size_t> p(g.num_vertices);
  std::iota(p.begin(), p.end(), 0);
  auto f = [&](std::size_t x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  };
  std::size_t c = g.num_vertices;
  for (const auto& e : g.edges) {
    const auto a = f(e.u), b = f(e.v);
    if (a != b) {
      p[a] = b;
      --c;
    }
  }
  return c;
}

template <class W>
std::vector<std::size_t> degrees(const EdgeList<W>& g) {
  std::vector<std::size_t> d(g.num_vertices, 0);
  for (const auto& e : g.edges) {
    ++d[e.u];
    ++d[e.v];
  }
  return d;
}

/// Test-local O(n^2) Dijkstra that also counts shortest paths (exact for integral weights).
struct Paths {
  std::vector<double> dist;
  std::vector<double> count;
  std::vector<std::size_t> hops;  // hop count along one shortest path tree
};

template <class W>
Paths shortest_paths(const EdgeList<W>& el, Vertex s) {
  const auto g = sssp::to_csr(el);
  const auto n = g.num_vertices();
  const double inf = std::numeric_limits<double>::infinity();
  Paths p{std::vector<double>(n, inf), std::vector<double>(n, 0), std::vector<std::size_t>(n, 0)};
  std::vector<bool> done(n, false);
  p.dist[s] = 0;
  p.count[s] = 1;
  for (std::size_t it = 0; it < n; ++it) {
    Vertex u = sssp::kNoVertex;
    for (Vertex v = 0; v < n; ++v) {
      if (!done[v] && p.dist[v] < inf && (u == sssp::kNoVertex || p.dist[v] < p.dist[u])) u = v;
    }
    if (u == sssp::kNoVertex) break;
    done[u] = true;
    for (const auto& a : g.arcs(u)) {
      if (done[a.to]) continue;
      const double nd = p.dist[u] + static_cast<double>(a.weight);
      if (nd < p.dist[a.to]) {
        p.dist[a.to] = nd;
        p.count[a.to] = p.count[u];
        p.hops[a.to] = p.hops[u] + 1;
      } else if (nd == p.dist[a.to]) {
        p.count[a.to] += p.count[u];
      }
    }
  }
  return p;
}

}  // namespace

TEMPLATE_TEST_CASE("small graph suite: unique names, valid graphs, valid sources", "[graph][small]",
                   std::uint64_t, double) {
  const auto suite = sssp::gen::small_graph_suite<TestType>();
  std::set<std::string> names;
  for (const auto& g : suite) {
    INFO(g.name);
    REQUIRE(names.insert(g.name).second);
    REQUIRE_NOTHROW(sssp::to_csr(g.graph));
    REQUIRE(g.source < g.graph.num_vertices);
    REQUIRE_FALSE(g.note.empty());
  }
  for (const char* required :
       {"single_vertex", "no_edges", "two_vertices", "path", "star", "cycle", "complete",
        "binary_tree", "disconnected", "multi_edges", "self_loops", "zero_weights",
        "all_equal_weights", "mixed_magnitudes", "extreme_ratio", "long_shortest_path",
        "many_equal_paths"}) {
    REQUIRE(names.count(required) == 1);
  }
}

TEMPLATE_TEST_CASE("small graph suite: structural shapes", "[graph][small]", std::uint64_t,
                   double) {
  const auto suite = sssp::gen::small_graph_suite<TestType>();
  {
    const auto& g = find(suite, "single_vertex").graph;
    REQUIRE(g.num_vertices == 1);
    REQUIRE(g.edges.empty());
  }
  {
    const auto& g = find(suite, "no_edges").graph;
    REQUIRE(g.num_vertices >= 2);
    REQUIRE(g.edges.empty());
  }
  {
    const auto& g = find(suite, "two_vertices").graph;
    REQUIRE(g.num_vertices == 2);
    REQUIRE(g.edges.size() == 1);
  }
  {
    const auto& g = find(suite, "path").graph;
    REQUIRE(g.edges.size() == g.num_vertices - 1);
    REQUIRE(components(g) == 1);
    for (auto d : degrees(g)) REQUIRE(d <= 2);
  }
  {
    const auto& g = find(suite, "star").graph;
    const auto d = degrees(g);
    REQUIRE(g.edges.size() == g.num_vertices - 1);
    REQUIRE(*std::max_element(d.begin(), d.end()) == g.num_vertices - 1);
  }
  {
    const auto& g = find(suite, "cycle").graph;
    REQUIRE(g.edges.size() == g.num_vertices);
    REQUIRE(components(g) == 1);
    for (auto d : degrees(g)) REQUIRE(d == 2);
  }
  {
    const auto& g = find(suite, "complete").graph;
    const auto n = g.num_vertices;
    REQUIRE(g.edges.size() == n * (n - 1) / 2);
    std::set<std::pair<Vertex, Vertex>> pairs;
    for (const auto& e : g.edges) {
      REQUIRE(e.u != e.v);
      pairs.insert({std::min(e.u, e.v), std::max(e.u, e.v)});
    }
    REQUIRE(pairs.size() == g.edges.size());
  }
  {
    const auto& g = find(suite, "binary_tree").graph;
    REQUIRE(g.edges.size() == g.num_vertices - 1);
    REQUIRE(components(g) == 1);
    for (auto d : degrees(g)) REQUIRE(d <= 3);
  }
  {
    const auto& g = find(suite, "disconnected").graph;
    REQUIRE(components(g) >= 2);
    const auto& ng = find(suite, "disconnected");
    const auto p = shortest_paths(g, ng.source);
    REQUIRE(std::any_of(p.dist.begin(), p.dist.end(), [](double d) { return std::isinf(d); }));
  }
  {
    const auto& g = find(suite, "multi_edges").graph;
    std::map<std::pair<Vertex, Vertex>, std::set<TestType>> w;
    for (const auto& e : g.edges) w[{std::min(e.u, e.v), std::max(e.u, e.v)}].insert(e.weight);
    REQUIRE(std::any_of(w.begin(), w.end(), [](const auto& kv) { return kv.second.size() >= 2; }));
  }
  {
    const auto& g = find(suite, "self_loops").graph;
    REQUIRE(std::any_of(g.edges.begin(), g.edges.end(), [](const auto& e) { return e.u == e.v; }));
  }
}

TEMPLATE_TEST_CASE("small graph suite: weight patterns", "[graph][small]", std::uint64_t, double) {
  const auto suite = sssp::gen::small_graph_suite<TestType>();
  {
    const auto& g = find(suite, "zero_weights").graph;
    std::size_t zeros = 0;
    for (const auto& e : g.edges) zeros += (e.weight == TestType{0});
    REQUIRE(zeros >= 1);
    REQUIRE(zeros < g.edges.size());
  }
  {
    const auto& g = find(suite, "all_equal_weights").graph;
    REQUIRE(g.edges.size() >= 3);
    for (const auto& e : g.edges) REQUIRE(e.weight == g.edges.front().weight);
  }
  auto ratio = [](const EdgeList<TestType>& g) {
    double mn = std::numeric_limits<double>::infinity(), mx = 0;
    for (const auto& e : g.edges) {
      mn = std::min(mn, static_cast<double>(e.weight));
      mx = std::max(mx, static_cast<double>(e.weight));
    }
    REQUIRE(mn > 0);
    return mx / mn;
  };
  REQUIRE(ratio(find(suite, "mixed_magnitudes").graph) >= 1e6);
  REQUIRE(ratio(find(suite, "extreme_ratio").graph) >= std::ldexp(1.0, 40));
}

TEMPLATE_TEST_CASE("small graph suite: path-structure cases", "[graph][small]", std::uint64_t,
                   double) {
  const auto suite = sssp::gen::small_graph_suite<TestType>();
  {
    const auto& ng = find(suite, "long_shortest_path");
    const auto n = ng.graph.num_vertices;
    const auto p = shortest_paths(ng.graph, ng.source);
    const auto far = std::max_element(p.hops.begin(), p.hops.end());
    REQUIRE(*far == n - 1);  // some shortest path uses every vertex
    REQUIRE(ng.graph.edges.size() > n - 1);
  }
  {
    // diamond chain of k stages: 2^k shortest paths from the source to the last vertex
    const auto g = sssp::gen::small::many_equal_paths<TestType>(10);
    REQUIRE(g.num_vertices == 31);
    const auto p = shortest_paths(g, 0);
    REQUIRE(p.count[30] == 1024.0);
    REQUIRE(p.dist[30] == 20.0);
    const auto& ng = find(suite, "many_equal_paths");
    const auto q = shortest_paths(ng.graph, ng.source);
    REQUIRE(*std::max_element(q.count.begin(), q.count.end()) >= 16.0);
  }
  {
    // extreme_ratio: the light long route must beat the single heavy edge
    const auto& ng = find(suite, "extreme_ratio");
    const auto p = shortest_paths(ng.graph, ng.source);
    const auto n = ng.graph.num_vertices;
    REQUIRE(p.hops[n - 1] > 1);
  }
}

TEST_CASE("small graph builders take sizes and weights", "[graph][small]") {
  namespace sg = sssp::gen::small;
  REQUIRE(sg::path<std::uint64_t>(5, 3).edges.size() == 4);
  REQUIRE(sg::star<double>(6, 1.5).edges.size() == 5);
  REQUIRE(sg::cycle<std::uint64_t>(4, 1).edges.size() == 4);
  REQUIRE(sg::complete<std::uint64_t>(6, 2).edges.size() == 15);
  REQUIRE(sg::binary_tree<double>(7, 1.0).edges.size() == 6);
  REQUIRE(sg::no_edges<double>(3).num_vertices == 3);
  for (const auto& e : sg::path<std::uint64_t>(5, 3).edges) REQUIRE(e.weight == 3);
}
