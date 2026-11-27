#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"

using sssp::EdgeList;
using sssp::Vertex;
using sssp::gen::GenOptions;
using sssp::gen::MultiEdgePolicy;
using sssp::gen::WeightDist;
using sssp::gen::WeightSpec;

namespace {

// ---- independent helpers (do not use the library's own component counter) ----

template <class W>
std::size_t components(const EdgeList<W>& g) {
  std::vector<std::size_t> p(g.num_vertices);
  std::iota(p.begin(), p.end(), 0);
  auto find = [&](std::size_t x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  };
  std::size_t c = g.num_vertices;
  for (const auto& e : g.edges) {
    const auto a = find(e.u), b = find(e.v);
    if (a != b) {
      p[a] = b;
      --c;
    }
  }
  return c;
}

template <class W>
bool is_simple(const EdgeList<W>& g) {
  std::set<std::pair<Vertex, Vertex>> seen;
  for (const auto& e : g.edges) {
    if (e.u == e.v) return false;
    if (!seen.insert({std::min(e.u, e.v), std::max(e.u, e.v)}).second) return false;
  }
  return true;
}

/// Degree with the usual convention: a self-loop contributes 2.
template <class W>
std::vector<std::size_t> degrees(const EdgeList<W>& g) {
  std::vector<std::size_t> d(g.num_vertices, 0);
  for (const auto& e : g.edges) {
    ++d[e.u];
    ++d[e.v];
  }
  return d;
}

template <class W>
bool same_topology(const EdgeList<W>& a, const EdgeList<W>& b) {
  if (a.num_vertices != b.num_vertices || a.edges.size() != b.edges.size()) return false;
  for (std::size_t i = 0; i < a.edges.size(); ++i) {
    if (a.edges[i].u != b.edges[i].u || a.edges[i].v != b.edges[i].v) return false;
  }
  return true;
}

template <class W>
bool identical(const EdgeList<W>& a, const EdgeList<W>& b) {
  if (!same_topology(a, b)) return false;
  for (std::size_t i = 0; i < a.edges.size(); ++i) {
    if (a.edges[i].weight != b.edges[i].weight) return false;
  }
  return true;
}

/// O(n^2) array Dijkstra used only as an instrument: counts strict improvements of an
/// already-finite tentative distance (= decrease-key operations of a heap-based Dijkstra).
struct DijkstraTrace {
  std::vector<double> dist;
  std::vector<Vertex> parent;
  std::size_t decrease_keys = 0;
};

template <class W>
DijkstraTrace trace_dijkstra(const EdgeList<W>& el, Vertex s) {
  const auto g = sssp::to_csr(el);
  const auto n = g.num_vertices();
  const double inf = std::numeric_limits<double>::infinity();
  DijkstraTrace t{std::vector<double>(n, inf), std::vector<Vertex>(n, sssp::kNoVertex), 0};
  std::vector<bool> done(n, false);
  t.dist[s] = 0;
  for (std::size_t it = 0; it < n; ++it) {
    Vertex u = sssp::kNoVertex;
    for (Vertex v = 0; v < n; ++v) {
      if (!done[v] && t.dist[v] < inf && (u == sssp::kNoVertex || t.dist[v] < t.dist[u])) u = v;
    }
    if (u == sssp::kNoVertex) break;
    done[u] = true;
    for (const auto& a : g.arcs(u)) {
      const double nd = t.dist[u] + static_cast<double>(a.weight);
      if (!done[a.to] && nd < t.dist[a.to]) {
        if (t.dist[a.to] < inf) ++t.decrease_keys;
        t.dist[a.to] = nd;
        t.parent[a.to] = u;
      }
    }
  }
  return t;
}

GenOptions opts(std::uint64_t seed, WeightSpec w = {}) {
  return GenOptions{.weights = w, .seed = seed, .connected = false};
}

template <class W>
void require_weights_in(const EdgeList<W>& g, double lo, double hi) {
  for (const auto& e : g.edges) {
    REQUIRE(static_cast<double>(e.weight) >= lo);
    REQUIRE(static_cast<double>(e.weight) <= hi);
  }
}

}  // namespace

// ------------------------------------------------------------------ Erdős–Rényi

TEMPLATE_TEST_CASE("G(n,m) has n vertices, round(n*d/2) simple edges and in-range weights",
                   "[graph][generators][er]", std::uint64_t, double) {
  const auto g = sssp::gen::erdos_renyi<TestType>(1000, 4.0, opts(1, {.lo = 1, .hi = 50}));
  REQUIRE(g.num_vertices == 1000);
  REQUIRE(g.edges.size() == 2000);
  REQUIRE(is_simple(g));
  require_weights_in(g, 1, 50);
  REQUIRE_NOTHROW(sssp::to_csr(g));
}

TEST_CASE("G(n,m) is deterministic in (parameters, seed) and seed-sensitive",
          "[graph][generators][er][determinism]") {
  const auto a = sssp::gen::erdos_renyi<std::uint64_t>(500, 3.0, opts(7));
  const auto b = sssp::gen::erdos_renyi<std::uint64_t>(500, 3.0, opts(7));
  const auto c = sssp::gen::erdos_renyi<std::uint64_t>(500, 3.0, opts(8));
  REQUIRE(identical(a, b));
  REQUIRE_FALSE(same_topology(a, c));
}

TEST_CASE("topology does not depend on the weight distribution or weight type",
          "[graph][generators][determinism]") {
  const auto a = sssp::gen::erdos_renyi<std::uint64_t>(300, 4.0, opts(3, {.lo = 1, .hi = 9}));
  const auto b = sssp::gen::erdos_renyi<std::uint64_t>(
      300, 4.0, opts(3, {.dist = WeightDist::kWideRange, .wide_log2 = 40}));
  const auto c =
      sssp::gen::erdos_renyi<double>(300, 4.0, opts(3, {.dist = WeightDist::kUniformReal}));
  REQUIRE(same_topology(a, b));
  REQUIRE(a.edges.size() == c.edges.size());
  for (std::size_t i = 0; i < a.edges.size(); ++i) {
    REQUIRE(a.edges[i].u == c.edges[i].u);
    REQUIRE(a.edges[i].v == c.edges[i].v);
  }
}

TEST_CASE("dense G(n,m) up to the complete graph stays simple", "[graph][generators][er]") {
  const auto full = sssp::gen::erdos_renyi_gnm<std::uint64_t>(30, 435, opts(5));
  REQUIRE(full.edges.size() == 435);
  REQUIRE(is_simple(full));
  const auto dense = sssp::gen::erdos_renyi_gnm<std::uint64_t>(40, 700, opts(5));
  REQUIRE(dense.edges.size() == 700);
  REQUIRE(is_simple(dense));
  REQUIRE_THROWS_AS(sssp::gen::erdos_renyi_gnm<std::uint64_t>(30, 436, opts(5)),
                    std::invalid_argument);
}

TEST_CASE("connected option joins every component with a spanning tree of components",
          "[graph][generators][connectivity]") {
  auto o = opts(11);
  const auto sparse = sssp::gen::erdos_renyi<std::uint64_t>(2000, 0.5, o);
  REQUIRE(components(sparse) > 1);
  o.connected = true;
  const auto joined = sssp::gen::erdos_renyi<std::uint64_t>(2000, 0.5, o);
  REQUIRE(components(joined) == 1);
  // the original edges are kept as a prefix; exactly (components - 1) edges are appended
  REQUIRE(joined.edges.size() == sparse.edges.size() + components(sparse) - 1);
  for (std::size_t i = 0; i < sparse.edges.size(); ++i) {
    REQUIRE(joined.edges[i].u == sparse.edges[i].u);
    REQUIRE(joined.edges[i].v == sparse.edges[i].v);
  }
}

TEST_CASE("connect_components on an edgeless graph yields a tree", "[graph][generators]") {
  EdgeList<double> g{50, {}};
  sssp::gen::connect_components(g, WeightSpec{.dist = WeightDist::kConstant, .lo = 2}, 1);
  REQUIRE(g.edges.size() == 49);
  REQUIRE(components(g) == 1);
  REQUIRE(sssp::gen::count_components(g) == 1);
  for (const auto& e : g.edges) REQUIRE(e.weight == 2.0);
}

// ------------------------------------------------------------------ random regular

TEST_CASE("random regular graph with switching is simple and exactly d-regular",
          "[graph][generators][regular]") {
  for (std::size_t d : {3u, 4u, 7u}) {
    const auto g =
        sssp::gen::random_regular<std::uint64_t>(200, d, MultiEdgePolicy::kSwitchToSimple, opts(2));
    REQUIRE(g.edges.size() == 200 * d / 2);
    REQUIRE(is_simple(g));
    for (auto deg : degrees(g)) REQUIRE(deg == d);
  }
}

TEST_CASE("configuration model keeps or drops loops and multi-edges as requested",
          "[graph][generators][regular]") {
  const auto keep =
      sssp::gen::random_regular<std::uint64_t>(30, 6, MultiEdgePolicy::kKeep, opts(4));
  REQUIRE(keep.edges.size() == 90);
  for (auto deg : degrees(keep)) REQUIRE(deg == 6);  // self-loop counts 2
  const auto drop =
      sssp::gen::random_regular<std::uint64_t>(30, 6, MultiEdgePolicy::kDrop, opts(4));
  REQUIRE(is_simple(drop));
  REQUIRE(drop.edges.size() <= 90);
  for (auto deg : degrees(drop)) REQUIRE(deg <= 6);
  REQUIRE_THROWS_AS(sssp::gen::random_regular<std::uint64_t>(5, 3, MultiEdgePolicy::kKeep, opts(1)),
                    std::invalid_argument);  // n*d odd
  REQUIRE_THROWS_AS(
      sssp::gen::random_regular<std::uint64_t>(4, 4, MultiEdgePolicy::kSwitchToSimple, opts(1)),
      std::invalid_argument);  // simple d-regular needs d < n
}

// ------------------------------------------------------------------ Barabási–Albert

TEST_CASE("Barabasi-Albert graph is simple, connected, with m edges per new vertex",
          "[graph][generators][ba]") {
  const std::size_t n = 2000, m = 3;
  const auto g = sssp::gen::barabasi_albert<std::uint64_t>(n, m, opts(6));
  const std::size_t m0 = m + 1;  // seed clique
  REQUIRE(g.edges.size() == m0 * (m0 - 1) / 2 + (n - m0) * m);
  REQUIRE(is_simple(g));
  REQUIRE(components(g) == 1);
  const auto deg = degrees(g);
  REQUIRE(*std::min_element(deg.begin(), deg.end()) >= m);
  // preferential attachment produces hubs far above the mean degree (~2m)
  REQUIRE(*std::max_element(deg.begin(), deg.end()) >= 10 * m);
}

// ------------------------------------------------------------------ grid / road-like

TEMPLATE_TEST_CASE("2D grid has 4-neighbour edges only and is connected",
                   "[graph][generators][grid]", std::uint64_t, double) {
  const std::size_t r = 7, c = 11;
  const auto g = sssp::gen::grid_2d<TestType>(r, c, opts(1));
  REQUIRE(g.num_vertices == r * c);
  REQUIRE(g.edges.size() == r * (c - 1) + c * (r - 1));
  REQUIRE(is_simple(g));
  REQUIRE(components(g) == 1);
  for (const auto& e : g.edges) {
    const long du = static_cast<long>(e.u / c) - static_cast<long>(e.v / c);
    const long dv = static_cast<long>(e.u % c) - static_cast<long>(e.v % c);
    REQUIRE(std::labs(du) + std::labs(dv) == 1);
  }
}

TEMPLATE_TEST_CASE("road-like graph is connected with positive local geometric weights",
                   "[graph][generators][road]", std::uint64_t, double) {
  sssp::gen::RoadLikeOptions ro;
  const auto g = sssp::gen::road_like<TestType>(30, 40, ro, opts(9));
  REQUIRE(g.num_vertices == 1200);
  REQUIRE(components(g) == 1);
  REQUIRE(is_simple(g));
  const double scale = std::is_integral_v<TestType> ? ro.integer_scale : 1.0;
  // Euclidean length of a local edge: at most the shortcut radius (diagonal) plus jitter
  const double max_len = (ro.shortcut_radius * std::sqrt(2.0) + 2 * ro.jitter * std::sqrt(2.0));
  for (const auto& e : g.edges) {
    REQUIRE(static_cast<double>(e.weight) > 0.0);
    REQUIRE(static_cast<double>(e.weight) <= max_len * scale + 1);
    const long dr = static_cast<long>(e.u / 40) - static_cast<long>(e.v / 40);
    const long dc = static_cast<long>(e.u % 40) - static_cast<long>(e.v % 40);
    REQUIRE(std::labs(dr) <= ro.shortcut_radius);
    REQUIRE(std::labs(dc) <= ro.shortcut_radius);
  }
  // some grid edges are removed and some shortcuts are added
  std::size_t non_grid = 0;
  for (const auto& e : g.edges) {
    const long dr = static_cast<long>(e.u / 40) - static_cast<long>(e.v / 40);
    const long dc = static_cast<long>(e.u % 40) - static_cast<long>(e.v % 40);
    if (std::labs(dr) + std::labs(dc) != 1) ++non_grid;
  }
  REQUIRE(non_grid > 0);
  REQUIRE(g.edges.size() - non_grid < 30 * 39 + 40 * 29);
}

TEST_CASE("road-like weights follow the triangle inequality on grid neighbours",
          "[graph][generators][road]") {
  // Weight equals the Euclidean distance of perturbed points; jitter < 1/2 keeps every grid
  // edge length within [1 - 2*jitter, sqrt(1 + (2*jitter)^2)].
  sssp::gen::RoadLikeOptions ro{.jitter = 0.2};
  const auto g = sssp::gen::road_like<double>(10, 10, ro, opts(3));
  for (const auto& e : g.edges) {
    const long dr = static_cast<long>(e.u / 10) - static_cast<long>(e.v / 10);
    const long dc = static_cast<long>(e.u % 10) - static_cast<long>(e.v % 10);
    if (std::labs(dr) + std::labs(dc) == 1) {
      REQUIRE(e.weight >= 1.0 - 2 * 0.2 - 1e-12);
      REQUIRE(e.weight <= std::sqrt(1.0 + 0.4 * 0.4) + 0.4 + 1e-12);
    }
  }
}

// ------------------------------------------------------------------ trees

TEST_CASE("random tree has n-1 edges and is connected", "[graph][generators][tree]") {
  for (std::size_t n : {1u, 2u, 3u, 10u, 1000u}) {
    const auto t = sssp::gen::random_tree<std::uint64_t>(n, opts(n));
    REQUIRE(t.num_vertices == n);
    REQUIRE(t.edges.size() == n - 1);
    REQUIRE(components(t) == 1);
  }
}

TEST_CASE("random tree is uniform over labelled trees (n=3: three trees)",
          "[graph][generators][tree]") {
  // Cayley: 3^(3-2) = 3 labelled trees on 3 vertices, identified by their centre vertex.
  std::map<Vertex, int> centre_count;
  constexpr int kTrials = 3000;
  for (int s = 0; s < kTrials; ++s) {
    const auto t = sssp::gen::random_tree<std::uint64_t>(3, opts(static_cast<std::uint64_t>(s)));
    const auto d = degrees(t);
    const auto c = static_cast<Vertex>(std::max_element(d.begin(), d.end()) - d.begin());
    ++centre_count[c];
  }
  REQUIRE(centre_count.size() == 3);
  for (const auto& [c, k] : centre_count) {
    REQUIRE(k > 900);
    REQUIRE(k < 1100);
  }
}

// ------------------------------------------------------------------ scale-free

TEST_CASE("Chung-Lu graph is simple with a heavy-tailed degree sequence",
          "[graph][generators][chung-lu]") {
  const auto g = sssp::gen::chung_lu<std::uint64_t>(5000, 4.0, 2.2, opts(12));
  REQUIRE(g.num_vertices == 5000);
  REQUIRE(g.edges.size() == 10000);
  REQUIRE(is_simple(g));
  const auto deg = degrees(g);
  REQUIRE(*std::max_element(deg.begin(), deg.end()) >= 40);  // >= 10x the mean degree
  // low-index vertices have the largest expected weights
  std::size_t head = 0, tail = 0;
  for (std::size_t i = 0; i < 50; ++i) head += deg[i];
  for (std::size_t i = 4950; i < 5000; ++i) tail += deg[i];
  REQUIRE(head > 5 * tail);
  REQUIRE_THROWS_AS(sssp::gen::chung_lu<std::uint64_t>(100, 4.0, 2.0, opts(1)),
                    std::invalid_argument);  // exponent must be > 2
}

// ------------------------------------------------------------------ adversarial

TEST_CASE("decrease-key adversary forces Theta(m) decrease-keys in Dijkstra",
          "[graph][generators][adversarial]") {
  const std::size_t k = 40;
  const auto g = sssp::gen::adversarial_decrease_key<std::uint64_t>(k);
  REQUIRE(g.num_vertices == k + 1);
  REQUIRE(g.edges.size() == k + k * (k - 1) / 2);
  const auto t = trace_dijkstra(g, 0);
  // by construction dist(v_j) = j, and settling v_i lowers the key of every v_j, j > i
  for (std::size_t j = 1; j <= k; ++j) REQUIRE(t.dist[j] == static_cast<double>(j));
  REQUIRE(t.decrease_keys == k * (k - 1) / 2);
}

TEMPLATE_TEST_CASE("long-path adversary has a shortest path with n-1 hops",
                   "[graph][generators][adversarial]", std::uint64_t, double) {
  const std::size_t n = 200;
  const auto g = sssp::gen::adversarial_long_path<TestType>(n);
  REQUIRE(g.num_vertices == n);
  const auto t = trace_dijkstra(g, 0);
  REQUIRE(t.dist[n - 1] == static_cast<double>(n - 1));
  std::size_t hops = 0;
  for (Vertex v = static_cast<Vertex>(n - 1); v != 0; v = t.parent[v]) ++hops;
  REQUIRE(hops == n - 1);
  REQUIRE(g.edges.size() > n - 1);  // the decoy skip edges are present
}

// ------------------------------------------------------------------ registry / generic entry

TEST_CASE("generator catalog lists every type and generate() dispatches deterministically",
          "[graph][generators][registry]") {
  const auto cat = sssp::gen::generator_catalog();
  std::set<std::string> names;
  for (const auto& info : cat) names.insert(std::string(info.name));
  for (const char* t : {"er", "regular", "ba", "grid", "road", "tree", "chung_lu",
                        "adversarial_decrease_key", "adversarial_long_path"}) {
    REQUIRE(names.count(t) == 1);
  }
  for (const auto& info : cat) {
    sssp::gen::GraphSpec spec;
    spec.type = std::string(info.name);
    spec.n = 64;
    spec.options.seed = 5;
    spec.options.connected = true;
    const auto a = sssp::gen::generate<std::uint64_t>(spec);
    const auto b = sssp::gen::generate<std::uint64_t>(spec);
    INFO(info.name);
    REQUIRE(identical(a, b));
    REQUIRE(a.num_vertices > 0);
    REQUIRE(components(a) == 1);
    REQUIRE_NOTHROW(sssp::to_csr(a));
    const auto d = sssp::gen::generate<double>(spec);
    REQUIRE_NOTHROW(sssp::to_csr(d));
  }
  sssp::gen::GraphSpec bad;
  bad.type = "nope";
  REQUIRE_THROWS_AS(sssp::gen::generate<std::uint64_t>(bad), std::invalid_argument);
}
