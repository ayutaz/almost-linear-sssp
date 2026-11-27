// Random geometric graph generator (P2-GRAPH-GEOMETRIC).
//
// The edge predicate is part of the documented contract: {u, v} (u < v) is an edge iff
// dx * dx + dy * dy <= r * r in IEEE double with dx = x_u - x_v, dy = y_u - y_v. The brute-force
// reference below evaluates exactly that expression over all pairs, so it must be compiled
// without FMA contraction, like the library.
#if defined(__clang__)
#pragma clang fp contract(off)
#endif

#include <algorithm>
#include <bit>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"

using sssp::EdgeList;
using sssp::Vertex;
using sssp::gen::GenOptions;
using sssp::gen::GeometricOptions;
using sssp::gen::GeometricWeights;
using sssp::gen::Point2;
using sssp::gen::WeightDist;
using sssp::gen::WeightSpec;

namespace {

using PairList = std::vector<std::pair<Vertex, Vertex>>;

/// O(n^2) reference: every pair u < v in lexicographic order that satisfies the predicate.
PairList brute_force_pairs(const std::vector<Point2>& p, double r) {
  PairList out;
  for (std::size_t u = 0; u < p.size(); ++u) {
    for (std::size_t v = u + 1; v < p.size(); ++v) {
      const double dx = p[u].x - p[v].x;
      const double dy = p[u].y - p[v].y;
      if (dx * dx + dy * dy <= r * r) {
        out.emplace_back(static_cast<Vertex>(u), static_cast<Vertex>(v));
      }
    }
  }
  return out;
}

template <class W>
PairList pairs_of(const EdgeList<W>& g) {
  PairList out;
  out.reserve(g.edges.size());
  for (const auto& e : g.edges) out.emplace_back(e.u, e.v);
  return out;
}

template <class W>
bool identical(const EdgeList<W>& a, const EdgeList<W>& b) {
  if (a.num_vertices != b.num_vertices || a.edges.size() != b.edges.size()) return false;
  for (std::size_t i = 0; i < a.edges.size(); ++i) {
    const auto& x = a.edges[i];
    const auto& y = b.edges[i];
    if (x.u != y.u || x.v != y.v) return false;
    if constexpr (std::is_floating_point_v<W>) {
      if (std::bit_cast<std::uint64_t>(x.weight) != std::bit_cast<std::uint64_t>(y.weight)) {
        return false;
      }
    } else if (x.weight != y.weight) {
      return false;
    }
  }
  return true;
}

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

double euclid(const std::vector<Point2>& p, Vertex u, Vertex v) {
  const double dx = p[u].x - p[v].x;
  const double dy = p[u].y - p[v].y;
  return std::sqrt(dx * dx + dy * dy);
}

GenOptions opts(std::uint64_t seed, WeightSpec w = {}, bool connected = false) {
  return GenOptions{.weights = w, .seed = seed, .connected = connected};
}

const GeometricOptions kFromSpec{};
const GeometricOptions kEuclid{.weights = GeometricWeights::kEuclidean};

/// Probability that two independent uniform points of the unit square are within distance r,
/// for 0 <= r <= 1 (the distance distribution of "square line picking"):
/// F(r) = pi r^2 - 8/3 r^3 + 1/2 r^4.
double pair_probability(double r) {
  return std::numbers::pi * r * r - 8.0 / 3.0 * r * r * r + 0.5 * r * r * r * r;
}

}  // namespace

TEST_CASE("geometric points are uniform draws in the unit square, deterministic per seed",
          "[graph][generators][geometric]") {
  const auto a = sssp::gen::random_geometric_points(5000, 7);
  const auto b = sssp::gen::random_geometric_points(5000, 7);
  const auto c = sssp::gen::random_geometric_points(5000, 8);
  REQUIRE(a.size() == 5000);
  REQUIRE(sssp::gen::random_geometric_points(0, 7).empty());
  std::size_t same = 0;
  double sx = 0, sy = 0;
  std::size_t lower_left = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    REQUIRE(a[i].x >= 0.0);
    REQUIRE(a[i].x < 1.0);
    REQUIRE(a[i].y >= 0.0);
    REQUIRE(a[i].y < 1.0);
    REQUIRE(std::bit_cast<std::uint64_t>(a[i].x) == std::bit_cast<std::uint64_t>(b[i].x));
    REQUIRE(std::bit_cast<std::uint64_t>(a[i].y) == std::bit_cast<std::uint64_t>(b[i].y));
    if (a[i].x == c[i].x && a[i].y == c[i].y) ++same;
    sx += a[i].x;
    sy += a[i].y;
    if (a[i].x < 0.5 && a[i].y < 0.5) ++lower_left;
  }
  REQUIRE(same == 0);
  // Mean of 5000 U[0,1) draws: sd = sqrt(1/12 / 5000) ~ 0.0041, so 0.02 is ~5 sd. The lower-left
  // quadrant count is Binomial(5000, 1/4): sd ~ 30.6, so +-150 is ~5 sd. Both axes must be
  // independent draws for the quadrant check to hold.
  REQUIRE(std::abs(sx / 5000.0 - 0.5) < 0.02);
  REQUIRE(std::abs(sy / 5000.0 - 0.5) < 0.02);
  REQUIRE(lower_left > 1250 - 150);
  REQUIRE(lower_left < 1250 + 150);
}

TEMPLATE_TEST_CASE("geometric edges equal the all-pairs brute force, in (u,v) order",
                   "[graph][generators][geometric]", std::uint64_t, double) {
  const double radii[] = {0.0, 1e-3, 0.01, 0.05,           0.1, 0.2333, 0.31,
                          0.5, 0.7,  1.0,  std::sqrt(2.0), 3.0};
  for (const std::size_t n :
       {std::size_t{2}, std::size_t{3}, std::size_t{10}, std::size_t{57}, std::size_t{400}}) {
    for (const double r : radii) {
      for (std::uint64_t seed = 1; seed <= 3; ++seed) {
        INFO("n=" << n << " r=" << r << " seed=" << seed);
        const auto pts = sssp::gen::random_geometric_points(n, seed);
        const auto g = sssp::gen::random_geometric<TestType>(n, r, kFromSpec, opts(seed));
        REQUIRE(g.num_vertices == n);
        REQUIRE(pairs_of(g) == brute_force_pairs(pts, r));
      }
    }
  }
}

TEST_CASE("geometric brute-force agreement holds at a radius that is not a cell multiple",
          "[graph][generators][geometric]") {
  // Many radii between two consecutive 1/k values exercise the cell-size choice near its
  // rounding boundary (cells must never be narrower than r).
  const std::size_t n = 1500;
  const auto pts = sssp::gen::random_geometric_points(n, 11);
  for (int k = 2; k <= 40; k += 3) {
    for (const double r : {1.0 / k, std::nextafter(1.0 / k, 1.0), std::nextafter(1.0 / k, 0.0)}) {
      INFO("k=" << k << " r=" << r);
      const auto g = sssp::gen::random_geometric<double>(n, r, kFromSpec, opts(11));
      REQUIRE(pairs_of(g) == brute_force_pairs(pts, r));
    }
  }
}

TEST_CASE("geometric boundary cases: n = 0, n = 1, r = 0, r >= sqrt(2)",
          "[graph][generators][geometric]") {
  for (const double r : {0.0, 0.5, 2.0}) {
    const auto g0 = sssp::gen::random_geometric<std::uint64_t>(0, r, kFromSpec, opts(1));
    REQUIRE(g0.num_vertices == 0);
    REQUIRE(g0.edges.empty());
    const auto g1 = sssp::gen::random_geometric<std::uint64_t>(1, r, kEuclid, opts(1));
    REQUIRE(g1.num_vertices == 1);
    REQUIRE(g1.edges.empty());
  }
  // r = 0: an edge needs two identical points (probability ~ n^2 2^-106 with 53-bit coordinates).
  REQUIRE(sssp::gen::random_geometric<double>(2000, 0.0, kFromSpec, opts(3)).edges.empty());
  // r >= sqrt(2), the diameter of the unit square: the complete graph, in lexicographic order.
  for (const double r : {std::sqrt(2.0), 1.5, 1e300}) {
    const std::size_t n = 60;
    const auto g = sssp::gen::random_geometric<double>(n, r, kFromSpec, opts(4));
    REQUIRE(g.edges.size() == n * (n - 1) / 2);
    std::size_t i = 0;
    for (std::size_t u = 0; u < n; ++u) {
      for (std::size_t v = u + 1; v < n; ++v, ++i) {
        REQUIRE(g.edges[i].u == u);
        REQUIRE(g.edges[i].v == v);
      }
    }
  }
  // With connected, n = 0 and n = 1 stay edgeless and r = 0 becomes a tree.
  REQUIRE(sssp::gen::random_geometric<double>(1, 0.0, kFromSpec, opts(1, {}, true)).edges.empty());
  const auto t = sssp::gen::random_geometric<double>(500, 0.0, kEuclid, opts(9, {}, true));
  REQUIRE(t.edges.size() == 499);
  REQUIRE(components(t) == 1);
}

TEST_CASE("geometric radius from the average degree is sqrt(d / (pi n))",
          "[graph][generators][geometric]") {
  REQUIRE(sssp::gen::geometric_radius(10000, 8.0) == std::sqrt(8.0 / (std::numbers::pi * 10000)));
  REQUIRE(sssp::gen::geometric_radius(1000, 0.0) == 0.0);
  // Fewer than two points have no pairs; the radius is defined as 0 instead of dividing by 0.
  REQUIRE(sssp::gen::geometric_radius(0, 4.0) == 0.0);
  REQUIRE(sssp::gen::geometric_radius(1, 4.0) == 0.0);
  REQUIRE_THROWS_AS(sssp::gen::geometric_radius(10, -1.0), std::invalid_argument);
  REQUIRE_THROWS_AS(sssp::gen::geometric_radius(10, std::numeric_limits<double>::quiet_NaN()),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(sssp::gen::geometric_radius(10, std::numeric_limits<double>::infinity()),
                    std::invalid_argument);
}

TEST_CASE("geometric mean degree matches (n-1) F(r) and is close to the target degree",
          "[graph][generators][geometric]") {
  // E[2m/n] = (n-1) F(r) with F the unit-square pair-distance CDF. The radius sqrt(d/(pi n))
  // ignores the boundary, so the target d is missed by the relative boundary loss
  // 1 - F(r)/(pi r^2) ~ 8r/(3 pi) (about 1% for n = 2e4, d = 8; 0.3% for n = 2e5).
  //
  // Tolerances: m is a U-statistic of the n points; Var(m) ~ E[m] + n d^2 r / 4 (the second
  // term is the boundary-induced variation of the local degree). For n = 2e4, d = 8:
  // E[m] = 8e4, sd ~ 290 = 0.36%, so 2% is > 5 sd. For n = 2e5: sd ~ 900 = 0.11%, so 0.5%
  // is > 4 sd. The 4% bound on the target covers boundary loss + 5 sd.
  // The n = 2e5 instance also exercises the cell grid at a size where an all-pairs scan
  // (2e10 pair tests) would not finish within the test timeout.
  struct Case {
    std::size_t n;
    std::uint64_t seed;
    double rel_tol;
  };
  const double d = 8.0;
  for (const Case c :
       {Case{20000, 1, 0.02}, Case{20000, 2, 0.02}, Case{20000, 3, 0.02}, Case{200000, 1, 0.005}}) {
    INFO("n=" << c.n << " seed=" << c.seed);
    const double r = sssp::gen::geometric_radius(c.n, d);
    const auto g = sssp::gen::random_geometric<std::uint64_t>(c.n, r, kFromSpec, opts(c.seed));
    const double mean = 2.0 * static_cast<double>(g.edges.size()) / static_cast<double>(c.n);
    const double expected = static_cast<double>(c.n - 1) * pair_probability(r);
    REQUIRE(std::abs(mean - expected) <= c.rel_tol * expected);
    REQUIRE(std::abs(mean - d) <= 0.04 * d);
  }
}

TEST_CASE("geometric euclidean weights are the point distances (uint64: rounded and scaled)",
          "[graph][generators][geometric]") {
  const std::size_t n = 800;
  const double r = 0.08;
  const auto pts = sssp::gen::random_geometric_points(n, 21);
  const auto gd = sssp::gen::random_geometric<double>(n, r, kEuclid, opts(21));
  REQUIRE(!gd.edges.empty());
  for (const auto& e : gd.edges) {
    const double len = euclid(pts, e.u, e.v);
    REQUIRE(std::bit_cast<std::uint64_t>(e.weight) == std::bit_cast<std::uint64_t>(len));
    REQUIRE(e.weight <= r);
  }
  for (const double scale : {1000.0, 1e6, 1.0}) {
    INFO("scale=" << scale);
    const GeometricOptions geo{.weights = GeometricWeights::kEuclidean, .integer_scale = scale};
    const auto gi = sssp::gen::random_geometric<std::uint64_t>(n, r, geo, opts(21));
    REQUIRE(pairs_of(gi) == pairs_of(gd));
    for (const auto& e : gi.edges) {
      // Nearest integer to len * scale (ties up). Zero is a legal result for very close points.
      const double scaled = euclid(pts, e.u, e.v) * scale;
      REQUIRE(static_cast<double>(e.weight) <= scaled + 0.5);
      REQUIRE(static_cast<double>(e.weight) > scaled - 0.5);
    }
  }
  // integer_scale only affects integer weights.
  const GeometricOptions big{.weights = GeometricWeights::kEuclidean, .integer_scale = 1e9};
  REQUIRE(identical(sssp::gen::random_geometric<double>(n, r, big, opts(21)), gd));
}

TEST_CASE("geometric topology does not depend on the weight choice or weight type",
          "[graph][generators][geometric]") {
  const std::size_t n = 1000;
  const double r = sssp::gen::geometric_radius(n, 6.0);
  const auto base = sssp::gen::random_geometric<std::uint64_t>(n, r, kFromSpec, opts(5));
  const auto real = sssp::gen::random_geometric<double>(
      n, r, kFromSpec, opts(5, {.dist = WeightDist::kUniformReal, .lo = 2, .hi = 3}));
  const auto eu = sssp::gen::random_geometric<double>(n, r, kEuclid, opts(5));
  REQUIRE(pairs_of(base) == pairs_of(real));
  REQUIRE(pairs_of(base) == pairs_of(eu));
  for (const auto& e : base.edges) {
    REQUIRE(e.weight >= 1);
    REQUIRE(e.weight <= 100);
  }
  for (const auto& e : real.edges) {
    REQUIRE(e.weight >= 2.0);
    REQUIRE(e.weight < 3.0);
  }
  // Weights come from the weight sub-stream in edge order, not from a single repeated draw.
  bool varied = false;
  for (const auto& e : base.edges) varied = varied || e.weight != base.edges.front().weight;
  REQUIRE(varied);
}

TEMPLATE_TEST_CASE("geometric graphs are deterministic in (parameters, seed) and seed-sensitive",
                   "[graph][generators][geometric]", std::uint64_t, double) {
  const std::size_t n = 3000;
  const double r = sssp::gen::geometric_radius(n, 5.0);
  for (const auto& geo : {kFromSpec, kEuclid}) {
    const auto a = sssp::gen::random_geometric<TestType>(n, r, geo, opts(42));
    const auto b = sssp::gen::random_geometric<TestType>(n, r, geo, opts(42));
    const auto c = sssp::gen::random_geometric<TestType>(n, r, geo, opts(43));
    REQUIRE(identical(a, b));
    REQUIRE(pairs_of(a) != pairs_of(c));
  }
}

TEMPLATE_TEST_CASE("geometric connected option appends a spanning tree of the components",
                   "[graph][generators][geometric]", std::uint64_t, double) {
  // Average degree 1.5 leaves many components (isolated points included).
  const std::size_t n = 2000;
  const double r = sssp::gen::geometric_radius(n, 1.5);
  for (const auto& geo : {kFromSpec, kEuclid}) {
    const auto plain = sssp::gen::random_geometric<TestType>(n, r, geo, opts(8));
    const auto conn = sssp::gen::random_geometric<TestType>(n, r, geo, opts(8, {}, true));
    const std::size_t comps = components(plain);
    REQUIRE(comps > 1);
    REQUIRE(components(conn) == 1);
    REQUIRE(conn.edges.size() == plain.edges.size() + comps - 1);
    for (std::size_t i = 0; i < plain.edges.size(); ++i) {
      REQUIRE(conn.edges[i].u == plain.edges[i].u);
      REQUIRE(conn.edges[i].v == plain.edges[i].v);
      REQUIRE(conn.edges[i].weight == plain.edges[i].weight);
    }
    if (geo.weights == GeometricWeights::kEuclidean) {
      // Joining edges are weighted like every other edge: by the length between their points.
      const auto pts = sssp::gen::random_geometric_points(n, 8);
      for (std::size_t i = plain.edges.size(); i < conn.edges.size(); ++i) {
        const auto& e = conn.edges[i];
        const double len = euclid(pts, e.u, e.v);
        if constexpr (std::is_floating_point_v<TestType>) {
          REQUIRE(e.weight == len);
        } else {
          REQUIRE(e.weight == static_cast<TestType>(std::floor(len * geo.integer_scale + 0.5)));
        }
      }
    }
  }
}

TEST_CASE("geometric generator rejects invalid parameters", "[graph][generators][geometric]") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  for (const double r : {-0.1, nan, inf}) {
    REQUIRE_THROWS_AS(sssp::gen::random_geometric<double>(10, r, kFromSpec, opts(1)),
                      std::invalid_argument);
  }
  for (const double s : {0.0, -1.0, nan, inf, 0x1p53}) {
    const GeometricOptions geo{.weights = GeometricWeights::kEuclidean, .integer_scale = s};
    REQUIRE_THROWS_AS(sssp::gen::random_geometric<std::uint64_t>(10, 0.5, geo, opts(1)),
                      std::invalid_argument);
  }
  // The weight spec is validated when it is used ...
  const WeightSpec bad{.dist = WeightDist::kUniformInt, .lo = 5, .hi = 1};
  REQUIRE_THROWS_AS(sssp::gen::random_geometric<double>(10, 0.5, kFromSpec, opts(1, bad)),
                    std::invalid_argument);
  // ... and ignored with euclidean weights, also when joining components.
  REQUIRE_NOTHROW(sssp::gen::random_geometric<double>(10, 0.1, kEuclid, opts(1, bad, true)));
  REQUIRE(sssp::gen::parse_geometric_weights("random") == GeometricWeights::kFromSpec);
  REQUIRE(sssp::gen::parse_geometric_weights("euclidean") == GeometricWeights::kEuclidean);
  REQUIRE_THROWS_AS(sssp::gen::parse_geometric_weights("euclid"), std::invalid_argument);
}

TEST_CASE("generate() dispatches type geometric with an explicit or a degree-derived radius",
          "[graph][generators][geometric][registry]") {
  bool listed = false;
  for (const auto& info : sssp::gen::generator_catalog())
    listed = listed || info.name == "geometric";
  REQUIRE(listed);
  sssp::gen::GraphSpec spec;
  spec.type = "geometric";
  spec.n = 1200;
  spec.avg_degree = 7.0;
  spec.options = opts(6);
  REQUIRE(identical(sssp::gen::generate<std::uint64_t>(spec),
                    sssp::gen::random_geometric<std::uint64_t>(
                        1200, sssp::gen::geometric_radius(1200, 7.0), kFromSpec, opts(6))));
  spec.radius = 0.04;
  spec.geometric = kEuclid;
  spec.geometric.integer_scale = 5000.0;
  REQUIRE(
      identical(sssp::gen::generate<std::uint64_t>(spec),
                sssp::gen::random_geometric<std::uint64_t>(1200, 0.04, spec.geometric, opts(6))));
  REQUIRE(identical(sssp::gen::generate<double>(spec),
                    sssp::gen::random_geometric<double>(1200, 0.04, spec.geometric, opts(6))));
}
