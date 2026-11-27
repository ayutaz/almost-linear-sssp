// Lemmas 6.4 and 6.7 (Sec 6.2, txt l.1509 and l.1560) of arXiv 2611.34911v1, exhaustively on
// small samples [P4-PROPERTY-LEMMAS].
//
// Already covered elsewhere (not repeated here): test_estimate_property.cpp checks Lemma 6.4 on
// the quotients recorded in EstimateTrace, Z >= d_S through the library's check_estimate, the
// windows and table size, and the mean of Z / d_S over all samples and demands, on random
// std::uint64_t samples.
// New here:
//   * every simple sample on t = 2, 3, 4 pivots with weights from a multi-scale set (several
//     bands, contractions and windows), and [slow] t = 4 with five weights and t = 5 with two,
//     for std::uint64_t and quarter-integer double, every ordered demand, one seed per sample
//     (the seed varies with the sample);
//   * the quotient Q_a of Lemma 6.4 is rebuilt here from the Sec 6.2 definition (contract the
//     edges of weight < a / T^2, retain [a / T^2, T^2 a] with distinct contracted endpoints,
//     keep the incident vertices) and its distance d comes from Floyd-Warshall on it, so the
//     chain a <= beta <= d <= d_S <= d + a is checked on the paper's object, not on the trace;
//     the trace must describe the same quotient (vertex count, edge count, images);
//   * Z = 4 * 2^i a + a is formed exactly here (unsigned __int128, doubles times 4), and
//     d_S <= Z, the disconnection flags and the top-scale bound Z <= (4 * 2^l + 1) a are checked
//     without check_estimate;
//   * Lemma 6.7(3): every edge lies in at most five windows, sum |V(Q_a)| <= 10 m and the table
//     has at most sum over the processed bands of (l + 1) |V(Q_a)| entries;
//   * Lemma 6.7(2) conditional on S: for fixed samples the mean of Z / d_S over many seeds is
//     checked per demand (E[Z | S] <= 256 u^2 d_S), not pooled over samples;
//   * the constants of the proof of Lemma 6.7 for t = 1 .. 4096: T, l and H by their
//     definitions, T < 8t, l + 1 <= 7u, H_t <= 2u, sum of scales < 4 t T^2 a, sum Delta / H <= a;
//   * seeded random multigraphs (parallel edges) through the shrinking driver.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/pivots/estimate.hpp"
#include "property/paper/shrink_support.hpp"
#include "sssp/graph.hpp"
#include "support/enumerate.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

using sssp::Edge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::model::WordRng;
using sssp::paper::pivots::Demand;
using sssp::paper::pivots::EstimateTrace;
using sssp::paper::pivots::kDisconnected;
namespace pv = sssp::paper::pivots;
namespace sh = sssp::test::shrink;

namespace {

using Exact = unsigned __int128;

template <class W>
constexpr W kInf = sssp::WeightTraits<W>::infinity();

/// Exact integer image of a length: x for integers, 4 x for quarter-integer doubles.
template <class W>
Exact exact(W x) {
  if constexpr (std::is_floating_point_v<W>) {
    const W y = x * W{4};
    REQUIRE(y == std::floor(y));
    REQUIRE(y < W{0x1p100});
    return static_cast<Exact>(y);
  } else {
    return static_cast<Exact>(x);
  }
}

unsigned ceil_log2(Exact x) {
  unsigned k = 0;
  while ((Exact{1} << k) < x) ++k;
  return k;
}

struct Constants {
  unsigned tau;  // T = 2^tau, least power of two >= 4t
  unsigned ell;  // least l with 2^l >= t T^2
  unsigned eta;  // H = 2^eta, least power of two >= (16t)^5
};

Constants paper_constants(std::size_t t) {
  const Exact tt = t;
  const unsigned tau = ceil_log2(4 * tt);
  const unsigned ell = ceil_log2(tt << (2 * tau));
  const Exact s = 16 * tt;
  return {tau, ell, ceil_log2(s * s * s * s * s)};
}

std::vector<std::size_t> uf_init(std::size_t n) {
  std::vector<std::size_t> p(n);
  for (std::size_t i = 0; i < n; ++i) p[i] = i;
  return p;
}
std::size_t uf_find(std::vector<std::size_t>& p, std::size_t x) {
  while (p[x] != x) x = p[x] = p[p[x]];
  return x;
}

/// The Sec 6.2 quotient Q_a, rebuilt from the definition.
struct Quotient {
  std::vector<std::size_t> comp;     // per pivot: contracted component (union-find root)
  std::map<std::size_t, Vertex> id;  // component root -> quotient vertex, incident ones only
  std::vector<Edge<Exact>> edges;    // retained edges on quotient vertices (exact lengths)
};

template <class W>
Quotient quotient(std::size_t t, const std::vector<Edge<W>>& e, Exact a, unsigned tau) {
  Quotient q;
  auto p = uf_init(t);
  for (const auto& x : e) {
    if ((exact(x.weight) << (2 * tau)) < a) p[uf_find(p, x.u)] = uf_find(p, x.v);
  }
  q.comp.resize(t);
  for (std::size_t v = 0; v < t; ++v) q.comp[v] = uf_find(p, v);
  for (const auto& x : e) {
    const Exact w = exact(x.weight);
    if ((w << (2 * tau)) < a || w > (a << (2 * tau))) continue;
    const std::size_t cu = q.comp[x.u];
    const std::size_t cv = q.comp[x.v];
    if (cu == cv) continue;
    const auto iu = q.id.try_emplace(cu, static_cast<Vertex>(q.id.size())).first->second;
    const auto iv = q.id.try_emplace(cv, static_cast<Vertex>(q.id.size())).first->second;
    q.edges.push_back({iu, iv, w});
  }
  return q;
}

std::vector<std::vector<Exact>> exact_fw(std::size_t n, const std::vector<Edge<Exact>>& e) {
  const Exact inf = ~Exact{0};
  std::vector<std::vector<Exact>> d(n, std::vector<Exact>(n, inf));
  for (std::size_t i = 0; i < n; ++i) d[i][i] = 0;
  for (const auto& x : e) {
    d[x.u][x.v] = std::min(d[x.u][x.v], x.weight);
    d[x.v][x.u] = std::min(d[x.v][x.u], x.weight);
  }
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        if (d[i][k] != inf && d[k][j] != inf) d[i][j] = std::min(d[i][j], d[i][k] + d[k][j]);
      }
    }
  }
  return d;
}

struct Coverage {
  std::uint64_t samples = 0;
  std::uint64_t demands = 0;
  std::uint64_t disconnected = 0;
  std::uint64_t contracted = 0;        // demands whose quotient contracts at least one edge
  std::uint64_t multi_band = 0;        // samples with >= 2 bands
  std::uint64_t lifted_gap = 0;        // demands with d < d_S (the + a of Lemma 6.4 is used)
  std::uint64_t scale_above_zero = 0;  // demands answered above the finest scale
};

/// Checks Lemmas 6.4 and 6.7(1), 6.7(3) on one sample and seed; std::nullopt if they hold.
template <class W>
std::optional<std::string> estimate_violation(const sh::Case<W>& c, std::uint64_t seed,
                                              Coverage* cov) {
  const std::size_t t = c.n;
  if (t < 2) return std::nullopt;
  for (const auto& e : c.edges) {
    if (e.u == e.v || !(W{0} < e.weight)) return std::nullopt;  // outside the precondition
  }
  std::vector<Demand> q;
  for (Vertex x = 0; x < t; ++x) {
    for (Vertex y = 0; y < t; ++y) {
      if (x != y) q.push_back({x, y});
    }
  }
  WordRng rng(seed);
  EstimateTrace<W> tr;
  const auto r = pv::estimate<W>(t, c.edges, q, rng, {}, &tr);
  std::ostringstream err;
  const auto fail = [&](const std::string& what) {
    err << what;
    return std::optional<std::string>(err.str());
  };

  const Constants k = paper_constants(t);
  if (tr.constants.tau != k.tau || tr.constants.ell != k.ell || tr.constants.eta != k.eta) {
    return fail("constants T, l, H differ from Sec 6.2");
  }
  // Exact sample distances and minimax bottlenecks (Kruskal sweep written here).
  std::vector<Edge<Exact>> ex;
  for (const auto& e : c.edges) ex.push_back({e.u, e.v, exact(e.weight)});
  const auto dist = exact_fw(t, ex);
  const Exact inf = ~Exact{0};
  std::vector<std::vector<Exact>> beta(t, std::vector<Exact>(t, inf));
  {
    auto order = ex;
    std::stable_sort(order.begin(), order.end(),
                     [](const auto& a, const auto& b) { return a.weight < b.weight; });
    auto p = uf_init(t);
    for (const auto& e : order) {
      const std::size_t a = uf_find(p, e.u);
      const std::size_t b = uf_find(p, e.v);
      if (a == b) continue;
      for (std::size_t x = 0; x < t; ++x) {
        for (std::size_t y = 0; y < t; ++y) {
          if (uf_find(p, x) == a && uf_find(p, y) == b) beta[x][y] = beta[y][x] = e.weight;
        }
      }
      p[a] = b;
    }
  }
  // Band anchors (Sec 6.2): the least unassigned weight a starts a band of the weights < T a.
  std::vector<Exact> sorted;
  for (const auto& e : ex) sorted.push_back(e.weight);
  std::sort(sorted.begin(), sorted.end());
  std::map<Exact, Exact> anchor_of;
  std::vector<Exact> anchors;
  for (const Exact w : sorted) {
    if (anchors.empty() || w >= (anchors.back() << k.tau)) anchors.push_back(w);
    anchor_of[w] = anchors.back();
  }
  // Lemma 6.7(3): every edge is retained in at most five windows; sum |V(Q_a)| <= 10 m.
  std::map<Exact, Quotient> quot;
  std::map<Exact, std::vector<std::vector<Exact>>> qdist;
  std::size_t vertex_sum = 0;
  for (const Exact a : anchors) {
    quot.emplace(a, quotient(t, c.edges, a, k.tau));
    qdist.emplace(a, exact_fw(quot.at(a).id.size(), quot.at(a).edges));
    vertex_sum += quot.at(a).id.size();
  }
  for (const Exact w : sorted) {
    std::size_t windows = 0;
    for (const Exact a : anchors) windows += ((w << (2 * k.tau)) >= a && w <= (a << (2 * k.tau)));
    if (windows > 5) return fail("an edge lies in more than five windows");
  }
  if (vertex_sum > 10 * c.edges.size()) return fail("sum of quotient sizes above 10 m");

  if (tr.demand.size() != q.size() || r.index.size() != q.size()) return fail("sizes");
  std::set<Exact> processed;
  for (std::size_t i = 0; i < q.size(); ++i) {
    const auto [x, y] = q[i];
    err.str("");
    err << "demand " << x << "-" << y << ": ";
    if (dist[x][y] == inf) {
      if (r.index[i] != kDisconnected) return fail("disconnected pair not flagged");
      if (cov != nullptr) ++cov->disconnected;
      continue;
    }
    if (r.index[i] == kDisconnected || r.index[i] >= r.table.size()) {
      return fail("connected pair without a table index");
    }
    // Lemma 6.4 on the quotient rebuilt here.
    const Exact b = beta[x][y];
    const Exact a = anchor_of.at(b);
    const Quotient& qa = quot.at(a);
    processed.insert(a);
    const auto ix = qa.id.find(qa.comp[x]);
    const auto iy = qa.id.find(qa.comp[y]);
    if (ix == qa.id.end() || iy == qa.id.end()) return fail("an image is not present in Q_a");
    if (ix->second == iy->second) return fail("the images coincide in Q_a");
    const Exact d = qdist.at(a)[ix->second][iy->second];
    if (d == inf) return fail("the images are not connected in Q_a");
    const Exact ds = dist[x][y];
    if (!(a <= b && b <= d && d <= ds && ds <= d + a)) {
      return fail("a <= beta <= d <= d_S <= d + a fails");
    }
    // The trace describes the same quotient.
    const auto& dt = tr.demand[i];
    if (dt.band >= tr.bands.size()) return fail("trace band out of range");
    const auto& bt = tr.bands[dt.band];
    if (exact(bt.anchor) != a || exact(dt.beta) != b) return fail("trace anchor or beta differs");
    if (bt.num_quotient_vertices != qa.id.size() || bt.quotient_edges.size() != qa.edges.size()) {
      return fail("trace quotient size differs from Q_a");
    }
    for (Vertex u = 0; u < t; ++u) {
      const bool present = qa.id.count(qa.comp[u]) != 0;
      if ((bt.image[u] != kNoVertex) != present) return fail("trace image presence differs");
      for (Vertex v = 0; v < t && present; ++v) {
        if (qa.id.count(qa.comp[v]) == 0) continue;
        if ((bt.image[u] == bt.image[v]) != (qa.comp[u] == qa.comp[v])) {
          return fail("trace images do not follow the contraction");
        }
      }
    }
    // Lemma 6.7(1) with Z formed exactly.
    const auto& z = r.table[r.index[i]];
    if (z.scale > k.ell) return fail("scale above the top scale");
    const Exact za = exact(z.anchor);
    if (za != a) return fail("table anchor is not the band anchor");
    const Exact zz = ((Exact{4} << z.scale) + 1) * za;
    if (zz < ds) return fail("Z < d_S");
    if (cov != nullptr) {
      ++cov->demands;
      if (qa.id.size() < t) ++cov->contracted;
      if (d < ds) ++cov->lifted_gap;
      if (z.scale > 0) ++cov->scale_above_zero;
    }
  }
  // Lemma 6.7(3): at most (l + 1) parts per quotient vertex over the processed bands.
  std::size_t bound = 0;
  for (const Exact a : processed) bound += (k.ell + 1) * quot.at(a).id.size();
  if (r.table.size() > bound) return fail("table has more than (l + 1) sum |V(Q_a)| entries");
  if (tr.bands.size() != processed.size()) return fail("trace bands differ from processed bands");
  if (cov != nullptr) {
    ++cov->samples;
    if (anchors.size() >= 2) ++cov->multi_band;
  }
  return std::nullopt;
}

template <class W>
std::vector<W> scaled(const std::vector<std::uint64_t>& xs) {
  std::vector<W> out;
  for (const auto x : xs) {
    if constexpr (std::is_floating_point_v<W>) {
      out.push_back(static_cast<W>(x) / W{4});
    } else {
      out.push_back(static_cast<W>(x));
    }
  }
  return out;
}

template <class W>
void run_exhaustive(std::size_t t, const std::vector<W>& weights, Coverage& cov) {
  std::uint64_t graph = 0;
  sssp::test::for_each_small_graph<W>(t, weights, [&](const std::vector<Edge<W>>& edges) {
    ++graph;
    {
      const std::uint64_t seed = graph + 1000 * t;
      const sh::Case<W> c{t, edges, 1};
      if (const auto v = estimate_violation(c, seed, &cov)) {
        sh::Case<W> m = c;
        const auto res = sh::shrink_case(
            m, [seed](const sh::Case<W>& x) { return estimate_violation<W>(x, seed, nullptr); });
        FAIL("seed " << seed << ", " << sh::describe(c) << ": " << *v
                     << "\n  minimal: " << sh::describe(m) << ": " << res.message);
      }
    }
  });
}

void require_coverage(const Coverage& cov) {
  WARN("samples " << cov.samples << ", demands " << cov.demands << ", disconnected "
                  << cov.disconnected << ", contracted quotients " << cov.contracted
                  << ", multi-band samples " << cov.multi_band << ", d < d_S " << cov.lifted_gap
                  << ", answered above scale 0 " << cov.scale_above_zero);
  CHECK(cov.disconnected > 0);
  CHECK(cov.contracted > 0);
  CHECK(cov.multi_band > 0);
  CHECK(cov.lifted_gap > 0);
  CHECK(cov.scale_above_zero > 0);
}

}  // namespace

TEMPLATE_TEST_CASE("Lemmas 6.4 and 6.7 on every simple sample with t <= 4",
                   "[property][estimate][lemma64]", std::uint64_t, double) {
  using W = TestType;
  // T = 8 (t = 2) or 16 (t = 3, 4): 1, 20, 300, 6000 form three or four bands, and the window of
  // 6000 contracts 1 and 20 (< 6000 / 256).
  Coverage cov;
  for (std::size_t t = 2; t <= 4; ++t) run_exhaustive<W>(t, scaled<W>({1, 20, 300, 6000}), cov);
  require_coverage(cov);
}

TEMPLATE_TEST_CASE("Lemmas 6.4 and 6.7 on every simple sample with t = 4 and t = 5 [slow]",
                   "[property][estimate][lemma64][slow]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  run_exhaustive<W>(4, scaled<W>({1, 20, 300, 6000, 1000000}), cov);
  run_exhaustive<W>(5, scaled<W>({1, 300}), cov);
  require_coverage(cov);
}

TEMPLATE_TEST_CASE("Lemmas 6.4 and 6.7 on seeded random multigraphs with shrinking",
                   "[property][estimate][lemma64]", std::uint64_t, double) {
  using W = TestType;
  const auto weights = scaled<W>({1, 2, 3, 37, 1369, 50653, 1874161});
  const auto generate = [&](std::uint64_t seed) {
    sssp::test::SplitMix64 g(seed);
    sh::Case<W> c;
    c.n = 2 + g.below(11);
    const std::size_t m = 1 + g.below(2 * c.n);
    for (std::size_t i = 0; i < m; ++i) {
      const auto u = static_cast<Vertex>(g.below(c.n));
      auto v = static_cast<Vertex>(g.below(c.n));
      if (u == v) v = static_cast<Vertex>((v + 1) % c.n);
      c.edges.push_back({u, v, weights[g.below(weights.size())]});
    }
    return c;
  };
  Coverage cov;
  const auto out = sh::run_seeded<W>(
      1, 400, generate, [&](const sh::Case<W>& c) { return estimate_violation(c, 77, &cov); },
      weights);
  INFO(out.report());
  CHECK_FALSE(out.failing_seed.has_value());
  require_coverage(cov);
}

TEST_CASE("Lemma 6.7(2): E[Z | S] <= 256 u^2 d_S per demand on fixed samples",
          "[property][estimate][lemma67]") {
  using W = std::uint64_t;
  // Fixed samples S; the expectation is over the estimator's randomness only (conditional on
  // S), measured per demand.
  std::vector<sh::Case<W>> samples;
  samples.push_back({4, {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}}, 1});             // path
  samples.push_back({5, {{0, 1, 3}, {0, 2, 3}, {0, 3, 3}, {0, 4, 3}}, 1});  // star
  samples.push_back({6, {{0, 1, 1}, {1, 2, 20}, {2, 3, 300}, {3, 4, 6000}, {4, 5, 7}}, 1});
  samples.push_back({5, {{0, 1, 2}, {1, 2, 2}, {2, 3, 2}, {3, 4, 2}, {4, 0, 2}, {0, 2, 3}}, 1});
  {
    sssp::test::SplitMix64 g(4242);
    sh::Case<W> c{12, {}, 1};
    for (Vertex v = 1; v < c.n; ++v)
      c.edges.push_back({static_cast<Vertex>(g.below(v)), v, 1 + g.below(50)});
    for (int i = 0; i < 8; ++i) {
      const auto u = static_cast<Vertex>(g.below(c.n));
      const auto v = static_cast<Vertex>((u + 1 + g.below(c.n - 1)) % c.n);
      c.edges.push_back({u, v, 1 + g.below(50)});
    }
    samples.push_back(c);
  }
  constexpr std::uint64_t kRuns = 300;
  double worst = 0;
  for (std::size_t si = 0; si < samples.size(); ++si) {
    const auto& s = samples[si];
    const std::size_t t = s.n;
    const auto dist = sssp::test::floyd_warshall<W>(t, s.edges);
    std::vector<Demand> q;
    for (Vertex x = 0; x < t; ++x) {
      for (Vertex y = x + 1; y < t; ++y) q.push_back({x, y});
    }
    std::vector<double> sum(q.size(), 0.0);
    for (std::uint64_t run = 0; run < kRuns; ++run) {
      WordRng rng(run * 2654435761ULL + si);
      const auto r = pv::estimate<W>(t, s.edges, q, rng);
      for (std::size_t i = 0; i < q.size(); ++i) {
        REQUIRE(r.index[i] != kDisconnected);
        const auto& z = r.table.at(r.index[i]);
        const Exact zz = ((Exact{4} << z.scale) + 1) * z.anchor;
        REQUIRE(zz >= dist[q[i].x][q[i].y]);  // deterministic part, every run
        sum[i] += static_cast<double>(zz) / static_cast<double>(dist[q[i].x][q[i].y]);
      }
    }
    const unsigned u = ceil_log2(t + 2);
    const double allowed = 256.0 * u * u;
    for (std::size_t i = 0; i < q.size(); ++i) {
      const double mean = sum[i] / static_cast<double>(kRuns);
      INFO("sample " << si << " demand " << q[i].x << "-" << q[i].y << ": mean Z/d_S " << mean
                     << " vs 256 u^2 = " << allowed);
      CHECK(mean <= allowed);
      worst = std::max(worst, mean / allowed);
    }
  }
  WARN("Lemma 6.7(2): largest per-demand mean of Z / d_S as a fraction of 256 u^2: " << worst);
}

TEST_CASE("Lemma 6.7: the constants of Sec 6.2 and the inequalities of its proof",
          "[property][estimate][lemma67]") {
  for (std::size_t t = 1; t <= 4096; ++t) {
    INFO("t = " << t);
    const auto c = pv::estimate_constants(t);
    const Constants k = paper_constants(t);
    REQUIRE(c.tau == k.tau);
    REQUIRE(c.ell == k.ell);
    REQUIRE(c.eta == k.eta);
    const Exact tt = t;
    const Exact big_t = Exact{1} << k.tau;
    const unsigned u = ceil_log2(tt + 2);
    double ht = 0;
    for (std::size_t i = 1; i <= t; ++i) ht += 1.0 / static_cast<double>(i);
    CHECK(big_t < 8 * tt);                                          // T < 8t
    CHECK(k.ell + 1 <= 7 * u);                                      // l_a + 1 <= 7u
    CHECK(ht <= 2.0 * u);                                           // H_t <= 2u
    CHECK((Exact{1} << (k.ell + 1)) - 1 < 4 * tt * big_t * big_t);  // sum of scales < 4 t T^2 a
    CHECK((Exact{1} << (k.ell + 1)) <= (Exact{1} << k.eta));        // sum Delta / H <= a
    CHECK((Exact{1} << k.eta) >= 4 * tt * big_t * big_t);           // H >= 4 t T^2 (Q-42)
  }
}
