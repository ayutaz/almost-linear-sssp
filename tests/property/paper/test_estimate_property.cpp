// Sec 6.2, Algorithm 6 (Estimate), Lemmas 6.4-6.7 of arXiv 2611.34911v1 on random integer
// samples over many seeds. Every invariant is checked against references computed in this file
// (Floyd-Warshall on the sample and on each quotient, minimax bottlenecks, band anchors by the
// Sec 6.2 rule, the first covering center by brute force).
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/pivots/estimate.hpp"
#include "sssp/graph.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

using sssp::Edge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::model::WordRng;
using sssp::paper::pivots::check_estimate;
using sssp::paper::pivots::Demand;
using sssp::paper::pivots::estimate;
using sssp::paper::pivots::EstimateTrace;
using sssp::paper::pivots::kDisconnected;
using sssp::paper::pivots::within_radius;

namespace {

using W = std::uint64_t;
constexpr W kInf = WeightTraits<W>::infinity();

struct Sample {
  std::size_t t = 0;
  std::vector<Edge<W>> edges;
};

// O(t) edges on t pivots, loopless, possibly disconnected; weights spread over several decades
// so that several bands and contractions occur.
Sample random_sample(std::uint64_t seed) {
  sssp::test::SplitMix64 g(seed);
  Sample s;
  s.t = 2 + g.below(30);
  const std::size_t m = s.t + g.below(2 * s.t);
  for (std::size_t i = 0; i < m; ++i) {
    const auto u = static_cast<Vertex>(g.below(s.t));
    auto v = static_cast<Vertex>(g.below(s.t));
    if (u == v) v = static_cast<Vertex>((v + 1) % s.t);
    const unsigned decade = static_cast<unsigned>(g.below(4));
    W w = 1 + g.below(9);
    for (unsigned k = 0; k < decade; ++k) w *= 37;
    s.edges.push_back({u, v, w});
  }
  return s;
}

std::vector<Demand> all_pairs(std::size_t t) {
  std::vector<Demand> q;
  for (Vertex x = 0; x < t; ++x) {
    for (Vertex y = x + 1; y < t; ++y) q.push_back({x, y});
  }
  return q;
}

// Minimax (bottleneck) distance: least possible maximum edge weight on an x-y path.
std::vector<std::vector<W>> minimax(std::size_t t, const std::vector<Edge<W>>& e) {
  std::vector<std::vector<W>> b(t, std::vector<W>(t, kInf));
  for (std::size_t i = 0; i < t; ++i) b[i][i] = 0;
  for (const auto& x : e) {
    b[x.u][x.v] = std::min(b[x.u][x.v], x.weight);
    b[x.v][x.u] = std::min(b[x.v][x.u], x.weight);
  }
  for (std::size_t k = 0; k < t; ++k) {
    for (std::size_t i = 0; i < t; ++i) {
      for (std::size_t j = 0; j < t; ++j) {
        b[i][j] = std::min(b[i][j], std::max(b[i][k], b[k][j]));
      }
    }
  }
  return b;
}

// Sec 6.2: the least unassigned weight a starts a band containing the remaining weights < T a.
W band_anchor_of(W beta, std::vector<W> weights, unsigned tau) {
  std::sort(weights.begin(), weights.end());
  W anchor = 0;
  bool open = false;
  for (const W w : weights) {
    if (!open || w >= (anchor << tau)) {
      anchor = w;
      open = true;
    }
    if (w == beta) return anchor;
  }
  FAIL("merge weight is not a sample weight");
  return 0;
}

}  // namespace

TEST_CASE("estimate: Z >= d_S and disconnection flags on random samples (Lemma 6.7(1))",
          "[estimate][property]") {
  for (std::uint64_t seed = 1; seed <= 150; ++seed) {
    const auto s = random_sample(seed);
    const auto q = all_pairs(s.t);
    const auto dist = sssp::test::floyd_warshall<W>(s.t, s.edges);
    WordRng rng(seed * 7919);
    const auto r = estimate<W>(s.t, s.edges, q, rng);
    const auto v = check_estimate<W>(s.t, s.edges, q, r, dist);
    INFO("seed " << seed << (v.empty() ? "" : ": " + v.front()));
    REQUIRE(v.empty());
  }
}

TEST_CASE("estimate: quotients satisfy Lemma 6.4 on random samples", "[estimate][property]") {
  for (std::uint64_t seed = 1; seed <= 80; ++seed) {
    const auto s = random_sample(seed);
    const auto q = all_pairs(s.t);
    const auto dist = sssp::test::floyd_warshall<W>(s.t, s.edges);
    const auto mm = minimax(s.t, s.edges);
    std::vector<W> weights;
    for (const auto& e : s.edges) weights.push_back(e.weight);
    WordRng rng(seed);
    EstimateTrace<W> tr;
    const auto r = estimate<W>(s.t, s.edges, q, rng, {}, &tr);
    INFO("seed " << seed);
    REQUIRE(tr.demand.size() == q.size());
    std::vector<std::vector<std::vector<W>>> qdist;
    for (const auto& b : tr.bands) {
      qdist.push_back(sssp::test::floyd_warshall<W>(b.num_quotient_vertices, b.quotient_edges));
    }
    for (std::size_t k = 0; k < q.size(); ++k) {
      const auto [x, y] = q[k];
      if (dist[x][y] == kInf) {
        CHECK(r.index[k] == kDisconnected);
        continue;
      }
      const auto& d = tr.demand[k];
      REQUIRE(d.band < tr.bands.size());
      const auto& band = tr.bands[d.band];
      const W a = band.anchor;
      // beta is the minimum path bottleneck and a its band anchor (Sec 6.2, proof of Lemma 6.4).
      CHECK(d.beta == mm[x][y]);
      CHECK(a == band_anchor_of(d.beta, weights, tr.constants.tau));
      // Images are distinct, present and connected.
      REQUIRE(d.qx != kNoVertex);
      REQUIRE(d.qy != kNoVertex);
      CHECK(d.qx == band.image[x]);
      CHECK(d.qy == band.image[y]);
      CHECK(d.qx != d.qy);
      const W dq = qdist[d.band][d.qx][d.qy];
      REQUIRE(dq != kInf);
      // a <= beta <= d <= d_S <= d + a.
      CHECK(a <= d.beta);
      CHECK(d.beta <= dq);
      CHECK(dq <= dist[x][y]);
      CHECK(dist[x][y] <= dq + a);
    }
  }
}

TEST_CASE("estimate: windows, quotient sizes and table size (Lemma 6.7(3))",
          "[estimate][property]") {
  for (std::uint64_t seed = 1; seed <= 80; ++seed) {
    const auto s = random_sample(seed);
    const auto q = all_pairs(s.t);
    WordRng rng(seed + 1000);
    EstimateTrace<W> tr;
    const auto r = estimate<W>(s.t, s.edges, q, rng, {}, &tr);
    INFO("seed " << seed);
    const unsigned tau = tr.constants.tau;
    std::vector<unsigned> windows(s.edges.size(), 0);
    std::size_t vertex_sum = 0;
    for (const auto& b : tr.bands) {
      vertex_sum += b.num_quotient_vertices;
      CHECK(b.quotient_edges.size() == b.retained.size());
      for (std::size_t k = 0; k < b.retained.size(); ++k) {
        const auto& e = s.edges.at(b.retained[k]);
        ++windows[b.retained[k]];
        // Retained edges lie in [a / T^2, T^2 a] (Sec 6.2).
        CHECK((e.weight << (2 * tau)) >= b.anchor);
        CHECK(e.weight <= (b.anchor << (2 * tau)));
        // Their endpoints map to the quotient edge's endpoints.
        const auto& qe = b.quotient_edges[k];
        CHECK(((qe.u == b.image[e.u] && qe.v == b.image[e.v]) ||
               (qe.u == b.image[e.v] && qe.v == b.image[e.u])));
        CHECK(qe.u != qe.v);
        CHECK(qe.weight == e.weight);
      }
    }
    for (const unsigned c : windows) CHECK(c <= 5);
    CHECK(r.table.size() <= (tr.constants.ell + 1) * vertex_sum);
    CHECK(vertex_sum <= 10 * s.edges.size());
  }
}

TEST_CASE("estimate: centers, partition trees and table indices (Lemma 6.6, Q-41)",
          "[estimate][property]") {
  for (std::uint64_t seed = 1; seed <= 60; ++seed) {
    const auto s = random_sample(seed);
    const auto q = all_pairs(s.t);
    WordRng rng(seed * 31 + 5);
    EstimateTrace<W> tr;
    const auto r = estimate<W>(s.t, s.edges, q, rng, {}, &tr);
    INFO("seed " << seed);
    for (const auto& b : tr.bands) {
      const std::size_t nq = b.num_quotient_vertices;
      const auto qd = sssp::test::floyd_warshall<W>(nq, b.quotient_edges);
      REQUIRE(b.scales.size() == tr.constants.ell + 1);
      for (std::size_t si = 0; si < b.scales.size(); ++si) {
        const auto& sc = b.scales[si];
        CHECK(sc.scale == tr.constants.ell - si);  // coarsest first
        REQUIRE(sc.order.size() == nq);
        REQUIRE(sc.center.size() == nq);
        // order is a permutation.
        std::vector<std::uint32_t> sorted = sc.order;
        std::sort(sorted.begin(), sorted.end());
        for (std::size_t z = 0; z < nq; ++z) CHECK(sorted[z] == z);
        // center[z] = first center in order within R (Lemma 6.6, by brute force).
        for (Vertex z = 0; z < nq; ++z) {
          Vertex first = kNoVertex;
          for (const auto c : sc.order) {
            if (qd[c][z] != kInf &&
                within_radius(qd[c][z], b.anchor, sc.scale, tr.constants.eta, sc.j)) {
              first = c;
              break;
            }
          }
          CHECK(sc.center[z] == first);
        }
        // Parts: same part => same center (invariant of Alg 6); parts refine coarser parts;
        // the table entry of a part stores (anchor, scale).
        for (Vertex x = 0; x < nq; ++x) {
          const auto& e = r.table.at(sc.table_index[x]);
          CHECK(e.anchor == b.anchor);
          CHECK(e.scale == sc.scale);
          for (Vertex y = 0; y < nq; ++y) {
            if (sc.table_index[x] != sc.table_index[y]) continue;
            CHECK(sc.center[x] == sc.center[y]);
            if (si > 0) CHECK(b.scales[si - 1].table_index[x] == b.scales[si - 1].table_index[y]);
          }
        }
      }
      // Top scale: one part per quotient component (Lemma 6.7(1), Q-41).
      const auto& top = b.scales.front().table_index;
      for (Vertex x = 0; x < nq; ++x) {
        for (Vertex y = 0; y < nq; ++y) {
          CHECK((top[x] == top[y]) == (qd[x][y] != kInf));
          CHECK((b.component[x] == b.component[y]) == (qd[x][y] != kInf));
        }
      }
    }
    // Each connected demand's index is that of its lowest common part (Alg 6 l.10).
    for (std::size_t k = 0; k < q.size(); ++k) {
      if (r.index[k] == kDisconnected) continue;
      const auto& d = tr.demand[k];
      const auto& b = tr.bands.at(d.band);
      std::uint32_t lowest = kDisconnected;
      for (const auto& sc : b.scales) {
        if (sc.table_index[d.qx] == sc.table_index[d.qy]) lowest = sc.table_index[d.qx];
      }
      CHECK(r.index[k] == lowest);
    }
  }
}

TEST_CASE("estimate: measured E[Z / d_S] and scans per vertex (Lemma 6.6, 6.7(2))",
          "[estimate][property]") {
  // Lemma 6.7(2): E[Z | S] <= 256 u^2 d_S. Lemma 6.6: expected adjacency scans per quotient
  // vertex and scale <= H_t. Both are measured here over seeds, not assumed.
  double ratio_sum = 0;
  std::size_t ratio_n = 0;
  double worst_allowed = std::numeric_limits<double>::infinity();
  std::uint64_t scans = 0;
  double scan_allowed = 0;
  for (std::uint64_t seed = 1; seed <= 60; ++seed) {
    const auto s = random_sample(seed);
    const auto q = all_pairs(s.t);
    const auto dist = sssp::test::floyd_warshall<W>(s.t, s.edges);
    WordRng rng(seed ^ 0xabcdefULL);
    EstimateTrace<W> tr;
    const auto r = estimate<W>(s.t, s.edges, q, rng, {}, &tr);
    std::uint64_t u = 0;
    while ((std::uint64_t{1} << u) < s.t + 2) ++u;
    worst_allowed = std::min(worst_allowed, 256.0 * static_cast<double>(u * u));
    for (std::size_t k = 0; k < q.size(); ++k) {
      if (r.index[k] == kDisconnected) continue;
      const auto& b = r.table[r.index[k]];
      const double z =
          static_cast<double>(b.anchor) * (4.0 * static_cast<double>(1ULL << b.scale) + 1.0);
      ratio_sum += z / static_cast<double>(dist[q[k].x][q[k].y]);
      ++ratio_n;
    }
    double ht = 0;
    for (std::size_t i = 1; i <= s.t; ++i) ht += 1.0 / static_cast<double>(i);
    scans += r.diag.counters.at("estimate.adjacency_scans");
    scan_allowed += ht * static_cast<double>(r.diag.counters.at("estimate.vertex_scales"));
  }
  REQUIRE(ratio_n > 0);
  const double mean = ratio_sum / static_cast<double>(ratio_n);
  INFO("mean Z/d_S = " << mean << ", 256u^2 >= " << worst_allowed);
  CHECK(mean <= worst_allowed);
  INFO("scans = " << scans << ", H_t-weighted vertex-scales = " << scan_allowed);
  CHECK(static_cast<double>(scans) <= scan_allowed);
}
