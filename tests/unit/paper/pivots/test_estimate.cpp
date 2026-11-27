// Sec 6.2, Algorithm 6 (Estimate) and Lemmas 6.4-6.7 of arXiv 2611.34911v1, on hand-computed
// examples. Expected constants and quotients are derived from the definitions in Sec 6.2
// ("Bands and quotients", "Centers, partitions, and table indices"), not from program output.
// Distances are checked against the Floyd-Warshall reference in tests/support.
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/pivots/estimate.hpp"
#include "sssp/graph.hpp"
#include "support/floyd_warshall.hpp"

using sssp::Edge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::model::WordRng;
using sssp::paper::pivots::bound_at_least;
using sssp::paper::pivots::check_estimate;
using sssp::paper::pivots::Demand;
using sssp::paper::pivots::estimate;
using sssp::paper::pivots::estimate_constants;
using sssp::paper::pivots::EstimateParams;
using sssp::paper::pivots::EstimateTrace;
using sssp::paper::pivots::kDisconnected;
using sssp::paper::pivots::pow2_multiple_at_least;
using sssp::paper::pivots::scaled_at_least;
using sssp::paper::pivots::TableBound;
using sssp::paper::pivots::within_radius;

namespace {

using W = std::uint64_t;

// Z = 4 * 2^scale * anchor + anchor (Alg 6 l.9), computed directly for small values.
W table_value(const TableBound<W>& b) {
  REQUIRE(b.scale + 2 < 40);
  REQUIRE(b.anchor < (W{1} << 20));
  return (b.anchor << (b.scale + 2)) + b.anchor;
}

std::vector<Demand> all_pairs(std::size_t t) {
  std::vector<Demand> q;
  for (Vertex x = 0; x < t; ++x) {
    for (Vertex y = x + 1; y < t; ++y) q.push_back({x, y});
  }
  return q;
}

}  // namespace

TEST_CASE("estimate constants follow Sec 6.2: T, top scale and H", "[estimate]") {
  // t = 2: T = least 2^tau >= 8 -> tau = 3; t T^2 = 128 -> l = 7; (16t)^5 = 32^5 = 2^25.
  auto c = estimate_constants(2);
  CHECK(c.tau == 3);
  CHECK(c.ell == 7);
  CHECK(c.eta == 25);
  // t = 3: 4t = 12 -> tau = 4; t T^2 = 768 -> l = 10; 48^5 = 254803968, 2^27 < it <= 2^28.
  c = estimate_constants(3);
  CHECK(c.tau == 4);
  CHECK(c.ell == 10);
  CHECK(c.eta == 28);
  // t = 100: 4t = 400 -> tau = 9; t T^2 = 100 * 2^18 -> l = 7 + 18 = 25;
  // 1600^5 = 1.048576e16 lies in (2^53, 2^54].
  c = estimate_constants(100);
  CHECK(c.tau == 9);
  CHECK(c.ell == 25);
  CHECK(c.eta == 54);
  // t = 2^20: 4t = 2^22; t T^2 = 2^64 -> l = 64; (16t)^5 = 2^120 exactly (multiword H).
  c = estimate_constants(std::size_t{1} << 20);
  CHECK(c.tau == 22);
  CHECK(c.ell == 64);
  CHECK(c.eta == 120);
  // Mode B (D-9) may only replace H; T and the top scale stay the paper's.
  EstimateParams p;
  p.h_log2 = 5;
  c = estimate_constants(3, p);
  CHECK(c.tau == 4);
  CHECK(c.ell == 10);
  CHECK(c.eta == 5);
}

TEST_CASE("exact scaled comparisons used by Estimate and NameRoutes", "[estimate]") {
  // Z = 4 * 2 * 3 + 3 = 27 for anchor 3 at scale 1.
  const TableBound<W> z{3, 1};
  CHECK(bound_at_least(z, W{27}));
  CHECK_FALSE(bound_at_least(z, W{28}));
  // P a* >= Z with P = 2^k.
  CHECK(scaled_at_least(W{27}, 0, z));
  CHECK_FALSE(scaled_at_least(W{26}, 0, z));
  CHECK(scaled_at_least(W{7}, 2, z));  // 28 >= 27
  CHECK_FALSE(scaled_at_least(W{6}, 2, z));
  // lambda 2^k >= d.
  CHECK(pow2_multiple_at_least(W{5}, 3, W{40}));
  CHECK_FALSE(pow2_multiple_at_least(W{5}, 3, W{41}));
  // Values beyond 64 bits are compared exactly, never wrapped.
  const W big = std::numeric_limits<W>::max() - 1;
  const TableBound<W> huge{W{1} << 62, 10};  // about 2^74
  CHECK(bound_at_least(huge, big));
  CHECK(pow2_multiple_at_least(W{1} << 63, 5, big));
  CHECK_FALSE(pow2_multiple_at_least(W{1}, 5, big));
  CHECK(scaled_at_least(W{1} << 62, 13, huge));  // 2^75 >= 2^74 + 2^62
  CHECK_FALSE(scaled_at_least(W{1} << 62, 12, huge));
}

TEST_CASE("radius test H d <= (H + j) Delta is exact (Sec 6.2)", "[estimate]") {
  // anchor 3, scale 2: Delta = 12; H = 4, j = 1: R = 12 * 5 / 4 = 15.
  const std::vector<std::uint64_t> j1{1};
  CHECK(within_radius(W{15}, W{3}, 2, 2, j1));
  CHECK_FALSE(within_radius(W{16}, W{3}, 2, 2, j1));
  // j = 0: R = Delta.
  const std::vector<std::uint64_t> j0{0};
  CHECK(within_radius(W{12}, W{3}, 2, 2, j0));
  CHECK_FALSE(within_radius(W{13}, W{3}, 2, 2, j0));
  // H = 1 (eta = 0): j must be 0 and R = Delta.
  CHECK(within_radius(W{12}, W{3}, 2, 0, std::span<const std::uint64_t>{}));
  CHECK_FALSE(within_radius(W{13}, W{3}, 2, 0, std::span<const std::uint64_t>{}));
  // Multiword: H = 2^70, j = 2^69 = H / 2: R = 1.5 Delta; anchor 2, scale 0: R = 3.
  const std::vector<std::uint64_t> jh{0, std::uint64_t{1} << 5};
  CHECK(within_radius(W{3}, W{2}, 0, 70, jh));
  CHECK_FALSE(within_radius(W{4}, W{2}, 0, 70, jh));
  // j >= H is rejected.
  const std::vector<std::uint64_t> jbad{4};
  CHECK_THROWS_AS(within_radius(W{1}, W{1}, 0, 2, jbad), std::invalid_argument);
}

TEST_CASE("estimate rejects inputs outside Algorithm 6's contract", "[estimate]") {
  WordRng rng(1);
  const std::vector<Demand> q{{0, 1}};
  const std::vector<Edge<W>> loop{{1, 1, 3}};
  CHECK_THROWS_AS(estimate<W>(3, loop, q, rng), std::invalid_argument);
  const std::vector<Edge<W>> zero{{0, 1, 0}};
  CHECK_THROWS_AS(estimate<W>(3, zero, q, rng), std::invalid_argument);
  const std::vector<Edge<W>> inf{{0, 1, std::numeric_limits<W>::max()}};
  CHECK_THROWS_AS(estimate<W>(3, inf, q, rng), std::invalid_argument);
  const std::vector<Edge<W>> out{{0, 3, 1}};
  CHECK_THROWS_AS(estimate<W>(3, out, q, rng), std::out_of_range);
  const std::vector<Edge<W>> ok{{0, 1, 1}};
  const std::vector<Demand> same{{2, 2}};
  CHECK_THROWS_AS(estimate<W>(3, ok, same, rng), std::invalid_argument);
  const std::vector<Demand> far{{0, 7}};
  CHECK_THROWS_AS(estimate<W>(3, ok, far, rng), std::out_of_range);
}

TEST_CASE("estimate on a unit path: one band, one quotient, Z >= d_S", "[estimate]") {
  // t = 3, path 0-1-2 with unit weights: a single band with anchor 1 (1 < T = 16); nothing is
  // contracted (T^2 * 1 >= 1) and both edges lie in [1/256, 256].
  const std::vector<Edge<W>> s{{0, 1, 1}, {1, 2, 1}};
  const auto q = all_pairs(3);
  WordRng rng(7);
  EstimateTrace<W> tr;
  const auto r = estimate<W>(3, s, q, rng, {}, &tr);
  REQUIRE(r.index.size() == 3);
  REQUIRE(tr.bands.size() == 1);
  CHECK(tr.bands[0].anchor == 1);
  CHECK(tr.bands[0].num_quotient_vertices == 3);
  CHECK(tr.bands[0].quotient_edges.size() == 2);
  REQUIRE(tr.bands[0].scales.size() == 11);  // i = 10 .. 0
  CHECK(tr.bands[0].scales.front().scale == 10);
  CHECK(tr.bands[0].scales.back().scale == 0);
  const std::vector<W> ds{1, 2, 1};  // d_S(0,1), d_S(0,2), d_S(1,2)
  for (std::size_t k = 0; k < q.size(); ++k) {
    REQUIRE(r.index[k] != kDisconnected);
    const auto& b = r.table.at(r.index[k]);
    CHECK(b.anchor == 1);
    CHECK(b.scale <= 10);
    CHECK(table_value(b) >= ds[k]);
    CHECK(tr.demand[k].beta == 1);
    CHECK(tr.demand[k].band == 0);
  }
  // The top scale has one part per quotient component (Lemma 6.7(1), Q-41).
  const auto& top = tr.bands[0].scales.front().table_index;
  CHECK(top[0] == top[1]);
  CHECK(top[1] == top[2]);
}

TEST_CASE("estimate bands and quotients on a geometric path (Sec 6.2)", "[estimate]") {
  // t = 4: T = 16, T^2 = 256. Path 0-1 (1), 1-2 (100), 2-3 (10000).
  // Bands: anchors 1, 100 (>= 16), 10000 (>= 1600).
  // Q_1: contract 256w < 1 (none); retain 256w >= 1 and w <= 256 -> edges 1 and 100.
  // Q_100: contract 256w < 100 (none); retain w <= 25600 -> all three edges.
  // Q_10000: contract 256w < 10000 -> edge 1; retain 256w >= 10000 -> edges 100 and 10000.
  const std::vector<Edge<W>> s{{0, 1, 1}, {1, 2, 100}, {2, 3, 10000}};
  const std::vector<Demand> q{{0, 1}, {1, 2}, {2, 3}, {0, 3}, {0, 2}};
  WordRng rng(11);
  EstimateTrace<W> tr;
  const auto r = estimate<W>(4, s, q, rng, {}, &tr);
  REQUIRE(tr.bands.size() == 3);
  CHECK(tr.bands[0].anchor == 1);
  CHECK(tr.bands[1].anchor == 100);
  CHECK(tr.bands[2].anchor == 10000);
  // Q_1
  CHECK(tr.bands[0].num_quotient_vertices == 3);
  CHECK(tr.bands[0].quotient_edges.size() == 2);
  CHECK(tr.bands[0].image[3] == kNoVertex);
  // Q_100
  CHECK(tr.bands[1].num_quotient_vertices == 4);
  CHECK(tr.bands[1].quotient_edges.size() == 3);
  // Q_10000: 0 and 1 share an image.
  CHECK(tr.bands[2].num_quotient_vertices == 3);
  CHECK(tr.bands[2].quotient_edges.size() == 2);
  CHECK(tr.bands[2].image[0] == tr.bands[2].image[1]);
  CHECK(tr.bands[2].image[1] != tr.bands[2].image[2]);
  // Merge weights and bands of the demands.
  const std::vector<W> beta{1, 100, 10000, 10000, 100};
  const std::vector<std::uint32_t> band{0, 1, 2, 2, 1};
  const std::vector<W> ds{1, 100, 10000, 10101, 101};
  for (std::size_t k = 0; k < q.size(); ++k) {
    CHECK(tr.demand[k].beta == beta[k]);
    CHECK(tr.demand[k].band == band[k]);
    REQUIRE(r.index[k] != kDisconnected);
    const auto& b = r.table.at(r.index[k]);
    CHECK(b.anchor == tr.bands[band[k]].anchor);
    CHECK(bound_at_least(b, ds[k]));
  }
  CHECK(r.diag.counters.at("estimate.bands_processed") == 3);
  CHECK(r.diag.counters.at("estimate.table_entries") == r.table.size());
  CHECK(r.diag.counters.at("estimate.max_windows_per_edge") <= 5);
}

TEST_CASE("estimate flags exactly the disconnected demands", "[estimate]") {
  const std::vector<Edge<W>> s{{0, 1, 5}, {2, 3, 7}};
  const std::vector<Demand> q{{0, 1}, {0, 2}, {1, 3}, {2, 3}, {0, 4}};
  WordRng rng(3);
  const auto r = estimate<W>(5, s, q, rng);
  REQUIRE(r.index.size() == 5);
  CHECK(r.index[0] != kDisconnected);
  CHECK(r.index[1] == kDisconnected);
  CHECK(r.index[2] == kDisconnected);
  CHECK(r.index[3] != kDisconnected);
  CHECK(r.index[4] == kDisconnected);
  CHECK(bound_at_least(r.table.at(r.index[0]), W{5}));
  CHECK(bound_at_least(r.table.at(r.index[3]), W{7}));
  CHECK(r.diag.counters.at("estimate.disconnected_demands") == 3);
}

TEST_CASE("estimate with no edges or fewer than two pivots returns no bounds", "[estimate]") {
  WordRng rng(5);
  const auto r0 = estimate<W>(1, std::span<const Edge<W>>{}, std::span<const Demand>{}, rng);
  CHECK(r0.table.empty());
  CHECK(r0.index.empty());
  const std::vector<Demand> q{{0, 1}};
  const auto r1 = estimate<W>(2, std::span<const Edge<W>>{}, q, rng);
  REQUIRE(r1.index.size() == 1);
  CHECK(r1.index[0] == kDisconnected);
}

TEST_CASE("estimate with the whole graph as sample bounds every pair", "[estimate]") {
  // A cycle with a chord and mixed weights; S is the whole graph and Q all pairs.
  const std::vector<Edge<W>> s{{0, 1, 3}, {1, 2, 40}, {2, 3, 2},   {3, 4, 900},
                               {4, 5, 1}, {5, 0, 77}, {1, 4, 5000}};
  const auto q = all_pairs(6);
  const auto dist = sssp::test::floyd_warshall<W>(6, s);
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    WordRng rng(seed);
    const auto r = estimate<W>(6, s, q, rng);
    const auto violations = check_estimate<W>(6, s, q, r, dist);
    INFO("seed " << seed << (violations.empty() ? "" : ": " + violations.front()));
    CHECK(violations.empty());
  }
}

TEST_CASE("estimate is deterministic for a fixed seed", "[estimate]") {
  const std::vector<Edge<W>> s{{0, 1, 3}, {1, 2, 4}, {2, 3, 5}, {3, 0, 6}, {0, 2, 20}};
  const auto q = all_pairs(4);
  WordRng a(42);
  WordRng b(42);
  const auto ra = estimate<W>(4, s, q, a);
  const auto rb = estimate<W>(4, s, q, b);
  REQUIRE(ra.table.size() == rb.table.size());
  for (std::size_t i = 0; i < ra.table.size(); ++i) {
    CHECK(ra.table[i].anchor == rb.table[i].anchor);
    CHECK(ra.table[i].scale == rb.table[i].scale);
  }
  CHECK(ra.index == rb.index);
}

TEST_CASE("check_estimate reports a bound below the sample distance", "[estimate]") {
  const std::vector<Edge<W>> s{{0, 1, 10}, {1, 2, 10}};
  const auto q = all_pairs(3);
  const auto dist = sssp::test::floyd_warshall<W>(3, s);
  WordRng rng(9);
  auto r = estimate<W>(3, s, q, rng);
  REQUIRE(check_estimate<W>(3, s, q, r, dist).empty());
  // Corrupt: an entry of value 4 * 1 + 1 = 5 < d_S(0, 2) = 20.
  r.table.push_back({1, 0});
  r.index[1] = static_cast<std::uint32_t>(r.table.size() - 1);
  CHECK_FALSE(check_estimate<W>(3, s, q, r, dist).empty());
  // Corrupt: a connected demand flagged as disconnected.
  auto r2 = estimate<W>(3, s, q, rng);
  r2.index[0] = kDisconnected;
  CHECK_FALSE(check_estimate<W>(3, s, q, r2, dist).empty());
}

TEST_CASE("estimate with Mode B H keeps Z >= d_S (D-9 changes only E[Z])", "[estimate]") {
  const std::vector<Edge<W>> s{{0, 1, 3}, {1, 2, 4}, {2, 3, 5}, {3, 4, 6}, {0, 4, 30}};
  const auto q = all_pairs(5);
  const auto dist = sssp::test::floyd_warshall<W>(5, s);
  EstimateParams p;
  p.h_log2 = 0;  // H = 1: R = Delta exactly
  for (std::uint64_t seed = 1; seed <= 10; ++seed) {
    WordRng rng(seed);
    EstimateTrace<W> tr;
    const auto r = estimate<W>(5, s, q, rng, p, &tr);
    CHECK(tr.constants.eta == 0);
    CHECK(check_estimate<W>(5, s, q, r, dist).empty());
  }
}

TEST_CASE("estimate on dyadic double weights bounds every pair", "[estimate]") {
  using D = double;
  const std::vector<Edge<D>> s{{0, 1, 0.5}, {1, 2, 1.25}, {2, 3, 64.0}, {0, 3, 100.0}};
  const auto q = all_pairs(4);
  const auto dist = sssp::test::floyd_warshall<D>(4, s);
  WordRng rng(17);
  const auto r = estimate<D>(4, s, q, rng);
  CHECK(check_estimate<D>(4, s, q, r, dist).empty());
}
