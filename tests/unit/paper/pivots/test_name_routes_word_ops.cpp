// Decision (b) of Q-51 applied to Algorithm 5 l.6 (arXiv 2611.34911v1, Sec 6.1): the partial
// Fisher-Yates shuffle that draws the sample S of U uses the ranges |U|, |U|-1, ..., so the bit
// count ceil(log2 h) of each draw (Sec 2, l.307-308) is read from an integer-logarithm table
// tabulated once over 1..|U| (Sec 2, l.303-304) instead of a shift loop per draw.
//  - Behaviour: the table requests exactly the same bits as the word logarithm, so every shuffle
//    (and hence every sample, estimate and certificate) is the one the per-draw word version
//    produced; checked against that version written out below as the reference.
//  - Cost: the shuffle reads the table it is given (a too-short table is rejected), and
//    name_routes tabulates once per call, |U| entries, however many trials run.
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/model/ilog.hpp"
#include "paper/model/random.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

using sssp::Edge;
using sssp::Vertex;
using sssp::paper::model::IlogTable;
using sssp::paper::model::uniform_below;
using sssp::paper::model::WordRng;
using sssp::paper::pivots::check_route_names;
using sssp::paper::pivots::name_routes;
using sssp::paper::pivots::NameRoutesParams;
using sssp::paper::pivots::partial_shuffle;
using sssp::paper::pivots::PivotRoutingInstance;

namespace {

using W = std::uint64_t;

// Alg 5 l.6 as partial_shuffle did it before the table: position i receives a uniform element
// of i..n-1, each draw computing ceil(log2 (n - i)) by the word logarithm.
template <class T>
void reference_partial_shuffle(std::vector<T>& items, std::size_t k, WordRng& rng) {
  for (std::size_t i = 0; i < k; ++i) {
    const auto j = i + static_cast<std::size_t>(uniform_below(rng, items.size() - i));
    std::swap(items[i], items[j]);
  }
}

std::vector<std::uint32_t> iota(std::size_t n) {
  std::vector<std::uint32_t> a(n);
  for (std::size_t i = 0; i < n; ++i) a[i] = static_cast<std::uint32_t>(i);
  return a;
}

// Sample sizes 0 <= k <= n: none, one, a third, half, all but one, all.
std::vector<std::size_t> sample_sizes(std::size_t n) {
  std::vector<std::size_t> ks{0, n / 3, n / 2, n};
  if (n >= 1) ks.push_back(1);
  if (n >= 1) ks.push_back(n - 1);
  return ks;
}

PivotRoutingInstance<W> complete(std::size_t t, std::uint64_t lo, std::uint64_t hi,
                                 std::uint64_t seed) {
  sssp::test::SplitMix64 r(seed);
  PivotRoutingInstance<W> g;
  g.t = t;
  for (Vertex x = 0; x < t; ++x) {
    for (Vertex y = x + 1; y < t; ++y) g.edges.push_back({x, y, r.between(lo, hi)});
  }
  return g;
}

std::uint64_t counter(const sssp::Diagnostics& d, const std::string& k) {
  const auto it = d.counters.find(k);
  return it == d.counters.end() ? 0 : it->second;
}

}  // namespace

TEST_CASE("partial_shuffle through an ilog table matches the per-draw word logarithm",
          "[name_routes][word]") {
  const std::vector<std::size_t> sizes{0, 1, 2, 3, 5, 31, 32, 33, 64, 65, 1000, 4097};
  for (std::uint64_t seed = 1; seed <= 40; ++seed) {
    for (const std::size_t n : sizes) {
      for (const std::size_t k : sample_sizes(n)) {
        INFO("seed " << seed << " n " << n << " k " << k);
        std::vector<std::uint32_t> expected = iota(n);
        WordRng ref(seed);
        reference_partial_shuffle(expected, k, ref);
        const std::uint64_t ref_next = ref.next_word();

        // The overload that tabulates its own range.
        std::vector<std::uint32_t> a = iota(n);
        WordRng rng(seed);
        partial_shuffle(std::span<std::uint32_t>(a), k, rng);
        REQUIRE(a == expected);
        REQUIRE(rng.next_word() == ref_next);  // the same number of words was consumed

        // A table over exactly 1..n and a wider one shared by several shuffles.
        for (const std::size_t range : {n, n + 4000}) {
          const IlogTable lg(range);
          std::vector<std::uint32_t> b = iota(n);
          WordRng rng2(seed);
          partial_shuffle(std::span<std::uint32_t>(b), k, rng2, lg);
          REQUIRE(b == expected);
          REQUIRE(rng2.next_word() == ref_next);
        }
      }
    }
  }
}

TEST_CASE("partial_shuffle reads its bit counts from the table it is given",
          "[name_routes][word]") {
  // The first draw has range n; a table over 1..n-1 cannot supply ceil(log2 n), so the shuffle
  // throws before it moves anything (Sec 2 l.303-304: the table covers the draw ranges).
  for (const std::size_t n : {std::size_t{2}, std::size_t{3}, std::size_t{17}, std::size_t{100}}) {
    const IlogTable short_table(n - 1);
    std::vector<std::uint32_t> a = iota(n);
    WordRng rng(7);
    REQUIRE_THROWS_AS(partial_shuffle(std::span<std::uint32_t>(a), 1, rng, short_table),
                      std::out_of_range);
    REQUIRE(a == iota(n));
  }
  // k = 0 draws nothing; ranges of size 1 (h = 1) need no logarithm.
  const IlogTable empty(0);
  std::vector<std::uint32_t> a = iota(5);
  WordRng rng(7);
  partial_shuffle(std::span<std::uint32_t>(a), 0, rng, empty);
  REQUIRE(a == iota(5));
  std::vector<std::uint32_t> one = iota(1);
  partial_shuffle(std::span<std::uint32_t>(one), 1, rng, empty);
  REQUIRE(one == iota(1));
  std::vector<std::uint32_t> b = iota(3);
  REQUIRE_THROWS_AS(partial_shuffle(std::span<std::uint32_t>(b), 4, rng, IlogTable(3)),
                    std::invalid_argument);
}

#if SSSP_DIAGNOSTICS
TEST_CASE("name_routes tabulates the sample-draw logarithms once per call", "[name_routes][word]") {
  // Mode A on t = 70, complete: |U| = 2415 > 32t, so Alg 5 l.5-14 runs. One table over 1..|U|
  // serves every trial (U only shrinks), so its build steps equal |U| whatever the trial count.
  {
    const auto g = complete(70, 1, 1000, 77);
    WordRng rng(123);
    const auto r = name_routes(g, rng);
    REQUIRE(counter(r.diag, "name_routes.trials") >= 1);
    CHECK(counter(r.diag, "name_routes.ilog_table_entries") ==
          counter(r.diag, "name_routes.representatives"));
  }
  // Mode B factors (sample t, residual 2t) on t = 30, complete (435 > 60): an accepted trial
  // removes at least half of U, so several trials run before |U| <= 2t.
  NameRoutesParams p;
  p.sample_factor = 1;
  p.residual_factor = 2;
  bool several = false;
  for (std::uint64_t seed = 1; seed <= 10; ++seed) {
    const auto g = complete(30, 1, 50, seed);
    WordRng rng(seed);
    const auto r = name_routes(g, rng, p);
    CHECK(check_route_names(g, r, sssp::test::floyd_warshall<W>(g.t, g.edges)).empty());
    REQUIRE(counter(r.diag, "name_routes.trials") >= 1);
    several = several || counter(r.diag, "name_routes.trials") >= 2;
    CHECK(counter(r.diag, "name_routes.ilog_table_entries") ==
          counter(r.diag, "name_routes.representatives"));
  }
  CHECK(several);
  // Without a trial (|U| <= 32t) nothing is tabulated.
  const auto small = complete(8, 1, 9, 3);
  WordRng rng(3);
  const auto r = name_routes(small, rng);
  CHECK(counter(r.diag, "name_routes.trials") == 0);
  CHECK(counter(r.diag, "name_routes.ilog_table_entries") == 0);
}
#else
TEST_CASE("name_routes table counter is compiled out without SSSP_DIAGNOSTICS",
          "[name_routes][word]") {
  const auto g = complete(70, 1, 1000, 77);
  WordRng rng(123);
  const auto r = name_routes(g, rng);
  CHECK(counter(r.diag, "name_routes.ilog_table_entries") == 0);
}
#endif
