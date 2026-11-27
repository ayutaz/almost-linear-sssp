// Mode B NameRoutes entry (src/experimental/paper_experimental/name_routes.hpp), task
// P5-EXP-NAMEROUTES. Hand examples and differential checks against the shared Mode A component.
//
// Expected constants are computed by hand from the paper (arXiv 2611.34911v1), not copied from
// a run:
//   * Alg 5 l.2-3: u = ceil(log2(t + 2)), g = 4u + 3, P = least power of two >= 8*256*g*u^2.
//     t = 3: u = 3, g = 15, 276480 -> P = 2^19.  t = 80: u = 7, g = 31, 3110912 -> P = 2^22.
//   * Sec 6.2 l.1473-1491: T = least power of two >= 4t, H = least power of two >= (16t)^5.
//     t = 3: T = 16, 48^5 = 254803968 -> H = 2^28.  t = 80: T = 512, 1280^5 ~ 3.44e15 -> 2^52.
//   * D-9 (Q-42): H = least power of two >= 4tT^2.
//     t = 3: 4*3*256 = 3072 -> 2^12.  t = 80: 4*80*2^18 = 83886080 -> 2^27.
//   * Alg 5 l.5-6: |S| = 16t, loop while |U| > 32t.
//
// "Mode A equality" compares every output field, every counter, the counted operations, the
// fallbacks and the number of random words drawn (SplitMix64 outputs are a bijection of its
// state, so the next output after a run identifies the number of words drawn).

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "experimental/paper_experimental/name_routes.hpp"
#include "experimental/paper_experimental/params.hpp"
#include "paper/model/random.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

namespace pe = sssp::paper_experimental;
namespace pv = sssp::paper::pivots;
using sssp::Vertex;
using sssp::paper::model::WordRng;

namespace {

using W = std::uint64_t;

std::uint64_t counter(const sssp::Diagnostics& d, const std::string& k) {
  const auto it = d.counters.find(k);
  return it == d.counters.end() ? 0 : it->second;
}

/// Number of words drawn from WordRng(seed) by a generator that now yields `next` next.
std::uint64_t words_drawn(std::uint64_t seed, std::uint64_t next) {
  WordRng probe(seed);
  for (std::uint64_t k = 0; k < (std::uint64_t{1} << 26); ++k) {
    if (probe.next_word() == next) return k;
  }
  FAIL("random word count not found");
  return 0;
}

template <class T>
pv::PivotRoutingInstance<T> complete(std::size_t t, std::uint64_t lo, std::uint64_t hi,
                                     std::uint64_t seed) {
  sssp::test::SplitMix64 r(seed);
  pv::PivotRoutingInstance<T> g;
  g.t = t;
  for (Vertex x = 0; x < t; ++x) {
    for (Vertex y = x + 1; y < t; ++y) g.edges.push_back({x, y, static_cast<T>(r.between(lo, hi))});
  }
  // Parallel copies, lighter or heavier than the first one (Alg 5 l.4 representatives).
  for (std::size_t i = 0; i < t; ++i) {
    const auto u = static_cast<Vertex>(r.below(t));
    const auto v = static_cast<Vertex>((u + 1 + r.below(t - 1)) % t);
    g.edges.push_back({v, u, static_cast<T>(r.between(lo, 2 * hi))});
  }
  return g;
}

bool same_ops(const sssp::OpCounters& a, const sssp::OpCounters& b) {
  return a.vertex_scans == b.vertex_scans && a.edge_scans == b.edge_scans &&
         a.relaxations == b.relaxations && a.comparisons == b.comparisons &&
         a.additions == b.additions && a.queue_pushes == b.queue_pushes &&
         a.queue_pops == b.queue_pops && a.decrease_keys == b.decrease_keys;
}

/// Runs the Mode B entry with `params` and the shared component with `reference` on the same
/// seed and requires identical results, counters, operations, fallbacks and random words.
template <class T>
void require_equal_runs(const pv::PivotRoutingInstance<T>& g, std::uint64_t seed,
                        const pe::PaperParams& params, const pv::NameRoutesParams& reference) {
  WordRng rb(seed);
  const auto b = pe::name_routes(g, rb, params);
  WordRng ra(seed);
  const auto a = pv::name_routes(g, ra, reference);
  INFO("t " << g.t << " M " << g.edges.size() << " seed " << seed);
  CHECK(b.names.dictionary == a.dictionary);
  CHECK(b.names.name == a.name);
  CHECK(b.names.p_log2 == a.p_log2);
  CHECK(b.names.diag.counters == a.diag.counters);
  CHECK(same_ops(b.names.diag.ops, a.diag.ops));
  REQUIRE(b.names.diag.fallbacks.size() == a.diag.fallbacks.size());
  for (std::size_t i = 0; i < a.diag.fallbacks.size(); ++i) {
    CHECK(b.names.diag.fallbacks[i].reason == a.diag.fallbacks[i].reason);
    CHECK(b.names.diag.fallbacks[i].vertices == a.diag.fallbacks[i].vertices);
    CHECK(b.names.diag.fallbacks[i].edges == a.diag.fallbacks[i].edges);
  }
  const std::uint64_t nb = rb.next_word();
  const std::uint64_t na = ra.next_word();
  CHECK(nb == na);
  CHECK(words_drawn(seed, nb) == words_drawn(seed, na));
}

pe::PaperParams with_outer(const pe::OuterParams& outer) {
  return pe::PaperParams({}, outer);
}

}  // namespace

TEST_CASE("Mode B NameRoutes setting: the paper's constants (Alg 5 l.2-6, Sec 6.2 l.1491)",
          "[paper_experimental][name_routes]") {
  const pe::PaperParams paper;
  auto s = pe::name_routes_setting(paper, 3);
  CHECK(s.t == 3);
  CHECK(s.p_log2 == 19);
  CHECK(s.paper_p_log2 == 19);
  CHECK(s.h_log2 == 28);
  CHECK(s.paper_h_log2 == 28);
  CHECK(s.sample_size == 48);
  CHECK(s.paper_sample_size == 48);
  CHECK(s.residual_limit == 96);
  CHECK(s.paper_residual_limit == 96);
  CHECK_FALSE(s.retry_cap.has_value());
  CHECK(s.is_paper());

  s = pe::name_routes_setting(paper, 80);
  CHECK(s.p_log2 == 22);
  CHECK(s.paper_p_log2 == 22);
  CHECK(s.h_log2 == 52);
  CHECK(s.paper_h_log2 == 52);
  CHECK(s.sample_size == 1280);
  CHECK(s.residual_limit == 2560);
  CHECK(s.is_paper());
}

TEST_CASE("Mode B NameRoutes setting: D-9 lowers H to pow2 >= 4tT^2, paper values kept beside",
          "[paper_experimental][name_routes]") {
  pe::OuterParams outer;
  outer.h_rule = pe::HRule::kProofMinimum;
  const auto params = with_outer(outer);
  auto s = pe::name_routes_setting(params, 3);
  CHECK(s.h_log2 == 12);
  CHECK(s.paper_h_log2 == 28);
  CHECK(s.p_log2 == 19);  // P, |S| and the residual limit stay the paper's
  CHECK(s.sample_size == 48);
  CHECK(s.residual_limit == 96);
  CHECK_FALSE(s.is_paper());
  s = pe::name_routes_setting(params, 80);
  CHECK(s.h_log2 == 27);
  CHECK(s.paper_h_log2 == 52);
  CHECK_FALSE(s.is_paper());
}

TEST_CASE("Mode B NameRoutes setting: D-17 factors and P, D-10 cap",
          "[paper_experimental][name_routes]") {
  pe::OuterParams outer;
  outer.sample_factor = 1;
  outer.residual_factor = 2;
  outer.p_log2 = 3;
  outer.retry_cap = 5;
  const auto s = pe::name_routes_setting(with_outer(outer), 10);
  CHECK(s.sample_size == 10);
  CHECK(s.paper_sample_size == 160);
  CHECK(s.residual_limit == 20);
  CHECK(s.paper_residual_limit == 320);
  CHECK(s.p_log2 == 3);
  // t = 10: u = ceil(log2 12) = 4, g = 19, 2048 * 19 * 16 = 622592 -> P = 2^20.
  CHECK(s.paper_p_log2 == 20);
  // t = 10: T = 64, H = pow2 >= 160^5 = 104857600000 -> 2^37 (2^36 = 68719476736).
  CHECK(s.h_log2 == 37);
  CHECK(s.retry_cap == std::optional<std::uint64_t>{5});
  CHECK_FALSE(s.is_paper());

  pe::OuterParams cap_only;
  cap_only.retry_cap = 64;  // a cap alone is a deviation (D-10), even if it is never reached
  CHECK_FALSE(pe::name_routes_setting(with_outer(cap_only), 10).is_paper());
}

TEST_CASE("Mode B NameRoutes with the paper's parameters equals Mode A exactly",
          "[paper_experimental][name_routes]") {
  const pe::PaperParams paper;
  // t >= 66 makes |U| = t(t-1)/2 exceed 32t, so the sampling loop of Alg 5 l.5 runs.
  for (std::uint64_t seed = 1; seed <= 4; ++seed) {
    const auto g = complete<W>(66 + 3 * seed, 1, 1000, seed);
    require_equal_runs(g, seed * 7919, paper, pv::NameRoutesParams{});
    WordRng rb(seed);
    const auto b = pe::name_routes(g, rb, paper);
    CHECK(counter(b.names.diag, "name_routes.trials") >= 1);
    CHECK_FALSE(b.retry_cap_reached);
    CHECK(b.fallback_edges == 0);
    CHECK(b.setting.is_paper());
  }
  // Below the residual threshold (no trial) and the double instantiation.
  require_equal_runs(complete<W>(9, 1, 50, 3), 11, paper, pv::NameRoutesParams{});
  require_equal_runs(complete<double>(70, 1, 900, 5), 13, paper, pv::NameRoutesParams{});
}

TEST_CASE("Mode B NameRoutes passes exactly the deviating values to the shared component",
          "[paper_experimental][name_routes]") {
  // The reference parameters are written from the hand-computed constants above, not through
  // PaperParams::name_routes_params.
  pe::OuterParams d9;
  d9.h_rule = pe::HRule::kProofMinimum;
  {
    const auto g = complete<W>(80, 1, 1000, 21);
    pv::NameRoutesParams ref;
    ref.estimate.h_log2 = 27;
    require_equal_runs(g, 31, with_outer(d9), ref);
  }
  pe::OuterParams all = d9;
  all.sample_factor = 1;
  all.residual_factor = 2;
  all.p_log2 = 3;
  all.retry_cap = 4;
  {
    const auto g = complete<W>(10, 1, 40, 22);
    pv::NameRoutesParams ref;
    ref.sample_factor = 1;
    ref.residual_factor = 2;
    ref.p_log2 = 3;
    ref.retry_cap = 4;
    // t = 10: T = 64, tT^2 = 40960 -> top scale 16; D-9: 4tT^2 -> 2^18.
    ref.estimate.h_log2 = 18;
    for (std::uint64_t seed = 1; seed <= 10; ++seed) {
      require_equal_runs(g, seed, with_outer(all), ref);
    }
  }
}

TEST_CASE("Mode B NameRoutes: D-9 is observable where the paper's H exceeds a word",
          "[paper_experimental][name_routes]") {
  // While both H fit one word, j in {0..H-1} is the top bits of one word either way, so D-9
  // barely changes the output. At t = 448 the paper's H = pow2 >= 7168^5 = 16807 * 2^50 is
  // 2^65 (two words per draw of j, Sec 2 l.307-308; Q-43), D-9's H = pow2 >= 4tT^2 with
  // T = 2048, tT^2 = 7 * 2^28 is 2^33 (one word). So the D-9 run must draw fewer words.
  const std::size_t t = 448;
  sssp::test::SplitMix64 r(448);
  pv::PivotRoutingInstance<W> g;
  g.t = t;
  for (std::size_t i = 0; i < 4 * t; ++i) {
    const auto u = static_cast<Vertex>(r.below(t));
    const auto v = static_cast<Vertex>((u + 1 + r.below(t - 1)) % t);
    g.edges.push_back({u, v, r.between(1, 5000)});
  }
  pe::OuterParams outer;
  outer.h_rule = pe::HRule::kProofMinimum;
  outer.sample_factor = 1;
  outer.residual_factor = 2;
  outer.retry_cap = 64;
  const auto params = with_outer(outer);
  const auto s = pe::name_routes_setting(params, t);
  CHECK(s.h_log2 == 33);
  CHECK(s.paper_h_log2 == 65);
  pv::NameRoutesParams ref;
  ref.sample_factor = 1;
  ref.residual_factor = 2;
  ref.retry_cap = 64;
  ref.estimate.h_log2 = 33;
  require_equal_runs(g, 5, params, ref);

  pv::NameRoutesParams paper_h = ref;
  paper_h.estimate.h_log2.reset();
  WordRng rb(5);
  const auto b = pe::name_routes(g, rb, params);
  WordRng ra(5);
  const auto a = pv::name_routes(g, ra, paper_h);
  REQUIRE(counter(b.names.diag, "name_routes.trials") >= 1);
  CHECK(words_drawn(5, rb.next_word()) < words_drawn(5, ra.next_word()));
  const auto dist = sssp::test::floyd_warshall<W>(g.t, g.edges);
  CHECK(pv::check_route_names(g, b.names, dist).empty());
  CHECK(pv::check_route_names(g, a, dist).empty());
}

TEST_CASE("Mode B NameRoutes: D-10 cap on a never-accepting instance is recorded and valid",
          "[paper_experimental][name_routes]") {
  // K_10 with unit weights. Every table value is Z = 4 * 2^i * a + a >= 5 (a is a sample weight,
  // Alg 6 l.9), so with P = 4 no a* with P a* >= Z and a* <= a_e = 1 exists (Alg 5 l.10): no edge
  // ever qualifies and every trial is rejected. Without D-10 the loop would not terminate.
  pv::PivotRoutingInstance<W> g;
  g.t = 10;
  for (Vertex x = 0; x < 10; ++x) {
    for (Vertex y = x + 1; y < 10; ++y) g.edges.push_back({x, y, 1});
  }
  pe::OuterParams outer;
  outer.sample_factor = 1;
  outer.residual_factor = 2;  // 45 > 20: the loop runs
  outer.p_log2 = 2;
  outer.retry_cap = 3;
  for (std::uint64_t seed = 1; seed <= 5; ++seed) {
    WordRng rng(seed);
    const auto r = pe::name_routes(g, rng, with_outer(outer));
    CHECK(r.retry_cap_reached);
    CHECK(r.fallback_edges == 45);
    REQUIRE(r.names.diag.fallbacks.size() == 1);
    CHECK(r.names.diag.fallbacks[0].vertices == 10);
    CHECK(r.names.diag.fallbacks[0].edges == 45);
    CHECK(counter(r.names.diag, "name_routes.trials") == 3);
    CHECK(counter(r.names.diag, "name_routes.rejected_trials") == 3);
    CHECK(counter(r.names.diag, "name_routes.accepted_rounds") == 0);
    CHECK(counter(r.names.diag, "name_routes.retry_cap_hits") == 1);
    CHECK(counter(r.names.diag, "name_routes.terminal_edges") == 45);
    // Alg 5 l.15: the residual certifies its own weight; one dictionary entry.
    CHECK(r.names.dictionary == std::vector<W>{1});
    CHECK(r.names.name == std::vector<std::uint32_t>(45, 0));
    const auto dist = sssp::test::floyd_warshall<W>(g.t, g.edges);
    CHECK(pv::check_route_names(g, r.names, dist).empty());
  }
}

TEST_CASE("Mode B NameRoutes: an unreached cap leaves no fallback",
          "[paper_experimental][name_routes]") {
  pe::OuterParams outer;
  outer.h_rule = pe::HRule::kProofMinimum;
  outer.retry_cap = 1000;
  const auto g = complete<W>(70, 1, 1000, 41);
  WordRng rng(41);
  const auto r = pe::name_routes(g, rng, with_outer(outer));
  CHECK(counter(r.names.diag, "name_routes.trials") >= 1);
  CHECK_FALSE(r.retry_cap_reached);
  CHECK(r.fallback_edges == 0);
  CHECK(r.names.diag.fallbacks.empty());
  CHECK(pv::check_route_names(g, r.names, sssp::test::floyd_warshall<W>(g.t, g.edges)).empty());
}

TEST_CASE("Mode B NameRoutes: instances without edges and invalid input",
          "[paper_experimental][name_routes]") {
  pe::OuterParams outer;
  outer.h_rule = pe::HRule::kProofMinimum;  // D-9 needs t; t < 2 has no edges (l.1422)
  const auto params = with_outer(outer);
  for (const std::size_t t : {std::size_t{0}, std::size_t{1}, std::size_t{5}}) {
    pv::PivotRoutingInstance<W> g;
    g.t = t;
    WordRng rng(1);
    const auto r = pe::name_routes(g, rng, params);
    CHECK(r.names.dictionary.empty());
    CHECK(r.names.name.empty());
    CHECK(r.setting.t == t);
    CHECK_FALSE(r.retry_cap_reached);
    CHECK(rng.next_word() == WordRng(1).next_word());  // nothing drawn
  }
  pv::PivotRoutingInstance<W> bad;
  bad.t = 3;
  bad.edges.push_back({0, 3, 1});
  WordRng rng(1);
  CHECK_THROWS_AS(pe::name_routes(bad, rng, params), std::out_of_range);
  bad.edges[0] = {1, 1, 1};
  CHECK_THROWS_AS(pe::name_routes(bad, rng, params), std::invalid_argument);
}
