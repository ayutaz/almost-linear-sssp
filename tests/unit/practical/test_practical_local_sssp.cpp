// Unit tests of the Mode C (practical-optimized) local-hierarchy SSSP, task P5-PRACTICAL-LOCAL:
// the option log (every chosen optimisation is named in the result), the FallbackRecords
// (reason, vertices, edges of the contracted source component) and the argument checks.
//
// The option log and the fallbacks are part of the result: none of these checks is guarded by
// SSSP_DIAGNOSTICS, so they also run in the release presets. Only path counters are guarded.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "experimental/practical/local_sssp.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::WeightTraits;
namespace pr = sssp::practical;

namespace {

/// Value logged for `key`, or "<missing>".
std::string logged(const std::vector<pr::OptionEntry>& log, const std::string& key) {
  const auto it = std::find_if(log.begin(), log.end(), [&](const auto& e) { return e.key == key; });
  return it == log.end() ? std::string("<missing>") : it->value;
}

pr::PracticalOptions always_hierarchy() {
  pr::PracticalOptions o;
  o.direct_search_max_vertices = 0;
  return o;
}

/// Path 0 - 1 - ... - (n-1) with weights 1, 2, ..., n-1.
template <class W>
CsrGraph<W> weighted_path(std::size_t n) {
  std::vector<Edge<W>> e;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    e.push_back({static_cast<Vertex>(i), static_cast<Vertex>(i + 1), static_cast<W>(i + 1)});
  }
  return CsrGraph<W>::from_edges(n, e);
}

}  // namespace

TEST_CASE("practical-local: algorithm name and mode", "[practical][local]") {
  CHECK(pr::kLocalSsspName == "practical-local");
  CHECK(pr::kLocalSsspMode == "practical-optimized");
}

TEST_CASE("practical-local: the option log names every choice", "[practical][local][options]") {
  const pr::PracticalOptions def;
  CHECK(def.assembly == pr::Assembly::kExactHeight);
  CHECK(def.separator == 0);
  CHECK(def.direct_search_max_vertices == 16);
  CHECK(def.supplemental_offers);

  const auto log = pr::option_log(def);
  // Tunable choices.
  CHECK(logged(log, "assembly") == "exact-height(D-11)");
  CHECK(logged(log, "separator") == "L_N");
  CHECK(logged(log, "direct_search_max_vertices") == "16");
  CHECK(logged(log, "supplemental_offers") == "on");
  // Fixed choices of this implementation, logged so that a timing carries them.
  CHECK(logged(log, "thresholds") == "none(N>=2)");
  CHECK(logged(log, "normalization") == "once-reused");
  CHECK(logged(log, "validation") == "skipped");
  CHECK(logged(log, "connectivity_check") == "skipped");
  CHECK(logged(log, "tree") == "lemma-8.5-restore");
  for (const auto& e : log) {
    INFO(e.key);
    CHECK(!e.key.empty());
    CHECK(!e.value.empty());
    CHECK(!e.note.empty());
  }
  // Keys are unique.
  for (std::size_t i = 0; i < log.size(); ++i) {
    for (std::size_t j = i + 1; j < log.size(); ++j) CHECK(log[i].key != log[j].key);
  }

  pr::PracticalOptions other;
  other.assembly = pr::Assembly::kLambdaGrid;
  other.separator = 3;
  other.direct_search_max_vertices = 0;
  other.supplemental_offers = false;
  const auto log2 = pr::option_log(other);
  CHECK(logged(log2, "assembly") == "lambda-grid(D-8)");
  CHECK(logged(log2, "separator") == "3");
  CHECK(logged(log2, "direct_search_max_vertices") == "0");
  CHECK(logged(log2, "supplemental_offers") == "off(D-14)");
}

TEST_CASE("practical-local: record_options writes practical.<key>=<value> counters",
          "[practical][local][options]") {
  sssp::Diagnostics d;
  pr::PracticalOptions o;
  o.separator = 5;
  pr::record_options(o, d);
  CHECK(d.counters.at("practical.assembly=exact-height(D-11)") == 1);
  CHECK(d.counters.at("practical.separator=5") == 1);
  CHECK(d.counters.at("practical.direct_search_max_vertices=16") == 1);
  CHECK(d.counters.at("practical.supplemental_offers=on") == 1);
  CHECK(d.counters.size() == pr::option_log(o).size());
  CHECK(d.fallbacks.empty());  // an option is not a fallback
  pr::record_options(o, d);    // idempotent
  CHECK(d.counters.at("practical.separator=5") == 1);
  CHECK(d.counters.size() == pr::option_log(o).size());
}

TEMPLATE_TEST_CASE("practical-local: every result carries its option log",
                   "[practical][local][options]", std::uint64_t, double) {
  using W = TestType;
  const auto g = weighted_path<W>(6);
  for (const std::uint64_t t : {std::uint64_t{0}, std::uint64_t{16}}) {
    pr::PracticalOptions o;
    o.direct_search_max_vertices = t;
    o.assembly = pr::Assembly::kLambdaGrid;
    o.separator = 2;
    const auto r = pr::local_sssp(g, 0, o);
    for (const auto& e : pr::option_log(o)) {
      INFO(e.key << "=" << e.value);
      REQUIRE(r.diag.counters.count("practical." + e.key + "=" + e.value) == 1);
    }
  }
}

TEMPLATE_TEST_CASE("practical-local: a small component is a recorded direct search",
                   "[practical][local][fallback]", std::uint64_t, double) {
  using W = TestType;
  const auto g = weighted_path<W>(6);  // N = 6, m' = 5
  CHECK(pr::direct_search_reason(16) ==
        "practical-local: N <= 16 (direct search, binary-heap Dijkstra)");

  // N <= threshold: one record with the contracted size.
  pr::PracticalOptions direct;
  direct.direct_search_max_vertices = 6;
  const auto a = pr::local_sssp(g, 0, direct);
  REQUIRE(a.diag.fallbacks.size() == 1);
  CHECK(a.diag.fallbacks[0].reason == pr::direct_search_reason(6));
  CHECK(a.diag.fallbacks[0].vertices == 6);
  CHECK(a.diag.fallbacks[0].edges == 5);

  // N > threshold: the hierarchy handles every vertex, no record.
  pr::PracticalOptions built;
  built.direct_search_max_vertices = 5;
  const auto b = pr::local_sssp(g, 0, built);
  CHECK(b.diag.fallbacks.empty());

  // Both paths give the distances 0, 1, 3, 6, 10, 15 and a valid tree.
  const std::vector<W> want{W{0}, W{1}, W{3}, W{6}, W{10}, W{15}};
  CHECK(a.dist == want);
  CHECK(b.dist == want);
  CHECK(sssp::verify_sssp(g, 0, a).ok);
  CHECK(sssp::verify_sssp(g, 0, b).ok);
#if SSSP_DIAGNOSTICS
  CHECK(a.diag.counters.at("practical_local_calls") == 1);
  CHECK(a.diag.counters.at("practical_local_direct_searches") == 1);
  CHECK(a.diag.counters.count("practical_local_hierarchies") == 0);
  CHECK(b.diag.counters.at("practical_local_hierarchies") == 1);
  CHECK(b.diag.counters.count("practical_local_direct_searches") == 0);
  CHECK(b.diag.counters.at("local_assemble_calls") == 1);
  CHECK(b.diag.counters.at("normalize_calls") == 1);
#endif
}

TEMPLATE_TEST_CASE("practical-local: fallback sizes are those of the contracted source component",
                   "[practical][local][fallback]", std::uint64_t, double) {
  using W = TestType;
  // Zero classes {0,1,2} and {3,4}; a loop at 5; {6,7} unreachable from 0.
  // Contracted source component: classes A={0,1,2}, B={3,4}, C={5}; N = 3.
  // Edges between different classes: 2-3 (4), 0-3 (9), 4-5 (2), 0-5 (7): m' = 4.
  const std::vector<Edge<W>> e{{0, 1, W{0}}, {1, 2, W{0}}, {2, 3, W{4}}, {0, 3, W{9}}, {3, 4, W{0}},
                               {4, 5, W{2}}, {5, 5, W{1}}, {0, 5, W{7}}, {6, 7, W{1}}};
  const auto g = CsrGraph<W>::from_edges(8, e);
  const W inf = WeightTraits<W>::infinity();
  const std::vector<W> want{W{0}, W{0}, W{0}, W{4}, W{4}, W{6}, inf, inf};

  pr::PracticalOptions direct;  // default threshold 16 >= 3
  const auto a = pr::local_sssp(g, 0, direct);
  REQUIRE(a.diag.fallbacks.size() == 1);
  CHECK(a.diag.fallbacks[0].reason == pr::direct_search_reason(16));
  CHECK(a.diag.fallbacks[0].vertices == 3);
  CHECK(a.diag.fallbacks[0].edges == 4);
  CHECK(a.dist == want);
  CHECK(sssp::verify_sssp(g, 0, a).ok);

  const auto b = pr::local_sssp(g, 0, always_hierarchy());
  CHECK(b.diag.fallbacks.empty());
  CHECK(b.dist == want);
  CHECK(sssp::verify_sssp(g, 0, b).ok);
  CHECK(b.parent[6] == kNoVertex);
  CHECK(b.parent[7] == kNoVertex);
}

TEMPLATE_TEST_CASE("practical-local: a one-vertex contracted component is recorded as trivial",
                   "[practical][local][fallback]", std::uint64_t, double) {
  using W = TestType;
  CHECK(pr::kTrivialComponentReason ==
        "practical-local: N = 1 (single contracted vertex, nothing to construct)");
  // n = 1; an isolated source; a source whose component is all zero edges.
  const auto single = CsrGraph<W>::from_edges(1, {});
  const auto isolated = CsrGraph<W>::from_edges(3, std::vector<Edge<W>>{{1, 2, W{5}}});
  const auto zeros = CsrGraph<W>::from_edges(
      4, std::vector<Edge<W>>{{0, 1, W{0}}, {1, 2, W{0}}, {2, 0, W{0}}, {2, 2, W{3}}});
  for (const auto* g : {&single, &isolated, &zeros}) {
    for (const std::uint64_t t : {std::uint64_t{0}, std::uint64_t{1}, std::uint64_t{16}}) {
      pr::PracticalOptions o;
      o.direct_search_max_vertices = t;
      const auto r = pr::local_sssp(*g, 0, o);
      INFO("n=" << g->num_vertices() << " threshold=" << t);
      REQUIRE(r.diag.fallbacks.size() == 1);
      CHECK(r.diag.fallbacks[0].reason == pr::kTrivialComponentReason);
      CHECK(r.diag.fallbacks[0].vertices == 1);
      CHECK(r.diag.fallbacks[0].edges == 0);
      CHECK(r.dist[0] == W{0});
      CHECK(sssp::verify_sssp(*g, 0, r).ok);
    }
  }
  const auto r = pr::local_sssp(zeros, 0, always_hierarchy());
  CHECK(r.dist == std::vector<W>{W{0}, W{0}, W{0}, WeightTraits<W>::infinity()});
}

TEST_CASE("practical-local: invalid arguments", "[practical][local]") {
  const auto g = weighted_path<std::uint64_t>(3);
  CHECK_THROWS_AS(pr::local_sssp(g, 3, pr::PracticalOptions{}), std::out_of_range);
  CHECK_THROWS_AS(pr::local_sssp(g, kNoVertex, pr::PracticalOptions{}), std::out_of_range);
  const auto empty = CsrGraph<double>::from_edges(0, {});
  CHECK_THROWS_AS(pr::local_sssp(empty, 0, pr::PracticalOptions{}), std::out_of_range);
}
