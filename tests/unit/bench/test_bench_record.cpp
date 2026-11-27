// Runner building blocks: algorithm registry, seeded source / order selection, RSS units,
// --env-json loading and record construction (docs/design/benchmark_schema.md).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstdint>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "experimental/practical/local_sssp.hpp"
#include "sssp/baseline.hpp"
#include "sssp/bench_runner.hpp"
#include "sssp/graph.hpp"
#include "sssp/paper.hpp"
#include "sssp/result.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::Vertex;
using sssp::bench::Json;
using sssp::bench::PhaseTimer;
using sssp::bench::RunnerError;

namespace {

// The 41 record keys of schema_version "1", in documentation order (benchmark_schema.md §2).
const std::vector<std::string> kRecordKeys = {"schema_version",
                                              "run_id",
                                              "timestamp",
                                              "paper_version",
                                              "git_sha",
                                              "git_dirty",
                                              "algorithm",
                                              "mode",
                                              "params",
                                              "os",
                                              "cpu",
                                              "ram_bytes",
                                              "arch",
                                              "compiler",
                                              "build_type",
                                              "build_flags",
                                              "diagnostics_enabled",
                                              "power",
                                              "load_average",
                                              "graph_type",
                                              "n",
                                              "m",
                                              "avg_degree",
                                              "weight_dist",
                                              "weight_type",
                                              "seed",
                                              "source",
                                              "reachable",
                                              "trial",
                                              "warmup",
                                              "order_seed",
                                              "isolated",
                                              "status",
                                              "error",
                                              "times_ns",
                                              "peak_rss_bytes",
                                              "rss_after_graph_bytes",
                                              "ops",
                                              "counters",
                                              "fallback",
                                              "correctness"};

const char* const kEnvText = R"({
  "paper_version": "2611.34911v1", "git_sha": "abcdef0", "git_dirty": false,
  "os": "TestOS 1.0", "cpu": {"brand": null, "logical_cores": 8, "physical_cores": null,
  "perf_cores": null, "efficiency_cores": null}, "ram_bytes": 17179869184, "arch": "arm64",
  "compiler": null, "build_type": "Release", "build_flags": "-O3", "diagnostics_enabled": %DIAG%,
  "power": null, "load_average": [1.25, 1.0, 0.5], "extra_key": 1})";

std::string env_text(bool diag) {
  std::string t = kEnvText;
  t.replace(t.find("%DIAG%"), 6, diag ? "true" : "false");
  return t;
}

std::vector<std::string> keys_of(const Json& j) {
  std::vector<std::string> k;
  for (const auto& [name, value] : j.items()) k.push_back(name);
  return k;
}

template <class W>
CsrGraph<W> small_graph() {
  // 0-1 (4), 0-2 (1), 2-1 (2), 1-3 (5), 4 isolated.
  return CsrGraph<W>::from_edges(
      5, {Edge<W>{0, 1, W{4}}, Edge<W>{0, 2, W{1}}, Edge<W>{2, 1, W{2}}, Edge<W>{1, 3, W{5}}});
}

sssp::bench::GraphInfo graph_info() {
  return {"er", 1000, 2000, 4.0, "uniform_int", "uint64", 1};
}

}  // namespace

TEST_CASE("bench registry: baselines with schema names and modes", "[bench]") {
  const std::vector<std::string> expected = {"bellman-ford", "dijkstra-binary", "dijkstra-4ary",
                                             "dijkstra-fibonacci", "dijkstra-radix"};
  std::vector<std::string> names;  // the baseline entries of the registry
  for (const auto& e : sssp::bench::algorithm_registry()) {
    if (e.mode == "baseline") names.emplace_back(e.name);
    REQUIRE(sssp::bench::schema_mode(e.name).has_value());
    CHECK(*sssp::bench::schema_mode(e.name) == e.mode);
    CHECK(e.solve_u64 != nullptr);
  }
  std::sort(names.begin(), names.end());
  auto sorted = expected;
  std::sort(sorted.begin(), sorted.end());
  CHECK(names == sorted);
  for (const auto& n : expected) {
    const auto* e = sssp::bench::find_algorithm(n);
    REQUIRE(e != nullptr);
    CHECK(e->mode == "baseline");
    CHECK_FALSE(e->takes_params);
    CHECK(e->exact_relaxation);
    CHECK((e->solve_f64 == nullptr) == (n == "dijkstra-radix"));
    CHECK(e->max_work == (n == "bellman-ford" ? sssp::bench::kBellmanFordMaxWork : 0));
  }
  CHECK(sssp::bench::find_algorithm("paper") == nullptr);
  CHECK(sssp::bench::schema_mode("paper") == "paper-faithful");
  CHECK(sssp::bench::schema_mode("paper-experimental") == "paper-experimental");
  CHECK(sssp::bench::schema_mode("paper-optimized") == "practical-optimized");
  CHECK_FALSE(sssp::bench::schema_mode("nope").has_value());
}

TEST_CASE("bench registry: practical-local is a uint64-only Mode C entry", "[bench]") {
  // P5-PRACTICAL-LOCAL-BENCH. General doubles can make the tree restoration of Lemma 8.5 throw
  // (local_sssp.hpp, "Numerics"), so no double solver is registered: exit 3, never a fallback.
  const auto* e = sssp::bench::find_algorithm("practical-local");
  REQUIRE(e != nullptr);
  CHECK(e->mode == "practical-optimized");
  CHECK(sssp::bench::schema_mode("practical-local") == "practical-optimized");
  CHECK(e->solve_u64 != nullptr);
  CHECK(e->solve_f64 == nullptr);
  CHECK_FALSE(e->takes_params);
  CHECK_FALSE(e->exact_relaxation);
  CHECK(e->max_work == 0);
  CHECK(e->max_vertices == 0);
  CHECK(sssp::bench::algorithm_registry().size() == 8);
}

TEST_CASE("bench registry: practical-local records its fallback per record", "[bench]") {
  const auto* e = sssp::bench::find_algorithm("practical-local");
  REQUIRE(e != nullptr);
  const auto g = small_graph<std::uint64_t>();
  const auto oracle = sssp::baseline::bellman_ford(g, 0);
  PhaseTimer timer;
  const auto r = e->solve_u64(g, 0, Json::object(), timer);
  CHECK(r.dist == oracle.dist);
  REQUIRE(timer.ns(PhaseTimer::Phase::kSssp).has_value());
  // Vertex 4 is isolated: the contracted source component has N = 4 <= 16 vertices, so the
  // direct search runs and is recorded with the contracted N and m' (4 vertices, 4 edges).
  REQUIRE(r.diag.fallbacks.size() == 1);
  CHECK(r.diag.fallbacks[0].reason ==
        sssp::practical::direct_search_reason(
            sssp::practical::PracticalOptions{}.direct_search_max_vertices));
  CHECK(r.diag.fallbacks[0].reason.find("<= 16") != std::string::npos);
  CHECK(r.diag.fallbacks[0].vertices == 4);
  CHECK(r.diag.fallbacks[0].edges == 4);
}

TEST_CASE("bench registry: Mode A entry under the names of sssp/paper.hpp", "[bench]") {
  using sssp::paper::PaperVariant;
  for (const PaperVariant v : sssp::paper::kPaperVariants) {
    const std::string_view name = sssp::paper::algorithm_name(v);
    INFO(name);
    const auto* e = sssp::bench::find_algorithm(name);
    REQUIRE(e != nullptr);
    CHECK(e->name == name);
    CHECK(e->mode == "paper-faithful");
    CHECK(sssp::bench::schema_mode(name) == "paper-faithful");
    CHECK(e->solve_u64 != nullptr);
    CHECK(e->solve_f64 != nullptr);
    CHECK_FALSE(e->takes_params);
    // Not a relaxation baseline: double results are compared with for_path_hops(2(n-1)).
    CHECK_FALSE(e->exact_relaxation);
    CHECK(e->max_work == 0);
    // Only the scan base case (Theta(N^2) vertex scans) carries a vertex limit.
    CHECK(e->max_vertices ==
          (v == PaperVariant::kFaithful ? sssp::bench::kScanBaseCaseMaxVertices : 0));
  }
  for (const auto& e : sssp::bench::algorithm_registry()) {
    if (e.mode == "baseline") CHECK(e.max_vertices == 0);
  }
}

TEMPLATE_TEST_CASE("bench registry: Mode A solvers return the paper's FallbackRecord", "[bench]",
                   std::uint64_t, double) {
  const auto g = small_graph<TestType>();
  const auto oracle = sssp::baseline::bellman_ford(g, 0);
  for (const sssp::paper::PaperVariant v : sssp::paper::kPaperVariants) {
    const auto* e = sssp::bench::find_algorithm(sssp::paper::algorithm_name(v));
    REQUIRE(e != nullptr);
    PhaseTimer timer;
    const auto r = e->template solver<TestType>()(g, 0, Json::object(), timer);
    CHECK(r.dist == oracle.dist);
    const auto direct = sssp::paper::paper_sssp(g, 0, sssp::paper::PaperOptions{v});
    REQUIRE(r.diag.fallbacks.size() == 1);
    REQUIRE(direct.diag.fallbacks.size() == 1);
    CHECK(r.diag.fallbacks[0].reason == direct.diag.fallbacks[0].reason);
    CHECK(r.diag.fallbacks[0].vertices == direct.diag.fallbacks[0].vertices);
    CHECK(r.diag.fallbacks[0].edges == direct.diag.fallbacks[0].edges);
    CHECK(r.diag.fallbacks[0].reason.starts_with("outer: N <= L_N^64 (paper base case"));
    // Vertex 4 is isolated: the contracted source component has 4 vertices and 4 edges.
    CHECK(r.diag.fallbacks[0].vertices == 4);
    CHECK(r.diag.fallbacks[0].edges == 4);
  }
}

TEMPLATE_TEST_CASE("bench registry: solvers agree with Bellman-Ford and time only sssp", "[bench]",
                   std::uint64_t, double) {
  const auto g = small_graph<TestType>();
  const auto oracle = sssp::baseline::bellman_ford(g, 0);
  for (const auto& e : sssp::bench::algorithm_registry()) {
    const auto solve = e.template solver<TestType>();
    if (solve == nullptr) continue;
    INFO(e.name);
    PhaseTimer timer;
    const auto r = solve(g, 0, Json::object(), timer);
    CHECK(r.dist == oracle.dist);
    REQUIRE(timer.ns(PhaseTimer::Phase::kSssp).has_value());
    CHECK(*timer.ns(PhaseTimer::Phase::kSssp) >= 1);
    CHECK_FALSE(timer.ns(PhaseTimer::Phase::kPreprocess).has_value());
    CHECK_FALSE(timer.ns(PhaseTimer::Phase::kTree).has_value());
  }
}

TEST_CASE("bench phase timer: accumulates, never reports 0, resets to unmeasured", "[bench]") {
  PhaseTimer t;
  CHECK_FALSE(t.ns(PhaseTimer::Phase::kSssp).has_value());
  const int v = t.time(PhaseTimer::Phase::kPreprocess, [] { return 7; });
  CHECK(v == 7);
  REQUIRE(t.ns(PhaseTimer::Phase::kPreprocess).has_value());
  const auto first = *t.ns(PhaseTimer::Phase::kPreprocess);
  CHECK(first >= 1);
  t.time(PhaseTimer::Phase::kPreprocess, [] { return 0; });
  CHECK(*t.ns(PhaseTimer::Phase::kPreprocess) >= first + 1);
  t.reset();
  CHECK_FALSE(t.ns(PhaseTimer::Phase::kPreprocess).has_value());
}

TEST_CASE("bench sources: deterministic, distinct, in range", "[bench]") {
  const auto a = sssp::bench::select_sources(1000, 4, 1);
  CHECK(a == sssp::bench::select_sources(1000, 4, 1));
  CHECK(a.size() == 4);
  CHECK(std::set<Vertex>(a.begin(), a.end()).size() == 4);
  for (auto v : a) CHECK(v < 1000);
  CHECK(a != sssp::bench::select_sources(1000, 4, 2));
  // The first sources do not depend on how many are requested.
  const auto b = sssp::bench::select_sources(1000, 2, 1);
  CHECK(std::equal(b.begin(), b.end(), a.begin()));
  auto all = sssp::bench::select_sources(10, 10, 7);
  std::sort(all.begin(), all.end());
  CHECK(all == std::vector<Vertex>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
  CHECK(sssp::bench::select_sources(1, 1, 99) == std::vector<Vertex>{0});
  CHECK_THROWS_AS(sssp::bench::select_sources(3, 4, 1), std::invalid_argument);
}

TEST_CASE("bench order: per-round seeded permutation", "[bench]") {
  const auto o = sssp::bench::round_order(5, 12345, 0, 0);
  auto sorted = o;
  std::sort(sorted.begin(), sorted.end());
  CHECK(sorted == std::vector<std::size_t>{0, 1, 2, 3, 4});
  CHECK(o == sssp::bench::round_order(5, 12345, 0, 0));
  std::set<std::vector<std::size_t>> seen;
  for (std::uint64_t round = 0; round < 20; ++round) {
    seen.insert(sssp::bench::round_order(5, 12345, 0, round));
  }
  CHECK(seen.size() > 1);  // the order changes between rounds
  CHECK(sssp::bench::round_order(1, 3, 0, 0) == std::vector<std::size_t>{0});
  CHECK(sssp::bench::round_order(0, 3, 0, 0).empty());
}

TEST_CASE("bench rss: ru_maxrss is normalised to bytes", "[bench]") {
  CHECK(sssp::bench::maxrss_to_bytes(1000, true) == 1000);
  CHECK(sssp::bench::maxrss_to_bytes(1000, false) == 1024000);
  CHECK_FALSE(sssp::bench::maxrss_to_bytes(0, true).has_value());
  CHECK_FALSE(sssp::bench::maxrss_to_bytes(-1, false).has_value());
  const auto now = sssp::bench::peak_rss_bytes();
  REQUIRE(now.has_value());
  CHECK(*now > 1024);  // a running process uses more than a KiB
}

TEST_CASE("bench env: --env-json supplies exactly the env fields", "[bench]") {
  const Json env = sssp::bench::load_env(env_text(true), true);
  const std::vector<std::string> expected = {
      "paper_version",       "git_sha", "git_dirty",   "os",         "cpu",
      "ram_bytes",           "arch",    "compiler",    "build_type", "build_flags",
      "diagnostics_enabled", "power",   "load_average"};
  CHECK(keys_of(env) == expected);
  CHECK(env.find("ram_bytes")->number_text() == "17179869184");
  CHECK(env.find("load_average")->dump() == "[1.25,1.0,0.5]");

  auto code = [](const std::string& text, bool diag) {
    try {
      sssp::bench::load_env(text, diag);
    } catch (const RunnerError& e) {
      return e.code();
    }
    return 0;
  };
  CHECK(code(env_text(false), true) == sssp::bench::kExitUsage);  // other build's env
  CHECK(code(env_text(true), false) == sssp::bench::kExitUsage);
  CHECK(code("[1]", true) == sssp::bench::kExitUsage);
  CHECK(code("{not json", true) == sssp::bench::kExitUsage);
  std::string missing = env_text(true);
  missing.replace(missing.find("\"arch\""), 6, "\"arcX\"");
  CHECK(code(missing, true) == sssp::bench::kExitUsage);
  std::string wrong_paper = env_text(true);
  wrong_paper.replace(wrong_paper.find("2611.34911v1"), 12, "2611.34911v2");
  CHECK(code(wrong_paper, true) == sssp::bench::kExitUsage);
}

TEST_CASE("bench timestamp: ISO-8601 with an explicit UTC offset", "[bench]") {
  CHECK(sssp::bench::utc_timestamp(std::chrono::system_clock::time_point{}) ==
        "1970-01-01T00:00:00+00:00");
  CHECK(sssp::bench::utc_timestamp(std::chrono::system_clock::time_point{} +
                                   std::chrono::seconds{1764207000}) ==
        "2025-11-27T01:30:00+00:00");
}

TEST_CASE("bench record: all schema keys, env copied, null when unmeasured", "[bench]") {
  const Json env = sssp::bench::load_env(env_text(false), false);
  sssp::bench::TrialRecord t;
  t.algorithm = "dijkstra-binary";
  t.mode = "baseline";
  t.source = 17;
  t.trial = 2;
  t.warmup = true;
  t.order_seed = 12345;
  t.times.graph_gen = 10;
  t.times.sssp = 20;
  t.times.end_to_end = 25;
  t.reachable = 999;
  t.rss_after_graph = 4096;
  const Json r =
      sssp::bench::make_record(env, "RID", "2026-11-27T00:00:00+00:00", graph_info(), t, false);
  CHECK(keys_of(r) == kRecordKeys);
  CHECK(r.find("schema_version")->as_string() == "1");
  CHECK(r.find("run_id")->as_string() == "RID");
  CHECK(r.find("paper_version")->as_string() == "2611.34911v1");
  CHECK(r.find("cpu")->dump() == env.find("cpu")->dump());
  CHECK(r.find("load_average")->dump() == "[1.25,1.0,0.5]");
  CHECK(r.find("params")->dump() == "{}");
  CHECK(r.find("n")->dump() == "1000");
  CHECK(r.find("m")->dump() == "2000");
  CHECK(r.find("avg_degree")->dump() == "4");
  CHECK(r.find("source")->dump() == "17");
  CHECK(r.find("reachable")->dump() == "999");
  CHECK(r.find("trial")->dump() == "2");
  CHECK(r.find("warmup")->dump() == "true");
  CHECK(r.find("order_seed")->dump() == "12345");
  CHECK(r.find("isolated")->dump() == "false");
  CHECK(r.find("status")->as_string() == "ok");
  CHECK(r.find("error")->is_null());
  CHECK(
      r.find("times_ns")->dump() ==
      R"({"graph_gen":10,"preprocess":null,"sssp":20,"tree":null,"end_to_end":25,"verify":null})");
  CHECK(r.find("peak_rss_bytes")->is_null());
  CHECK(r.find("rss_after_graph_bytes")->dump() == "4096");
  CHECK(r.find("ops")->is_null());
  CHECK(r.find("counters")->is_null());
  CHECK(r.find("fallback")->is_null());  // no diagnostics object: nothing recorded
  CHECK(r.find("correctness")->dump() == R"({"verified":false,"ok":null,"method":null})");
}

TEST_CASE("bench record: diagnostics gate ops and counters, fallbacks always", "[bench]") {
  Diagnostics d;
  d.ops.vertex_scans = 5;
  d.ops.queue_pops = 3;
  d.counters["samplings"] = 2;
  d.fallbacks.push_back({"outer: base case", 10, 20});
  d.fallbacks.push_back({"local: small", 1, 2});
  sssp::bench::TrialRecord t;
  t.algorithm = "dijkstra-binary";
  t.mode = "baseline";
  t.times.sssp = 1;
  t.times.end_to_end = 1;
  t.diag = &d;

  const Json off = sssp::bench::make_record(sssp::bench::load_env(env_text(false), false), "R",
                                            "2026-11-27T00:00:00+00:00", graph_info(), t, false);
  CHECK(off.find("ops")->is_null());
  CHECK(off.find("counters")->is_null());
  const std::string fb =
      R"({"count":2,"records":[{"reason":"outer: base case","vertices":10,"edges":20},)"
      R"({"reason":"local: small","vertices":1,"edges":2}],"vertices":11,"edges":22})";
  CHECK(off.find("fallback")->dump() == fb);

  const Json on = sssp::bench::make_record(sssp::bench::load_env(env_text(true), true), "R",
                                           "2026-11-27T00:00:00+00:00", graph_info(), t, true);
  CHECK(on.find("ops")->dump() ==
        R"({"vertex_scans":5,"edge_scans":0,"relaxations":0,"comparisons":0,"additions":0,)"
        R"("queue_pushes":0,"queue_pops":3,"decrease_keys":0})");
  CHECK(on.find("counters")->dump() == R"({"samplings":2})");
  CHECK(on.find("fallback")->dump() == fb);
}

TEST_CASE("bench record: isolation, verification outcome and failures", "[bench]") {
  const Json env = sssp::bench::load_env(env_text(false), false);
  sssp::bench::TrialRecord t;
  t.algorithm = "dijkstra-4ary";
  t.mode = "baseline";
  t.times.sssp = 3;
  t.times.end_to_end = 4;
  t.peak_rss = 1 << 20;

  // Peak RSS belongs to one algorithm only when the process was isolated.
  CHECK(sssp::bench::make_record(env, "R", "2026-11-27T00:00:00+00:00", graph_info(), t, false)
            .find("peak_rss_bytes")
            ->is_null());
  t.isolated = true;
  CHECK(sssp::bench::make_record(env, "R", "2026-11-27T00:00:00+00:00", graph_info(), t, false)
            .find("peak_rss_bytes")
            ->dump() == "1048576");

  t.correctness = {true, true, "reference:dijkstra-binary"};
  t.times.verify = 9;
  Json r = sssp::bench::make_record(env, "R", "2026-11-27T00:00:00+00:00", graph_info(), t, false);
  CHECK(r.find("correctness")->dump() ==
        R"({"verified":true,"ok":true,"method":"reference:dijkstra-binary"})");
  CHECK(r.find("status")->as_string() == "ok");

  // A failed verification can never be reported as ok.
  t.correctness = {true, false, "reference:dijkstra-binary"};
  t.error = "distance mismatch";
  r = sssp::bench::make_record(env, "R", "2026-11-27T00:00:00+00:00", graph_info(), t, false);
  CHECK(r.find("status")->as_string() == "error");
  CHECK(r.find("error")->as_string() == "distance mismatch");

  // A run that threw (e.g. integer overflow) has no source-independent measurements.
  sssp::bench::TrialRecord f;
  f.algorithm = "dijkstra-binary";
  f.mode = "baseline";
  f.status = "error";
  f.error = "integer distance overflow";
  f.times.graph_gen = 5;
  r = sssp::bench::make_record(env, "R", "2026-11-27T00:00:00+00:00", graph_info(), f, false);
  CHECK(r.find("status")->as_string() == "error");
  CHECK(r.find("times_ns")->find("sssp")->is_null());
  CHECK(r.find("reachable")->is_null());

  sssp::bench::GraphInfo gi = graph_info();
  gi.avg_degree.reset();
  CHECK(sssp::bench::make_record(env, "R", "2026-11-27T00:00:00+00:00", gi, t, false)
            .find("avg_degree")
            ->is_null());
}
