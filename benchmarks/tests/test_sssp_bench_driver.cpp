// Tests of the Google Benchmark driver (benchmarks/sssp_bench_driver.hpp, P2-BENCH-CPP).
//
// The driver must (1) offer every algorithm of the bench_runner registry for every weight type it
// has a solver for, on shared inputs, and (2) refuse to time a case whose result is wrong. These
// tests never look at times: they run each case once with a one-iteration minimum and only check
// what was registered, run and verified.

#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/bench_runner.hpp"
#include "sssp/graph.hpp"
#include "sssp_bench_driver.hpp"

namespace {

using sssp::CsrGraph;
using sssp::SsspResult;
using sssp::Vertex;
using sssp::bench::algorithm_registry;
using sssp::bench::AlgorithmEntry;
using sssp::bench::Json;
using sssp::bench::PhaseTimer;
using sssp::bench::gbench::CaseCheck;
using sssp::bench::gbench::DriverConfig;
using sssp::bench::gbench::Suite;
using sssp::bench::gbench::WeightType;

bool starts_with(std::string_view s, std::string_view prefix) {
  return s.substr(0, prefix.size()) == prefix;
}

std::string_view type_name(WeightType t) {
  return t == WeightType::kU64 ? "uint64" : "double";
}

DriverConfig small_config() {
  DriverConfig c;
  c.graph_types = {"er", "grid"};
  c.sizes = {64};
  c.avg_degree = 4.0;
  c.weight_dists = {"uniform_int", "uniform_real"};
  c.seed = 7;
  c.num_sources = 2;
  return c;
}

// A path 0 - 1 - 2 - 3 with weights 1, 2, 3.
CsrGraph<std::uint64_t> path4() {
  return CsrGraph<std::uint64_t>::from_edges(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 3}});
}

SsspResult<std::uint64_t> wrong_distance(const CsrGraph<std::uint64_t>& g, Vertex s,
                                         const Json& /*params*/, PhaseTimer& /*timer*/) {
  auto r = sssp::baseline::dijkstra_binary_heap(g, s);
  r.dist.back() += 1;  // one distance off by one: the reference must catch it
  return r;
}

SsspResult<std::uint64_t> broken_tree(const CsrGraph<std::uint64_t>& g, Vertex s,
                                      const Json& /*params*/, PhaseTimer& /*timer*/) {
  auto r = sssp::baseline::dijkstra_binary_heap(g, s);
  r.parent.back() = sssp::kNoVertex;  // right distances, but no shortest-path tree
  return r;
}

SsspResult<std::uint64_t> overflowing(const CsrGraph<std::uint64_t>& /*g*/, Vertex /*s*/,
                                      const Json& /*params*/, PhaseTimer& /*timer*/) {
  throw std::overflow_error("distance overflow");
}

/// Collects the runs reported by Google Benchmark without printing anything.
class CollectingReporter : public benchmark::BenchmarkReporter {
 public:
  bool ReportContext(const Context& /*context*/) override { return true; }
  void ReportRuns(const std::vector<Run>& runs) override {
    for (const auto& r : runs) {
      names.push_back(r.benchmark_name());
      if (r.skipped != benchmark::internal::NotSkipped) errors.push_back(r.benchmark_name());
      iterations.push_back(r.iterations);
    }
  }
  std::vector<std::string> names;
  std::vector<std::string> errors;
  std::vector<benchmark::IterationCount> iterations;
};

}  // namespace

TEST_CASE("every registered algorithm gets a case for each weight type it supports",
          "[bench][gbench]") {
  const DriverConfig config = small_config();
  Suite suite(config);
  const std::size_t sources = config.num_sources;
  // Two graph types x one size; uniform_int -> uint64, uniform_real -> double.
  REQUIRE(suite.inputs_u64().size() == 2);
  REQUIRE(suite.inputs_f64().size() == 2);

  std::size_t expected_total = 0;
  for (const AlgorithmEntry& e : algorithm_registry()) {
    for (const WeightType t : {WeightType::kU64, WeightType::kF64}) {
      const bool has_solver =
          t == WeightType::kU64 ? e.solve_u64 != nullptr : e.solve_f64 != nullptr;
      const std::size_t inputs =
          t == WeightType::kU64 ? suite.inputs_u64().size() : suite.inputs_f64().size();
      const auto count = static_cast<std::size_t>(std::count_if(
          suite.cases().begin(), suite.cases().end(),
          [&](const Suite::Case& c) { return c.algorithm->name == e.name && c.weight_type == t; }));
      INFO(std::string(e.name) << " x " << type_name(t));
      CHECK(count == (has_solver ? inputs * sources : 0));
      if (has_solver) expected_total += inputs * sources;
    }
  }
  CHECK(suite.cases().size() == expected_total);

  // The radix heap has no double solver and must not appear with double weights.
  CHECK(std::none_of(suite.cases().begin(), suite.cases().end(), [](const Suite::Case& c) {
    return c.algorithm->name == "dijkstra-radix" && c.weight_type == WeightType::kF64;
  }));
}

TEST_CASE("all algorithms of one input share the same graph object and sources",
          "[bench][gbench]") {
  Suite suite(small_config());
  for (const auto& in : suite.inputs_u64()) {
    CHECK(in.sources == sssp::bench::select_sources(in.graph.num_vertices(), in.sources.size(), 7));
  }
  for (const auto& in : suite.inputs_f64()) {
    CHECK(in.sources == sssp::bench::select_sources(in.graph.num_vertices(), in.sources.size(), 7));
  }
  // A case names its input by index, so every algorithm of that input uses one CsrGraph.
  std::set<std::string> labels;
  for (const auto& c : suite.cases()) {
    const std::string& label = c.weight_type == WeightType::kU64
                                   ? suite.inputs_u64()[c.input].label
                                   : suite.inputs_f64()[c.input].label;
    CHECK(c.name ==
          std::string(c.algorithm->name) + "/" + label + "/src:" + std::to_string(c.source_index));
    labels.insert(label);
  }
  CHECK(labels.size() == 4);
}

TEST_CASE("every case agrees with binary-heap Dijkstra on small inputs", "[bench][gbench]") {
  DriverConfig config;
  config.graph_types = {"er", "grid",    "road",     "tree",
                        "ba", "regular", "chung_lu", "adversarial_long_path"};
  config.sizes = {50, 97};
  config.avg_degree = 4.0;
  config.weight_dists = {"uniform_int", "zero_mixed", "wide_range", "uniform_real", "exponential"};
  config.seed = 11;
  config.num_sources = 3;
  Suite suite(config);
  REQUIRE_FALSE(suite.cases().empty());

  std::size_t paper_cases = 0;
  for (std::size_t i = 0; i < suite.cases().size(); ++i) {
    const Suite::Case& c = suite.cases()[i];
    const CaseCheck& chk = suite.check(i);
    INFO(c.name << ": " << chk.message);
    CHECK(chk.status == CaseCheck::Status::kOk);
    if (c.algorithm->mode == "paper-faithful") {
      // Alg 8 l.5: the base case is always recorded as a fallback (CLAUDE.md fallback policy).
      ++paper_cases;
      CHECK(chk.fallback_records >= 1);
      CHECK(chk.fallback_vertices >= 1);
    }
  }
  CHECK(paper_cases > 0);
  CHECK_FALSE(suite.any_failed());

  std::ostringstream log;
  CHECK(suite.verify_all(log));
  CHECK(log.str().find("FAILED") == std::string::npos);
}

TEST_CASE("check_case rejects wrong results instead of timing them", "[bench][gbench]") {
  const auto g = path4();
  const auto reference = sssp::baseline::dijkstra_binary_heap(g, 0).dist;

  AlgorithmEntry wrong{"wrong-distance", "baseline", &wrong_distance, nullptr};
  const CaseCheck a = sssp::bench::gbench::check_case(wrong, g, 0, reference);
  CHECK(a.status == CaseCheck::Status::kFailed);
  CHECK_FALSE(a.message.empty());

  AlgorithmEntry tree{"broken-tree", "baseline", &broken_tree, nullptr};
  CHECK(sssp::bench::gbench::check_case(tree, g, 0, reference).status ==
        CaseCheck::Status::kFailed);

  AlgorithmEntry overflow{"overflow", "baseline", &overflowing, nullptr};
  const CaseCheck c = sssp::bench::gbench::check_case(overflow, g, 0, reference);
  CHECK(c.status == CaseCheck::Status::kFailed);
  CHECK(c.message.find("overflow") != std::string::npos);

  const AlgorithmEntry* dij = sssp::bench::find_algorithm("dijkstra-binary");
  REQUIRE(dij != nullptr);
  CHECK(sssp::bench::gbench::check_case(*dij, g, 0, reference).status == CaseCheck::Status::kOk);
}

TEST_CASE("size limits of the registry are refusals, not results", "[bench][gbench]") {
  AlgorithmEntry limited{"limited", "baseline", nullptr, nullptr};
  limited.max_vertices = 10;
  CHECK_FALSE(sssp::bench::gbench::refusal_reason(limited, 10, 100).has_value());
  CHECK(sssp::bench::gbench::refusal_reason(limited, 11, 1).has_value());

  AlgorithmEntry work{"work", "baseline", nullptr, nullptr};
  work.max_work = 100;
  CHECK_FALSE(sssp::bench::gbench::refusal_reason(work, 10, 10).has_value());
  CHECK(sssp::bench::gbench::refusal_reason(work, 10, 11).has_value());

  const AlgorithmEntry* faithful = sssp::bench::find_algorithm("paper-faithful");
  REQUIRE(faithful != nullptr);
  CHECK(sssp::bench::gbench::refusal_reason(*faithful, sssp::bench::kScanBaseCaseMaxVertices + 1, 0)
            .has_value());
}

TEST_CASE("registered benchmarks cover the registry and run once without error",
          "[bench][gbench]") {
  Suite suite(small_config());
  benchmark::ClearRegisteredBenchmarks();
  const std::vector<std::string> names = suite.register_benchmarks();
  REQUIRE(names.size() == suite.cases().size());
  CHECK(std::set<std::string>(names.begin(), names.end()).size() == names.size());
  for (const AlgorithmEntry& e : algorithm_registry()) {
    INFO(e.name);
    CHECK(std::any_of(names.begin(), names.end(), [&](const std::string& n) {
      return starts_with(n, std::string(e.name) + "/");
    }));
  }

  // One iteration per benchmark: this checks the wiring, not the speed.
  std::array<char, 6> prog{"tests"};
  std::array<char, 26> min_time{"--benchmark_min_time=1x"};
  std::array<char*, 2> argv{prog.data(), min_time.data()};
  int argc = static_cast<int>(argv.size());
  benchmark::Initialize(&argc, argv.data());
  CollectingReporter reporter;
  const std::size_t ran = benchmark::RunSpecifiedBenchmarks(&reporter);
  benchmark::ClearRegisteredBenchmarks();

  CHECK(ran == names.size());
  CHECK(reporter.errors.empty());
  for (const auto& n : names) {
    INFO(n);
    // Manual time appends "/manual_time" to the reported name.
    CHECK(std::find(reporter.names.begin(), reporter.names.end(), n + "/manual_time") !=
          reporter.names.end());
  }
  CHECK(std::all_of(reporter.iterations.begin(), reporter.iterations.end(),
                    [](benchmark::IterationCount it) { return it == 1; }));
  CHECK_FALSE(suite.any_failed());
}

TEST_CASE("driver arguments are parsed and validated", "[bench][gbench]") {
  const std::vector<std::string> args{"--sssp_graphs=er,road",
                                      "--sssp_sizes=100,1000",
                                      "--sssp_avg_degree=8",
                                      "--sssp_weights=uniform_real",
                                      "--sssp_seed=5",
                                      "--sssp_sources=3",
                                      "--sssp_algorithms=dijkstra-binary,paper-faithful-d1",
                                      "--sssp_verify_only"};
  const DriverConfig c = sssp::bench::gbench::parse_driver_args(args);
  CHECK(c.graph_types == std::vector<std::string>{"er", "road"});
  CHECK(c.sizes == std::vector<std::size_t>{100, 1000});
  CHECK(c.avg_degree == 8.0);
  CHECK(c.weight_dists == std::vector<std::string>{"uniform_real"});
  CHECK(c.seed == 5);
  CHECK(c.num_sources == 3);
  CHECK(c.algorithms == std::vector<std::string>{"dijkstra-binary", "paper-faithful-d1"});
  CHECK(c.verify_only);

  using sssp::bench::gbench::parse_driver_args;
  CHECK_THROWS_AS(parse_driver_args(std::vector<std::string>{"--sssp_bogus=1"}),
                  std::invalid_argument);
  CHECK_THROWS_AS(parse_driver_args(std::vector<std::string>{"--sssp_algorithms=paper"}),
                  std::invalid_argument);
  CHECK_THROWS_AS(parse_driver_args(std::vector<std::string>{"--sssp_weights=normal"}),
                  std::invalid_argument);
  CHECK_THROWS_AS(parse_driver_args(std::vector<std::string>{"--sssp_sizes=10,x"}),
                  std::invalid_argument);
  CHECK_THROWS_AS(parse_driver_args(std::vector<std::string>{"--sssp_sources=0"}),
                  std::invalid_argument);
  CHECK_THROWS_AS(parse_driver_args(std::vector<std::string>{"--sssp_graphs="}),
                  std::invalid_argument);

  CHECK(sssp::bench::gbench::weight_type_of("uniform_int") == "uint64");
  CHECK(sssp::bench::gbench::weight_type_of("builtin") == "uint64");
  CHECK(sssp::bench::gbench::weight_type_of("exponential") == "double");
  CHECK(sssp::bench::gbench::has_builtin_weights("road"));
  CHECK_FALSE(sssp::bench::gbench::has_builtin_weights("er"));
}
