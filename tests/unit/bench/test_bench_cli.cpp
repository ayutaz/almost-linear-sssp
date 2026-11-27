// Argument parsing of sssp_bench_runner (docs/design/benchmark_runner_cli.md §2) and the
// unsupported-combination check (exit code 3, §3.7 and §4).

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <initializer_list>
#include <string>
#include <vector>

#include "sssp/bench_runner.hpp"

using sssp::bench::check_supported;
using sssp::bench::kExitUnsupported;
using sssp::bench::kExitUsage;
using sssp::bench::nominal_avg_degree;
using sssp::bench::Options;
using sssp::bench::parse_args;
using sssp::bench::RunnerError;
using sssp::bench::VerifyMode;

namespace {

/// Minimal valid invocation plus `extra`.
std::vector<std::string> base_args(std::initializer_list<std::string> extra = {}) {
  std::vector<std::string> a = {
      "--type",          "er",         "--n",      "100",      "--algorithms",
      "dijkstra-binary", "--env-json", "env.json", "--run-id", "RID"};
  a.insert(a.end(), extra.begin(), extra.end());
  return a;
}

int exit_code_of(const std::vector<std::string>& args) {
  try {
    const Options o = parse_args(args);
    check_supported(o);
  } catch (const RunnerError& e) {
    return e.code();
  }
  return 0;
}

}  // namespace

TEST_CASE("bench cli: the documented example invocation", "[bench]") {
  const Options o = parse_args(std::vector<std::string>{
      "--type",    "er",          "--n",           "100000",        "--avg-degree",  "4",
      "--weights", "uniform_int", "--weight-type", "uint64",        "--seed",        "1",
      "--sources", "2",           "--algorithms",  "dijkstra-4ary", "--trials",      "5",
      "--warmup",  "1",           "--order-seed",  "12345",         "--timeout-sec", "600",
      "--verify",  "reference",   "--env-json",    "e.json",        "--run-id",      "R1",
      "--out",     "o.jsonl",     "--isolate"});
  CHECK(o.spec.type == "er");
  CHECK(o.spec.n == 100000);
  CHECK(o.spec.avg_degree == 4.0);
  CHECK(o.weight_dist == "uniform_int");
  CHECK_FALSE(o.use_double);
  CHECK(o.spec.options.seed == 1);
  CHECK(o.num_sources == 2);
  CHECK(o.algorithms == std::vector<std::string>{"dijkstra-4ary"});
  CHECK(o.trials == 5);
  CHECK(o.warmup == 1);
  CHECK(o.order_seed == 12345);
  REQUIRE(o.timeout_sec.has_value());
  CHECK(*o.timeout_sec == 600.0);
  CHECK(o.verify == VerifyMode::kReference);
  CHECK(o.env_json == "e.json");
  CHECK(o.run_id == "R1");
  CHECK(o.out == "o.jsonl");
  CHECK(o.isolate);
  CHECK(o.params.items().empty());
}

TEST_CASE("bench cli: defaults follow the CLI table", "[bench]") {
  const Options o = parse_args(base_args());
  CHECK(o.num_sources == 1);
  CHECK(o.source_list.empty());
  CHECK(o.trials == 5);
  CHECK(o.warmup == 1);
  CHECK(o.order_seed == 0);
  CHECK_FALSE(o.timeout_sec.has_value());
  CHECK(o.verify == VerifyMode::kReference);
  CHECK_FALSE(o.isolate);
  CHECK(o.out.empty());
  CHECK(o.weight_dist == "uniform_int");
}

TEST_CASE("bench cli: generator options use sssp_gen names", "[bench]") {
  const Options o = parse_args(base_args({"--weights",   "zero_mixed",  "--lo",
                                          "2",           "--hi",        "9",
                                          "--zero-prob", "0.5",         "--weight-type",
                                          "double",      "--connected", "--degree",
                                          "6",           "--attach",    "3",
                                          "--exponent",  "2.2",         "--wide-log2",
                                          "40",          "--mean",      "3",
                                          "--rows",      "4",           "--cols",
                                          "5",           "--multi",     "keep"}));
  CHECK(o.weight_dist == "zero_mixed");
  CHECK(o.spec.options.weights.lo == 2.0);
  CHECK(o.spec.options.weights.hi == 9.0);
  CHECK(o.spec.options.weights.zero_probability == 0.5);
  CHECK(o.spec.options.weights.wide_log2 == 40);
  CHECK(o.spec.options.weights.mean == 3.0);
  CHECK(o.use_double);
  CHECK(o.spec.options.connected);
  CHECK(o.spec.degree == 6);
  CHECK(o.spec.attach == 3);
  CHECK(o.spec.exponent == 2.2);
  CHECK(o.spec.rows == 4);
  CHECK(o.spec.cols == 5);
}

TEST_CASE("bench cli: multiple algorithms, source list and verify modes", "[bench]") {
  Options o = parse_args(std::vector<std::string>{
      "--type", "grid", "--n", "100", "--algorithms", "dijkstra-binary,dijkstra-radix,bellman-ford",
      "--source-list", "0,5,99", "--verify", "bellman-ford", "--env-json", "e", "--run-id", "r"});
  CHECK(o.algorithms ==
        std::vector<std::string>{"dijkstra-binary", "dijkstra-radix", "bellman-ford"});
  CHECK(o.source_list == std::vector<sssp::Vertex>{0, 5, 99});
  CHECK(o.verify == VerifyMode::kBellmanFord);
  o = parse_args(base_args({"--verify", "none"}));
  CHECK(o.verify == VerifyMode::kNone);
}

TEST_CASE("bench cli: --param values are JSON when possible, text otherwise", "[bench]") {
  const Options o = parse_args(std::vector<std::string>{
      "--type", "er", "--algorithms", "paper", "--param", "threshold=4", "--param", "name=abc",
      "--param", "flag=true", "--param", "ratio=0.5", "--env-json", "e", "--run-id", "r"});
  REQUIRE(o.params.items().size() == 4);
  CHECK(o.params.find("threshold")->dump() == "4");
  CHECK(o.params.find("name")->dump() == "\"abc\"");
  CHECK(o.params.find("flag")->dump() == "true");
  CHECK(o.params.find("ratio")->dump() == "0.5");
}

TEST_CASE("bench cli: argument errors exit with code 2", "[bench]") {
  CHECK(exit_code_of(base_args({"--bogus"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--seed"})) == kExitUsage);      // missing value
  CHECK(exit_code_of(base_args({"--n", "abc"})) == kExitUsage);  // not a number
  CHECK(exit_code_of(base_args({"--n", "-5"})) == kExitUsage);   // negative count
  CHECK(exit_code_of(base_args({"--n", "12x"})) == kExitUsage);  // trailing garbage
  CHECK(exit_code_of(base_args({"--weights", "gauss"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--weight-type", "float"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--verify", "maybe"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--trials", "0"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--sources", "0"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--timeout-sec", "0"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--sources", "2", "--source-list", "1,2"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--source-list", "1,,2"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--param", "noequals"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--param", "=3"})) == kExitUsage);
  // Baselines take no parameters.
  CHECK(exit_code_of(base_args({"--param", "k=1"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--weights", "uniform_int", "--lo", "5", "--hi", "1"})) ==
        kExitUsage);
  CHECK(exit_code_of(base_args({"--type", "hypercube"})) == kExitUsage);
  // --weights builtin is only for families whose weights are fixed by construction.
  CHECK(exit_code_of(base_args({"--weights", "builtin"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--type", "road", "--weights", "uniform_int"})) == kExitUsage);
}

TEST_CASE("bench cli: required options and --isolate", "[bench]") {
  CHECK(exit_code_of({"--type", "er", "--env-json", "e", "--run-id", "r"}) == kExitUsage);
  CHECK(exit_code_of({"--type", "er", "--algorithms", "dijkstra-binary", "--run-id", "r"}) ==
        kExitUsage);
  CHECK(exit_code_of({"--type", "er", "--algorithms", "dijkstra-binary", "--env-json", "e"}) ==
        kExitUsage);
  CHECK(exit_code_of({"--algorithms", "dijkstra-binary", "--env-json", "e", "--run-id", "r"}) ==
        kExitUsage);  // no --type and no --graph-file
  CHECK(exit_code_of(base_args({"--algorithms", "dijkstra-binary,dijkstra-binary"})) == kExitUsage);
  CHECK(exit_code_of(base_args({"--isolate"})) == 0);
  CHECK(exit_code_of(base_args({"--algorithms", "dijkstra-binary,dijkstra-4ary", "--isolate"})) ==
        kExitUsage);
  // --list / --version / --help need nothing else.
  CHECK(parse_args(std::vector<std::string>{"--list"}).list);
  CHECK(parse_args(std::vector<std::string>{"--version"}).version);
  CHECK(parse_args(std::vector<std::string>{"--help"}).help);
}

TEST_CASE("bench cli: builtin-weight families", "[bench]") {
  CHECK(parse_args(base_args({"--type", "road"})).weight_dist == "builtin");
  CHECK(parse_args(base_args({"--type", "adversarial_long_path", "--weights", "builtin"}))
            .weight_dist == "builtin");
  CHECK(parse_args(base_args({"--type", "tree"})).weight_dist == "uniform_int");
}

TEST_CASE("bench cli: unsupported combinations exit with code 3", "[bench]") {
  CHECK(exit_code_of(base_args({"--algorithms", "dijkstra-radix", "--weight-type", "double"})) ==
        kExitUnsupported);
  CHECK(exit_code_of(base_args({"--algorithms", "dijkstra-radix"})) == 0);
  CHECK(exit_code_of(base_args({"--algorithms", "dijkstra-binary", "--weight-type", "double"})) ==
        0);
  // Known to the schema but not implemented yet, and unknown names.
  CHECK(exit_code_of(base_args({"--algorithms", "paper"})) == kExitUnsupported);
  CHECK(exit_code_of(base_args({"--algorithms", "paper-optimized"})) == kExitUnsupported);
  CHECK(exit_code_of(base_args({"--algorithms", "dijkstra-binary,no-such-algo"})) ==
        kExitUnsupported);
  try {
    check_supported(parse_args(base_args({"--algorithms", "paper"})));
    FAIL("expected RunnerError");
  } catch (const RunnerError& e) {
    CHECK(std::string(e.what()).find("not implemented") != std::string::npos);
  }
}

TEST_CASE("bench cli: nominal average degree written to records", "[bench]") {
  auto deg = [](std::initializer_list<std::string> extra) {
    return nominal_avg_degree(parse_args(base_args(extra)));
  };
  CHECK(deg({"--avg-degree", "8"}) == 8.0);
  CHECK(deg({"--type", "chung_lu", "--avg-degree", "16"}) == 16.0);
  CHECK(deg({"--type", "regular", "--degree", "6"}) == 6.0);
  CHECK(deg({"--type", "ba", "--attach", "3"}) == 6.0);
  CHECK_FALSE(deg({"--type", "grid"}).has_value());
  CHECK_FALSE(deg({"--type", "road"}).has_value());
  CHECK_FALSE(deg({"--type", "tree"}).has_value());
  CHECK_FALSE(deg({"--type", "adversarial_long_path"}).has_value());
}
