// End-to-end behaviour of the runner (run_main) on small generated graphs:
// docs/design/benchmark_runner_cli.md §3 (behaviour) and §4 (exit codes).

#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "sssp/bench_runner.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph_io.hpp"

namespace fs = std::filesystem;
using sssp::bench::Json;

namespace {

/// Fresh scratch directory under the system temp dir (removed at scope exit).
class TempDir {
 public:
  TempDir() {
    static std::atomic<int> counter{0};
    path_ = fs::temp_directory_path() /
            ("sssp_bench_run_test_" + std::to_string(::getpid()) + "_" + std::to_string(counter++));
    fs::create_directories(path_);
  }
  ~TempDir() {
    std::error_code ec;
    fs::remove_all(path_, ec);
  }
  TempDir(const TempDir&) = delete;
  TempDir& operator=(const TempDir&) = delete;
  [[nodiscard]] std::string file(const std::string& name) const { return (path_ / name).string(); }

 private:
  fs::path path_;
};

std::string write_env(const TempDir& dir) {
  const std::string path = dir.file("env.json");
  std::ofstream(path) << R"({"paper_version": "2611.34911v1", "git_sha": "abcdef0",
    "git_dirty": false, "os": "TestOS", "cpu": {"brand": null, "logical_cores": null,
    "physical_cores": null, "perf_cores": null, "efficiency_cores": null}, "ram_bytes": null,
    "arch": "test", "compiler": null, "build_type": "Debug", "build_flags": null,
    "diagnostics_enabled": )"
                      << (sssp::bench::kDiagnosticsCompiled ? "true" : "false")
                      << R"(, "power": null, "load_average": null})";
  return path;
}

struct Run {
  int code = -1;
  std::string out;
  std::string err;
};

Run run(const std::vector<std::string>& args) {
  std::ostringstream out;
  std::ostringstream err;
  Run r;
  r.code = sssp::bench::run_main(args, out, err);
  r.out = out.str();
  r.err = err.str();
  return r;
}

std::vector<Json> read_records(const std::string& path) {
  std::vector<Json> records;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty()) records.push_back(Json::parse(line));
  }
  return records;
}

std::vector<std::string> args_for(const TempDir& dir, const std::string& out,
                                  std::initializer_list<std::string> extra) {
  std::vector<std::string> a = {"--env-json", write_env(dir), "--run-id", "TEST-RUN", "--out", out};
  a.insert(a.end(), extra.begin(), extra.end());
  return a;
}

std::uint64_t as_u64(const Json& j) {
  return std::stoull(j.number_text());
}

}  // namespace

TEST_CASE("bench run: records for every source x round x algorithm", "[bench]") {
  TempDir dir;
  const std::string out = dir.file("r.jsonl");
  const auto r =
      run(args_for(dir, out, {"--type",       "er",
                              "--n",          "300",
                              "--avg-degree", "4",
                              "--weights",    "uniform_int",
                              "--seed",       "3",
                              "--sources",    "2",
                              "--algorithms", "dijkstra-binary,dijkstra-4ary,dijkstra-fibonacci",
                              "--trials",     "2",
                              "--warmup",     "1",
                              "--order-seed", "12345"}));
  INFO(r.err);
  REQUIRE(r.code == 0);
  const auto recs = read_records(out);
  REQUIRE(recs.size() == 2 * (1 + 2) * 3);

  const auto g = sssp::gen::generate<std::uint64_t>([] {
    sssp::gen::GraphSpec s;
    s.type = "er";
    s.n = 300;
    s.avg_degree = 4;
    s.options.seed = 3;
    return s;
  }());
  const auto sources = sssp::bench::select_sources(300, 2, 3);
  const std::set<std::string> graph_gen_values = [&] {
    std::set<std::string> v;
    for (const auto& rec : recs) v.insert(rec.find("times_ns")->find("graph_gen")->dump());
    return v;
  }();
  CHECK(graph_gen_values.size() == 1);  // generated once, shared by every run

  std::size_t warmups = 0;
  for (std::size_t i = 0; i < recs.size(); ++i) {
    const Json& rec = recs[i];
    INFO(rec.dump());
    CHECK(rec.items().size() == 41);
    CHECK(rec.find("run_id")->as_string() == "TEST-RUN");
    CHECK(rec.find("mode")->as_string() == "baseline");
    CHECK(rec.find("graph_type")->as_string() == "er");
    CHECK(as_u64(*rec.find("n")) == 300);
    CHECK(as_u64(*rec.find("m")) == g.edges.size());
    CHECK(rec.find("avg_degree")->dump() == "4");
    CHECK(rec.find("weight_dist")->as_string() == "uniform_int");
    CHECK(rec.find("weight_type")->as_string() == "uint64");
    CHECK(as_u64(*rec.find("seed")) == 3);
    CHECK(as_u64(*rec.find("order_seed")) == 12345);
    CHECK(rec.find("isolated")->dump() == "false");
    CHECK(rec.find("peak_rss_bytes")->is_null());
    CHECK(rec.find("status")->as_string() == "ok");
    CHECK(rec.find("correctness")->dump() ==
          R"({"verified":true,"ok":true,"method":"reference:dijkstra-binary"})");
    const Json& t = *rec.find("times_ns");
    CHECK(t.find("preprocess")->is_null());
    CHECK(t.find("tree")->is_null());
    CHECK(as_u64(*t.find("end_to_end")) >= as_u64(*t.find("sssp")));
    CHECK(as_u64(*t.find("verify")) >= 1);
    CHECK(as_u64(*rec.find("reachable")) <= 300);
    CHECK(rec.find("ops")->is_null() == !sssp::bench::kDiagnosticsCompiled);
    CHECK(rec.find("counters")->is_null() == !sssp::bench::kDiagnosticsCompiled);
    CHECK(rec.find("fallback")->dump() == R"({"count":0,"records":[],"vertices":0,"edges":0})");
    // Source-major order: all rounds of source 0, then of source 1.
    CHECK(as_u64(*rec.find("source")) == sources[i / 9]);
    const std::size_t round = (i % 9) / 3;
    CHECK(rec.find("warmup")->as_bool() == (round == 0));
    CHECK(as_u64(*rec.find("trial")) == (round == 0 ? 0 : round - 1));
    if (rec.find("warmup")->as_bool()) ++warmups;
  }
  CHECK(warmups == 2 * 3);

  // Within each round the algorithms run in the seeded order.
  const std::vector<std::string> algos = {"dijkstra-binary", "dijkstra-4ary", "dijkstra-fibonacci"};
  for (std::size_t si = 0; si < 2; ++si) {
    for (std::size_t round = 0; round < 3; ++round) {
      const auto order = sssp::bench::round_order(3, 12345, si, round);
      for (std::size_t k = 0; k < 3; ++k) {
        CHECK(recs[si * 9 + round * 3 + k].find("algorithm")->as_string() == algos[order[k]]);
      }
    }
  }
}

TEST_CASE("bench run: output is appended, one record per line", "[bench]") {
  TempDir dir;
  const std::string out = dir.file("r.jsonl");
  const auto a = args_for(dir, out,
                          {"--type", "grid", "--n", "64", "--algorithms", "dijkstra-radix",
                           "--trials", "1", "--warmup", "0"});
  REQUIRE(run(a).code == 0);
  REQUIRE(run(a).code == 0);
  const auto recs = read_records(out);
  REQUIRE(recs.size() == 2);
  CHECK(recs[0].find("avg_degree")->is_null());  // grid: degree fixed by the family
  CHECK(as_u64(*recs[0].find("n")) == 64);
}

TEST_CASE("bench run: --isolate records one peak RSS for all records", "[bench]") {
  TempDir dir;
  const std::string out = dir.file("r.jsonl");
  const auto r =
      run(args_for(dir, out,
                   {"--type", "ba", "--n", "200", "--attach", "2", "--weights", "uniform_real",
                    "--weight-type", "double", "--algorithms", "dijkstra-fibonacci", "--trials",
                    "3", "--warmup", "1", "--isolate"}));
  INFO(r.err);
  REQUIRE(r.code == 0);
  const auto recs = read_records(out);
  REQUIRE(recs.size() == 4);
  std::set<std::string> peaks;
  for (const auto& rec : recs) {
    CHECK(rec.find("isolated")->dump() == "true");
    CHECK(rec.find("weight_type")->as_string() == "double");
    CHECK(rec.find("avg_degree")->dump() == "4");  // ba: 2 x attach
    CHECK(rec.find("status")->as_string() == "ok");
    REQUIRE_FALSE(rec.find("peak_rss_bytes")->is_null());
    CHECK(as_u64(*rec.find("peak_rss_bytes")) >= as_u64(*rec.find("rss_after_graph_bytes")));
    peaks.insert(rec.find("peak_rss_bytes")->dump());
  }
  CHECK(peaks.size() == 1);
}

TEST_CASE("bench run: verification modes", "[bench]") {
  TempDir dir;
  const std::string out = dir.file("r.jsonl");
  REQUIRE(run(args_for(dir, out,
                       {"--type", "tree", "--n", "100", "--algorithms", "dijkstra-4ary", "--trials",
                        "1", "--warmup", "0", "--verify", "none"}))
              .code == 0);
  REQUIRE(run(args_for(dir, out,
                       {"--type", "road", "--n", "100", "--algorithms", "bellman-ford", "--trials",
                        "1", "--warmup", "0", "--verify", "bellman-ford"}))
              .code == 0);
  const auto recs = read_records(out);
  REQUIRE(recs.size() == 2);
  CHECK(recs[0].find("correctness")->dump() == R"({"verified":false,"ok":null,"method":null})");
  CHECK(recs[0].find("times_ns")->find("verify")->is_null());
  CHECK(recs[1].find("correctness")->dump() ==
        R"({"verified":true,"ok":true,"method":"bellman-ford"})");
  CHECK(recs[1].find("weight_dist")->as_string() == "builtin");
}

TEST_CASE("bench run: integer overflow is an error record, not a crash", "[bench]") {
  TempDir dir;
  const std::string out = dir.file("r.jsonl");
  // A 5000-vertex path with weight 2^53: distances beyond 2048 hops exceed 2^64.
  const auto r =
      run(args_for(dir, out, {"--type",        "grid",     "--rows",       "1",
                              "--cols",        "5000",     "--n",          "5000",
                              "--weights",     "constant", "--lo",         "9007199254740992",
                              "--source-list", "0",        "--algorithms", "dijkstra-binary",
                              "--trials",      "1",        "--warmup",     "0"}));
  INFO(r.err);
  REQUIRE(r.code == 0);
  const auto recs = read_records(out);
  REQUIRE(recs.size() == 1);
  CHECK(recs[0].find("status")->as_string() == "error");
  CHECK(recs[0].find("error")->as_string().find("overflow") != std::string::npos);
  CHECK(recs[0].find("times_ns")->find("sssp")->is_null());
  CHECK(recs[0].find("times_ns")->find("end_to_end")->is_null());
  CHECK_FALSE(recs[0].find("times_ns")->find("graph_gen")->is_null());
  CHECK(recs[0].find("correctness")->find("verified")->dump() == "false");
}

TEST_CASE("bench run: --graph-file reads a generated edge list", "[bench]") {
  TempDir dir;
  sssp::gen::GraphSpec s;
  s.type = "chung_lu";
  s.n = 150;
  s.avg_degree = 6;
  s.options.seed = 9;
  const auto el = sssp::gen::generate<std::uint64_t>(s);
  const std::string gpath = dir.file("g.txt");
  sssp::io::write_edge_list_file(gpath, el);
  const std::string out = dir.file("r.jsonl");
  const auto r =
      run(args_for(dir, out,
                   {"--graph-file", gpath, "--type", "chung_lu", "--avg-degree", "6",
                    "--algorithms", "dijkstra-binary", "--trials", "1", "--warmup", "0"}));
  INFO(r.err);
  REQUIRE(r.code == 0);
  const auto recs = read_records(out);
  REQUIRE(recs.size() == 1);
  CHECK(as_u64(*recs[0].find("n")) == el.num_vertices);
  CHECK(as_u64(*recs[0].find("m")) == el.edges.size());
  CHECK(recs[0].find("graph_type")->as_string() == "chung_lu");
  CHECK(recs[0].find("avg_degree")->dump() == "6");
}

TEST_CASE("bench run: Mode A records carry the paper base-case FallbackRecord", "[bench]") {
  for (const bool use_double : {false, true}) {
    TempDir dir;
    const std::string out = dir.file("r.jsonl");
    auto args =
        args_for(dir, out, {"--type",        "er",
                            "--n",           "120",
                            "--avg-degree",  "3",
                            "--weights",     use_double ? "uniform_real" : "zero_mixed",
                            "--weight-type", use_double ? "double" : "uint64",
                            "--seed",        "5",
                            "--sources",     "2",
                            "--algorithms",  "paper-faithful,paper-faithful-d1,dijkstra-binary",
                            "--trials",      "1",
                            "--warmup",      "1"});
    const auto r = run(args);
    INFO(r.err);
    REQUIRE(r.code == 0);
    const auto recs = read_records(out);
    REQUIRE(recs.size() == 2 * 2 * 3);
    std::set<std::string> seen;
    for (const Json& rec : recs) {
      INFO(rec.dump());
      const std::string algo = rec.find("algorithm")->as_string();
      seen.insert(algo);
      CHECK(rec.find("status")->as_string() == "ok");
      CHECK(rec.find("correctness")->dump() ==
            R"({"verified":true,"ok":true,"method":"reference:dijkstra-binary"})");
      const Json& fb = *rec.find("fallback");
      if (algo == "dijkstra-binary") {
        CHECK(rec.find("mode")->as_string() == "baseline");
        CHECK(as_u64(*fb.find("count")) == 0);
        continue;
      }
      CHECK(rec.find("mode")->as_string() == "paper-faithful");
      // Exactly one record (Alg 8 l.5 base case), written also when diagnostics are off.
      REQUIRE(as_u64(*fb.find("count")) == 1);
      const Json& f = fb.find("records")->elements().at(0);
      const std::string expected_reason =
          algo == "paper-faithful"
              ? "outer: N <= L_N^64 (paper base case, scan Dijkstra)"
              : "outer: N <= L_N^64 (paper base case, binary-heap Dijkstra, deviation D-1)";
      CHECK(f.find("reason")->as_string() == expected_reason);
      // The fallback share is the contracted source component: at most (n, m), at least the
      // source; the totals equal the single record.
      CHECK(as_u64(*f.find("vertices")) >= 1);
      CHECK(as_u64(*f.find("vertices")) <= as_u64(*rec.find("n")));
      CHECK(as_u64(*f.find("edges")) <= as_u64(*rec.find("m")));
      CHECK(as_u64(*fb.find("vertices")) == as_u64(*f.find("vertices")));
      CHECK(as_u64(*fb.find("edges")) == as_u64(*f.find("edges")));
      // paper_sssp is one call: contraction, base case and tree restoration are all in sssp.
      const Json& t = *rec.find("times_ns");
      CHECK(t.find("preprocess")->is_null());
      CHECK(t.find("tree")->is_null());
      CHECK(as_u64(*t.find("sssp")) >= 1);
    }
    CHECK(seen == std::set<std::string>{"paper-faithful", "paper-faithful-d1", "dijkstra-binary"});
  }
}

TEST_CASE("bench run: the scan base case is refused above its vertex limit", "[bench]") {
  TempDir dir;
  const std::string out = dir.file("r.jsonl");
  const std::string too_many = std::to_string(sssp::bench::kScanBaseCaseMaxVertices + 1);
  const auto refused = run(args_for(
      dir, out,
      {"--type", "tree", "--n", too_many, "--algorithms", "paper-faithful-d1,paper-faithful",
       "--trials", "1", "--warmup", "0", "--verify", "none"}));
  CHECK(refused.code == sssp::bench::kExitUnsupported);
  CHECK(refused.err.find("paper-faithful is limited to n <= " +
                         std::to_string(sssp::bench::kScanBaseCaseMaxVertices)) !=
        std::string::npos);
  CHECK(refused.err.find(too_many) != std::string::npos);
  CHECK_FALSE(fs::exists(out));  // refused before the first record, nothing substituted

  // The binary-heap variant has no limit: the same graph runs.
  const auto d1 =
      run(args_for(dir, out,
                   {"--type", "tree", "--n", too_many, "--algorithms", "paper-faithful-d1",
                    "--trials", "1", "--warmup", "0", "--verify", "none"}));
  INFO(d1.err);
  CHECK(d1.code == 0);
  CHECK(read_records(out).size() == 1);
}

TEST_CASE("bench run: exit codes and nothing written on refusal", "[bench]") {
  TempDir dir;
  const std::string out = dir.file("r.jsonl");
  auto code = [&](std::initializer_list<std::string> extra) {
    const auto r = run(args_for(dir, out, extra));
    CHECK_FALSE(r.err.empty());  // every refusal explains itself
    return r.code;
  };
  using sssp::bench::kExitUnsupported;
  using sssp::bench::kExitUsage;
  CHECK(code({"--type", "er", "--n", "50", "--algorithms", "dijkstra-radix", "--weight-type",
              "double", "--weights", "uniform_real"}) == kExitUnsupported);
  CHECK(code({"--type", "er", "--n", "50", "--algorithms", "paper"}) == kExitUnsupported);
  CHECK(code({"--type", "er", "--n", "50", "--algorithms", "no-such"}) == kExitUnsupported);
  // Bellman-Ford is an oracle for small graphs only (n * m > kBellmanFordMaxWork).
  CHECK(code({"--type", "er", "--n", "40000", "--avg-degree", "2", "--algorithms", "bellman-ford",
              "--trials", "1", "--warmup", "0"}) == kExitUnsupported);
  CHECK(code({"--type", "er", "--n", "40000", "--avg-degree", "2", "--algorithms",
              "dijkstra-binary", "--verify", "bellman-ford", "--trials", "1", "--warmup", "0"}) ==
        kExitUnsupported);
  CHECK(code({"--type", "er", "--n", "50", "--algorithms", "dijkstra-binary", "--bogus"}) ==
        kExitUsage);
  CHECK(code({"--type", "er", "--n", "50", "--algorithms", "dijkstra-binary", "--source-list",
              "50"}) == kExitUsage);  // source out of range
  CHECK(code({"--type", "er", "--n", "5", "--avg-degree", "1", "--algorithms", "dijkstra-binary",
              "--sources", "6"}) == kExitUsage);  // more sources than vertices
  CHECK(code({"--type", "regular", "--n", "5", "--degree", "3", "--algorithms",
              "dijkstra-binary"}) == kExitUsage);  // n * d odd: no such graph
  CHECK(code({"--graph-file", dir.file("missing.txt"), "--algorithms", "dijkstra-binary"}) ==
        kExitUsage);
  CHECK_FALSE(fs::exists(out));

  const auto bad_env = run({"--env-json", dir.file("nope.json"), "--run-id", "R", "--out", out,
                            "--type", "er", "--algorithms", "dijkstra-binary"});
  CHECK(bad_env.code == kExitUsage);
  CHECK_FALSE(fs::exists(out));
}

TEST_CASE("bench run: --timeout-sec stops before starting new runs", "[bench]") {
  TempDir dir;
  const std::string out = dir.file("r.jsonl");
  const auto r = run(args_for(dir, out,
                              {"--type", "er", "--n", "100", "--algorithms", "dijkstra-binary",
                               "--trials", "1000", "--warmup", "0", "--timeout-sec", "1e-9"}));
  CHECK(r.code == 0);
  CHECK(r.err.find("timeout") != std::string::npos);
  CHECK(read_records(out).size() < 1000);
}

TEST_CASE("bench run: --list and --version", "[bench]") {
  const auto l = run({"--list"});
  REQUIRE(l.code == 0);
  const Json j = Json::parse(l.out);
  std::set<std::string> types;
  for (const auto& t : j.find("graph_types")->elements()) types.insert(t.as_string());
  for (const auto& g : sssp::gen::generator_catalog()) CHECK(types.count(std::string(g.name)) == 1);
  std::set<std::string> dists;
  for (const auto& d : j.find("weight_dists")->elements()) dists.insert(d.as_string());
  CHECK(dists == std::set<std::string>{"constant", "uniform_int", "uniform_real", "log_uniform",
                                       "exponential", "zero_mixed", "wide_range", "builtin"});
  const auto& algos = j.find("algorithms")->elements();
  CHECK(algos.size() == sssp::bench::algorithm_registry().size());
  for (const auto& a : algos) {
    const std::string name = a.find("name")->as_string();
    CHECK(a.find("mode")->as_string() == *sssp::bench::schema_mode(name));
    CHECK(a.find("weight_types")->dump() == ((name == "dijkstra-radix" || name == "practical-local")
                                                 ? R"(["uint64"])"
                                                 : R"(["uint64","double"])"));
    // Size limits (exit 3 above them); null when the algorithm has none.
    const auto* e = sssp::bench::find_algorithm(name);
    REQUIRE(e != nullptr);
    CHECK(a.find("max_work")->dump() ==
          (e->max_work == 0 ? std::string("null") : std::to_string(e->max_work)));
    CHECK(a.find("max_vertices")->dump() ==
          (e->max_vertices == 0 ? std::string("null") : std::to_string(e->max_vertices)));
  }
  std::set<std::string> names;
  for (const auto& a : algos) names.insert(a.find("name")->as_string());
  CHECK(names.count("paper-faithful") == 1);
  CHECK(names.count("paper-faithful-d1") == 1);
  CHECK(names.count("practical-local") == 1);

  const auto v = run({"--version"});
  CHECK(v.code == 0);
  CHECK(v.out.find("2611.34911v1") != std::string::npos);
  const auto h = run({"--help"});
  CHECK(h.code == 0);
  CHECK(h.out.find("--algorithms") != std::string::npos);
}
