// Driver of sssp_bench_runner (docs/design/benchmark_runner_cli.md §3):
// build the graph once, then for each source run W warm-ups and T trials of every algorithm in a
// seeded per-round order, verify outside the timed region, and write one record per run.

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iterator>
#include <new>
#include <optional>
#include <ostream>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/bench_runner.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/graph_io.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"
#include "sssp/version.hpp"

namespace sssp::bench {

namespace {

using Clock = std::chrono::steady_clock;

std::string read_text_file(const std::string& path, std::string_view what) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw RunnerError(kExitUsage, std::string(what) + ": cannot read '" + path + "'");
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

bool ends_with(std::string_view s, std::string_view suffix) {
  return s.size() >= suffix.size() && s.substr(s.size() - suffix.size()) == suffix;
}

/// Writes records as JSONL. Non-isolated runs append and flush line by line; isolated runs keep
/// the records and write them just before exit with the final ru_maxrss (CLI §3.6).
class RecordSink {
 public:
  RecordSink(std::string path, std::ostream& stdout_stream, bool isolate)
      : path_(std::move(path)), stdout_(stdout_stream), isolate_(isolate) {}

  void add(Json record) {
    ++count_;
    if (isolate_) {
      pending_.push_back(std::move(record));
    } else {
      write(record);
    }
  }

  void finish() {
    if (!isolate_) return;
    const auto peak = peak_rss_bytes();
    for (auto& r : pending_) {
      r.set("peak_rss_bytes", peak ? Json::integer(*peak) : Json::null());
      write(r);
    }
    pending_.clear();
  }

  [[nodiscard]] std::size_t count() const noexcept { return count_; }

 private:
  std::ostream& stream() {
    if (path_.empty()) return stdout_;
    if (!file_.is_open()) {
      file_.open(path_, std::ios::out | std::ios::app);
      if (!file_) throw RunnerError(kExitError, "cannot open --out '" + path_ + "' for appending");
    }
    return file_;
  }

  void write(const Json& record) {
    std::string line = record.dump();
    line += '\n';
    std::ostream& os = stream();
    os << line;
    os.flush();
    if (!os) throw RunnerError(kExitError, "failed to write a record to '" + path_ + "'");
  }

  std::string path_;
  std::ostream& stdout_;
  bool isolate_;
  std::ofstream file_;
  std::vector<Json> pending_;
  std::size_t count_ = 0;
};

template <class W>
EdgeList<W> load_graph(const Options& o) {
  if (o.graph_file.empty()) {
    try {
      return gen::generate<W>(o.spec);
    } catch (const std::invalid_argument& e) {
      throw RunnerError(kExitUsage, std::string("graph generation: ") + e.what());
    }
  }
  try {
    if (ends_with(o.graph_file, ".gr")) {
      // sssp_gen --format dimacs writes every edge as two arcs; merging inverts it exactly.
      return io::read_dimacs_file<W>(
          o.graph_file, {io::DimacsArcMode::kMergeReversePairs, io::DuplicateEdgePolicy::kKeepAll});
    }
    return io::read_edge_list_file<W>(o.graph_file);
  } catch (const std::runtime_error& e) {  // cannot open, or io::ParseError
    throw RunnerError(kExitUsage, "--graph-file '" + o.graph_file + "': " + e.what());
  }
}

std::string first_lines(const std::string& text, std::size_t max_chars) {
  std::string flat;
  for (const char c : text) {
    if (c == '\n') {
      if (!flat.empty() && flat.back() != ' ') flat += "; ";
    } else {
      flat += c;
    }
    if (flat.size() >= max_chars) {
      flat += "...";
      break;
    }
  }
  while (!flat.empty() && (flat.back() == ' ' || flat.back() == ';')) flat.pop_back();
  return flat;
}

/// Reference distances for one source, computed once (outside every timed region).
template <class W>
struct Reference {
  std::optional<std::vector<W>> dist;
  std::optional<std::string> error;
};

template <class W>
int run_benchmark(const Options& o, const Json& env, Clock::time_point start, std::ostream& out,
                  std::ostream& err) {
  constexpr bool kDouble = std::is_same_v<W, double>;

  // 1. Graph: generated (or read) once and shared by every algorithm and trial.
  const auto g0 = Clock::now();
  EdgeList<W> el = load_graph<W>(o);
  const CsrGraph<W> g = CsrGraph<W>::from_edges(el.num_vertices, std::move(el.edges));
  const std::uint64_t graph_gen_ns = PhaseTimer::to_ns(Clock::now() - g0);
  const auto rss_after_graph = peak_rss_bytes();

  const std::size_t n = g.num_vertices();
  const std::size_t m = g.num_edges();
  if (n == 0) throw RunnerError(kExitUsage, "the graph has no vertices");
  const GraphInfo info{
      o.spec.type,        n, m, nominal_avg_degree(o), o.weight_dist, kDouble ? "double" : "uint64",
      o.spec.options.seed};

  // 2. Everything that can be refused is refused before the first record is written.
  std::vector<const AlgorithmEntry*> algos;
  const std::uint64_t work = static_cast<std::uint64_t>(n) * static_cast<std::uint64_t>(m);
  for (const auto& name : o.algorithms) {
    const AlgorithmEntry* e = find_algorithm(name);
    if (e->max_work != 0 && work > e->max_work) {
      throw RunnerError(kExitUnsupported,
                        name + " is limited to n * m <= " + std::to_string(e->max_work) +
                            "; this graph has " + std::to_string(work));
    }
    if (e->max_vertices != 0 && n > e->max_vertices) {
      throw RunnerError(kExitUnsupported,
                        name + " is limited to n <= " + std::to_string(e->max_vertices) +
                            "; this graph has n = " + std::to_string(n));
    }
    algos.push_back(e);
  }
  if (o.verify == VerifyMode::kBellmanFord && work > kBellmanFordMaxWork) {
    throw RunnerError(kExitUnsupported, "--verify bellman-ford is limited to n * m <= " +
                                            std::to_string(kBellmanFordMaxWork) +
                                            "; use --verify reference");
  }
  std::vector<Vertex> sources = o.source_list;
  if (sources.empty()) {
    if (o.num_sources > n) {
      throw RunnerError(kExitUsage, "--sources " + std::to_string(o.num_sources) +
                                        " exceeds the number of vertices " + std::to_string(n));
    }
    sources = select_sources(n, o.num_sources, o.spec.options.seed);
  }
  for (const Vertex s : sources) {
    if (s >= n) {
      throw RunnerError(kExitUsage, "source " + std::to_string(s) +
                                        " is not a vertex (n = " + std::to_string(n) + ")");
    }
  }

  const std::string method =
      o.verify == VerifyMode::kBellmanFord ? "bellman-ford" : "reference:dijkstra-binary";
  RecordSink sink(o.out, out, o.isolate);
  const std::size_t rounds = o.warmup + o.trials;
  std::size_t planned = sources.size() * rounds * algos.size();
  bool timed_out = false;

  try {
    for (std::size_t si = 0; si < sources.size() && !timed_out; ++si) {
      const Vertex s = sources[si];
      Reference<W> ref;
      for (std::size_t round = 0; round < rounds && !timed_out; ++round) {
        const auto order = round_order(algos.size(), o.order_seed, si, round);
        for (const std::size_t idx : order) {
          if (o.timeout_sec &&
              std::chrono::duration<double>(Clock::now() - start).count() >= *o.timeout_sec) {
            timed_out = true;
            break;
          }
          const AlgorithmEntry& e = *algos[idx];
          TrialRecord t;
          t.algorithm = std::string(e.name);
          t.mode = std::string(e.mode);
          if (e.takes_params) t.params = o.params;
          t.source = s;
          t.warmup = round < o.warmup;
          t.trial = t.warmup ? round : round - o.warmup;
          t.order_seed = o.order_seed;
          t.isolated = o.isolate;
          t.times.graph_gen = graph_gen_ns;
          t.rss_after_graph = rss_after_graph;

          // 3. Timed region: preprocess + sssp + tree, nothing else.
          PhaseTimer timer;
          std::optional<SsspResult<W>> result;
          try {
            const auto e0 = Clock::now();
            SsspResult<W> r = e.template solver<W>()(g, s, o.params, timer);
            const auto e1 = Clock::now();
            t.times.end_to_end = PhaseTimer::to_ns(e1 - e0);
            result = std::move(r);
          } catch (const std::overflow_error& ex) {
            t.status = "error";
            t.error = std::string("integer distance overflow: ") + ex.what();
          }

          // 4. Outside the timed region: counts, verification, record.
          if (result) {
            t.times.preprocess = timer.ns(PhaseTimer::Phase::kPreprocess);
            t.times.sssp = timer.ns(PhaseTimer::Phase::kSssp);
            t.times.tree = timer.ns(PhaseTimer::Phase::kTree);
            t.reachable = static_cast<std::uint64_t>(
                std::count_if(result->dist.begin(), result->dist.end(),
                              [](W d) { return d != WeightTraits<W>::infinity(); }));
            t.diag = &result->diag;
            if (o.verify != VerifyMode::kNone) {
              if (!ref.dist && !ref.error) {
                try {
                  ref.dist = o.verify == VerifyMode::kBellmanFord
                                 ? baseline::bellman_ford(g, s).dist
                                 : baseline::dijkstra_binary_heap(g, s).dist;
                } catch (const std::overflow_error& ex) {
                  ref.error =
                      std::string("reference failed: integer distance overflow: ") + ex.what();
                }
              }
              if (ref.dist) {
                // Relaxation baselines agree bit for bit even for double (numerics.md §2);
                // other algorithms get the two-rounded-results bound (numerics.md §3.2).
                const Tolerance tol = (!kDouble || e.exact_relaxation)
                                          ? Tolerance::exact()
                                          : Tolerance::for_path_hops(2 * (n - 1));
                const auto v0 = Clock::now();
                const VerificationReport rep = verify_sssp(g, s, *result, &*ref.dist, tol);
                t.times.verify = PhaseTimer::to_ns(Clock::now() - v0);
                t.correctness = {true, rep.ok, method};
                if (!rep.ok) {
                  t.status = "error";
                  t.error = "verification failed: " + first_lines(rep.summary(), 400);
                }
              } else {
                t.status = "error";
                t.error = ref.error;
              }
            }
          }
          sink.add(make_record(env, o.run_id, utc_timestamp(std::chrono::system_clock::now()), info,
                               t, kDiagnosticsCompiled));
        }
      }
    }
  } catch (...) {
    try {
      sink.finish();  // keep what was measured before the failure
    } catch (...) {   // NOLINT(bugprone-empty-catch): the original exception is the one to report
    }
    throw;
  }
  sink.finish();
  if (timed_out) {
    err << "sssp_bench_runner: --timeout-sec " << *o.timeout_sec << " reached; wrote "
        << sink.count() << " of " << planned << " planned record(s) and started no new run\n";
  }
  return kExitOk;
}

}  // namespace

Json list_json() {
  Json types = Json::array();
  for (const auto& g : gen::generator_catalog()) types.push_back(Json::string(std::string(g.name)));
  Json dists = Json::array();
  for (const auto d :
       {gen::WeightDist::kUniformInt, gen::WeightDist::kUniformReal, gen::WeightDist::kLogUniform,
        gen::WeightDist::kExponential, gen::WeightDist::kConstant, gen::WeightDist::kZeroMixed,
        gen::WeightDist::kWideRange}) {
    dists.push_back(Json::string(std::string(gen::to_string(d))));
  }
  dists.push_back(Json::string("builtin"));
  Json algos = Json::array();
  for (const auto& e : algorithm_registry()) {
    Json a = Json::object();
    a.set("name", Json::string(std::string(e.name)));
    a.set("mode", Json::string(std::string(e.mode)));
    Json wt = Json::array();
    if (e.solve_u64 != nullptr) wt.push_back(Json::string("uint64"));
    if (e.solve_f64 != nullptr) wt.push_back(Json::string("double"));
    a.set("weight_types", std::move(wt));
    a.set("takes_params", Json::boolean(e.takes_params));
    a.set("max_work", e.max_work == 0 ? Json::null() : Json::integer(e.max_work));
    a.set("max_vertices", e.max_vertices == 0 ? Json::null() : Json::integer(e.max_vertices));
    algos.push_back(std::move(a));
  }
  Json j = Json::object();
  j.set("graph_types", std::move(types));
  j.set("weight_dists", std::move(dists));
  Json weight_types = Json::array();
  weight_types.push_back(Json::string("uint64"));
  weight_types.push_back(Json::string("double"));
  j.set("weight_types", std::move(weight_types));
  j.set("algorithms", std::move(algos));
  return j;
}

int run_main(std::span<const std::string> args, std::ostream& out, std::ostream& err) {
  const auto start = Clock::now();
  try {
    const Options o = parse_args(args);
    if (o.help) {
      out << usage_text();
      return kExitOk;
    }
    if (o.version) {
      out << "sssp_bench_runner " << version_string() << " (paper " << paper_arxiv_id() << ")\n";
      return kExitOk;
    }
    if (o.list) {
      out << list_json().dump() << "\n";
      return kExitOk;
    }
    check_supported(o);
    const Json env = load_env(read_text_file(o.env_json, "--env-json"), kDiagnosticsCompiled);
    return o.use_double ? run_benchmark<double>(o, env, start, out, err)
                        : run_benchmark<std::uint64_t>(o, env, start, out, err);
  } catch (const RunnerError& e) {
    err << "sssp_bench_runner: " << e.what() << "\n";
    if (e.code() == kExitUsage) err << "(see sssp_bench_runner --help)\n";
    return e.code();
  } catch (const std::bad_alloc&) {
    err << "sssp_bench_runner: out of memory (std::bad_alloc)\n";
    return kExitOutOfMemory;
  } catch (const std::exception& e) {
    err << "sssp_bench_runner: " << e.what() << "\n";
    return kExitError;
  }
}

}  // namespace sssp::bench
