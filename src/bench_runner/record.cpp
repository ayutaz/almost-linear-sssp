// Benchmark records (docs/design/benchmark_schema.md), --env-json loading and ru_maxrss.

#include <sys/resource.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <optional>
#include <span>
#include <string>
#include <string_view>

#include "sssp/bench_runner.hpp"
#include "sssp/version.hpp"

namespace sssp::bench {

// ------------------------------------------------------------------ memory

std::optional<std::uint64_t> maxrss_to_bytes(long ru_maxrss, bool unit_is_bytes) noexcept {
  if (ru_maxrss <= 0) return std::nullopt;
  const auto v = static_cast<std::uint64_t>(ru_maxrss);
  return unit_is_bytes ? v : v * 1024U;
}

std::optional<std::uint64_t> peak_rss_bytes() noexcept {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0) return std::nullopt;
#if defined(__APPLE__)
  constexpr bool kBytes = true;  // macOS: bytes
#else
  constexpr bool kBytes = false;  // Linux (and the BSDs): KiB
#endif
  return maxrss_to_bytes(usage.ru_maxrss, kBytes);
}

// ------------------------------------------------------------------ env

Json load_env(std::string_view text, bool diagnostics_compiled) {
  Json parsed;
  try {
    parsed = Json::parse(text);
  } catch (const JsonError& e) {
    throw RunnerError(kExitUsage, std::string("--env-json: ") + e.what());
  }
  if (parsed.kind() != Json::Kind::kObject) {
    throw RunnerError(kExitUsage, "--env-json: must contain a JSON object");
  }
  Json env = Json::object();
  for (const auto key : kEnvFields) {
    const Json* v = parsed.find(key);
    if (v == nullptr) {
      throw RunnerError(kExitUsage, "--env-json: missing field '" + std::string(key) +
                                        "' (write it with scripts/bench/envinfo.py)");
    }
    env.set(key, *v);
  }
  const Json& paper = *env.find("paper_version");
  if (paper.kind() != Json::Kind::kString || paper.as_string() != paper_arxiv_id()) {
    throw RunnerError(kExitUsage, "--env-json: paper_version " + paper.dump() +
                                      " does not match this build (" +
                                      std::string(paper_arxiv_id()) + ")");
  }
  const Json& diag = *env.find("diagnostics_enabled");
  if (diag.kind() != Json::Kind::kBool || diag.as_bool() != diagnostics_compiled) {
    throw RunnerError(kExitUsage, "--env-json: diagnostics_enabled " + diag.dump() +
                                      " does not match this build (SSSP_DIAGNOSTICS=" +
                                      (diagnostics_compiled ? "1" : "0") +
                                      "); was it captured for another preset?");
  }
  return env;
}

// ------------------------------------------------------------------ records

std::string utc_timestamp(std::chrono::system_clock::time_point t) {
  const std::time_t secs = std::chrono::system_clock::to_time_t(t);
  std::tm tm{};
  gmtime_r(&secs, &tm);
  std::array<char, 32> buf{};
  const std::size_t len = std::strftime(buf.data(), buf.size(), "%Y-%m-%dT%H:%M:%S+00:00", &tm);
  return {buf.data(), len};
}

namespace {

template <class T>
Json opt_int(const std::optional<T>& v) {
  return v ? Json::integer(static_cast<std::uint64_t>(*v)) : Json::null();
}

Json u64(std::uint64_t v) {
  return Json::integer(v);
}

Json ops_json(const OpCounters& o) {
  Json j = Json::object();
  j.set("vertex_scans", u64(o.vertex_scans));
  j.set("edge_scans", u64(o.edge_scans));
  j.set("relaxations", u64(o.relaxations));
  j.set("comparisons", u64(o.comparisons));
  j.set("additions", u64(o.additions));
  j.set("queue_pushes", u64(o.queue_pushes));
  j.set("queue_pops", u64(o.queue_pops));
  j.set("decrease_keys", u64(o.decrease_keys));
  return j;
}

Json fallback_json(const Diagnostics& d) {
  Json records = Json::array();
  std::uint64_t vertices = 0;
  std::uint64_t edges = 0;
  for (const auto& f : d.fallbacks) {
    Json r = Json::object();
    r.set("reason", Json::string(f.reason));
    r.set("vertices", u64(f.vertices));
    r.set("edges", u64(f.edges));
    records.push_back(std::move(r));
    vertices += f.vertices;
    edges += f.edges;
  }
  Json j = Json::object();
  j.set("count", u64(d.fallbacks.size()));
  j.set("records", std::move(records));
  j.set("vertices", u64(vertices));
  j.set("edges", u64(edges));
  return j;
}

}  // namespace

Json make_record(const Json& env, std::string_view run_id, std::string_view timestamp,
                 const GraphInfo& graph, const TrialRecord& t, bool diagnostics_compiled) {
  const bool verification_failed = t.correctness.verified && t.correctness.ok == false;
  Json r = Json::object();
  r.set("schema_version", Json::string("1"));
  r.set("run_id", Json::string(std::string(run_id)));
  r.set("timestamp", Json::string(std::string(timestamp)));
  for (const auto key : std::span(kEnvFields).first<3>()) r.set(key, *env.find(key));
  r.set("algorithm", Json::string(t.algorithm));
  r.set("mode", Json::string(t.mode));
  r.set("params", t.params);
  for (const auto key : std::span(kEnvFields).subspan<3>()) r.set(key, *env.find(key));
  r.set("graph_type", Json::string(graph.graph_type));
  r.set("n", u64(graph.n));
  r.set("m", u64(graph.m));
  r.set("avg_degree", graph.avg_degree ? Json::real(*graph.avg_degree) : Json::null());
  r.set("weight_dist", Json::string(graph.weight_dist));
  r.set("weight_type", Json::string(graph.weight_type));
  r.set("seed", u64(graph.seed));
  r.set("source", opt_int(t.source));
  r.set("reachable", opt_int(t.reachable));
  r.set("trial", opt_int(t.trial));
  r.set("warmup", Json::boolean(t.warmup));
  r.set("order_seed", opt_int(t.order_seed));
  r.set("isolated", Json::boolean(t.isolated));
  r.set("status", Json::string(verification_failed ? "error" : t.status));
  r.set("error", t.error ? Json::string(*t.error) : Json::null());
  Json times = Json::object();
  times.set("graph_gen", opt_int(t.times.graph_gen));
  times.set("preprocess", opt_int(t.times.preprocess));
  times.set("sssp", opt_int(t.times.sssp));
  times.set("tree", opt_int(t.times.tree));
  times.set("end_to_end", opt_int(t.times.end_to_end));
  times.set("verify", opt_int(t.times.verify));
  r.set("times_ns", std::move(times));
  // Peak RSS belongs to one algorithm only in an isolated process (schema §3).
  r.set("peak_rss_bytes", t.isolated ? opt_int(t.peak_rss) : Json::null());
  r.set("rss_after_graph_bytes", opt_int(t.rss_after_graph));
  // Counters are compiled out when SSSP_DIAGNOSTICS=0: null, never a fake 0.
  const bool counters = diagnostics_compiled && t.diag != nullptr;
  r.set("ops", counters ? ops_json(t.diag->ops) : Json::null());
  Json named = Json::object();
  if (counters) {
    for (const auto& [name, value] : t.diag->counters) named.set(name, u64(value));
  }
  r.set("counters", counters ? std::move(named) : Json::null());
  // Fallback records are kept in every build (CLAUDE.md fallback policy).
  r.set("fallback", t.diag != nullptr ? fallback_json(*t.diag) : Json::null());
  Json corr = Json::object();
  corr.set("verified", Json::boolean(t.correctness.verified));
  corr.set("ok", t.correctness.verified && t.correctness.ok ? Json::boolean(*t.correctness.ok)
                                                            : Json::null());
  corr.set("method", t.correctness.verified ? Json::string(t.correctness.method) : Json::null());
  r.set("correctness", std::move(corr));
  return r;
}

}  // namespace sssp::bench
