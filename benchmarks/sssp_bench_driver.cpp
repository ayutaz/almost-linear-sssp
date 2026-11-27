// Google Benchmark driver for the bench_runner registry (see sssp_bench_driver.hpp).

#include "sssp_bench_driver.hpp"

#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"

namespace sssp::bench::gbench {

namespace {

// ------------------------------------------------------------------ argument parsing

constexpr std::string_view kPrefix = "--sssp_";

std::vector<std::string> split_list(std::string_view key, std::string_view value) {
  std::vector<std::string> out;
  std::size_t start = 0;
  while (true) {
    const std::size_t comma = value.find(',', start);
    const std::string_view item = value.substr(
        start, comma == std::string_view::npos ? std::string_view::npos : comma - start);
    if (item.empty()) {
      throw std::invalid_argument("--sssp_" + std::string(key) + ": empty list item in '" +
                                  std::string(value) + "'");
    }
    out.emplace_back(item);
    if (comma == std::string_view::npos) break;
    start = comma + 1;
  }
  return out;
}

std::uint64_t parse_u64(std::string_view key, std::string_view text) {
  std::uint64_t v = 0;
  const auto [ptr, ec] = std::from_chars(text.data(), text.data() + text.size(), v);
  if (ec != std::errc{} || ptr != text.data() + text.size() || text.empty()) {
    throw std::invalid_argument("--sssp_" + std::string(key) + ": '" + std::string(text) +
                                "' is not a non-negative integer");
  }
  return v;
}

std::size_t parse_positive_size(std::string_view key, std::string_view text) {
  const std::uint64_t v = parse_u64(key, text);
  if (v == 0) throw std::invalid_argument("--sssp_" + std::string(key) + " must be at least 1");
  return static_cast<std::size_t>(v);
}

double parse_positive_double(std::string_view key, const std::string& text) {
  errno = 0;
  char* end = nullptr;
  const double v = std::strtod(text.c_str(), &end);
  if (text.empty() || end != text.c_str() + text.size() || errno != 0 || !std::isfinite(v) ||
      v <= 0.0) {
    throw std::invalid_argument("--sssp_" + std::string(key) + ": '" + text +
                                "' is not a positive finite number");
  }
  return v;
}

bool is_generator(std::string_view name) {
  const auto catalog = gen::generator_catalog();
  return std::any_of(catalog.begin(), catalog.end(),
                     [&](const gen::GeneratorInfo& g) { return g.name == name; });
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

std::string format_degree(double d) {
  std::ostringstream os;
  os << d;
  return os.str();
}

/// Nominal average degree written into the label: the value the generator was given (er,
/// chung_lu, geometric), the degree (regular), 2 * attach (ba); nullopt for fixed-degree families.
std::optional<double> nominal_degree(const gen::GraphSpec& spec) {
  if (spec.type == "er" || spec.type == "chung_lu" || spec.type == "geometric") {
    return spec.avg_degree;
  }
  if (spec.type == "regular") return static_cast<double>(spec.degree);
  if (spec.type == "ba") return 2.0 * static_cast<double>(spec.attach);
  return std::nullopt;
}

template <class W>
constexpr WeightType weight_type_v =
    std::is_same_v<W, double> ? WeightType::kF64 : WeightType::kU64;

// ------------------------------------------------------------------ the timed loop

template <class W>
void time_case(benchmark::State& state, const AlgorithmEntry& entry, const CsrGraph<W>& g,
               Vertex source) {
  using Clock = std::chrono::steady_clock;
  const SolveFn<W> solve = entry.solver<W>();
  const Json params = Json::object();
  for (auto _ : state) {
    PhaseTimer timer;
    const auto t0 = Clock::now();
    SsspResult<W> result = solve(g, source, params, timer);
    const auto t1 = Clock::now();
    benchmark::DoNotOptimize(result.dist.data());
    state.SetIterationTime(std::chrono::duration<double>(t1 - t0).count());
  }  // the result is destroyed here, outside the measured interval [t0, t1]
}

void run_case(benchmark::State& state, Suite& suite, std::size_t index) {
  const Suite::Case& c = suite.cases()[index];
  const CaseCheck& chk = suite.check(index);  // outside the timed loop
  if (chk.status == CaseCheck::Status::kRefused) {
    state.SkipWithMessage("refused: " + chk.message);
    return;
  }
  if (chk.status == CaseCheck::Status::kFailed) {
    state.SkipWithError("correctness check failed, not timed: " + chk.message);
    return;
  }
  std::size_t n = 0;
  std::size_t m = 0;
  if (c.weight_type == WeightType::kU64) {
    const auto& in = suite.inputs_u64()[c.input];
    n = in.graph.num_vertices();
    m = in.graph.num_edges();
    time_case<std::uint64_t>(state, *c.algorithm, in.graph, in.sources[c.source_index]);
  } else {
    const auto& in = suite.inputs_f64()[c.input];
    n = in.graph.num_vertices();
    m = in.graph.num_edges();
    time_case<double>(state, *c.algorithm, in.graph, in.sources[c.source_index]);
  }
  // Context for the report (from the checked run, not timed). Fallback share = fallback_vertices
  // / n and fallback_edges / m (CLAUDE.md: reports state it next to every timing).
  state.counters["n"] = benchmark::Counter(static_cast<double>(n));
  state.counters["m"] = benchmark::Counter(static_cast<double>(m));
  state.counters["fallback_records"] =
      benchmark::Counter(static_cast<double>(chk.fallback_records));
  state.counters["fallback_vertices"] =
      benchmark::Counter(static_cast<double>(chk.fallback_vertices));
  state.counters["fallback_edges"] = benchmark::Counter(static_cast<double>(chk.fallback_edges));
}

}  // namespace

// ------------------------------------------------------------------ configuration

std::string_view driver_usage() noexcept {
  return "sssp_bench [Google Benchmark flags] [driver options]\n"
         "\n"
         "Driver options (lists are comma separated):\n"
         "  --sssp_graphs=T,...      graph types (sssp_gen --type names; default er,grid)\n"
         "  --sssp_sizes=N,...       requested vertex counts (default 1000,10000)\n"
         "  --sssp_avg_degree=D      degree parameter of er/chung_lu/geometric/regular/ba "
         "(default 4)\n"
         "  --sssp_weights=W,...     weight distributions; the weight type follows from the\n"
         "                           distribution (default uniform_int,uniform_real). road and\n"
         "                           adversarial_long_path always use their builtin weights\n"
         "  --sssp_seed=S            graph and source-selection seed (default 1)\n"
         "  --sssp_sources=K         sources per graph (default 1)\n"
         "  --sssp_algorithms=A,...  registry names (default: every registered algorithm)\n"
         "  --sssp_verify_only       check every case against binary-heap Dijkstra and\n"
         "                           verify_sssp, time nothing (exit 1 on a mismatch)\n"
         "  --sssp_help              this text\n"
         "\n"
         "Only the solver call is timed (manual time). Graph generation, CSR construction, the\n"
         "reference run and verification happen before timing; a case that fails its check is\n"
         "reported as an error and never timed.\n";
}

std::string_view weight_type_of(std::string_view weight_dist) {
  static constexpr std::array<std::pair<std::string_view, std::string_view>, 8> kTypes{{
      {"constant", "uint64"},
      {"uniform_int", "uint64"},
      {"log_uniform", "uint64"},
      {"wide_range", "uint64"},
      {"zero_mixed", "uint64"},
      {"uniform_real", "double"},
      {"exponential", "double"},
      {"builtin", "uint64"},
  }};
  for (const auto& [name, type] : kTypes) {
    if (name == weight_dist) return type;
  }
  throw std::invalid_argument("unknown weight distribution '" + std::string(weight_dist) + "'");
}

bool has_builtin_weights(std::string_view graph_type) noexcept {
  return graph_type == "road" || graph_type == "adversarial_long_path";
}

DriverConfig parse_driver_args(std::span<const std::string> args) {
  DriverConfig c;
  for (const std::string& arg : args) {
    const std::string_view a = arg;
    if (a.substr(0, kPrefix.size()) != kPrefix) {
      throw std::invalid_argument("unknown argument '" + arg + "'");
    }
    const std::string_view body = a.substr(kPrefix.size());
    const std::size_t eq = body.find('=');
    const std::string_view key = body.substr(0, eq);
    const bool has_value = eq != std::string_view::npos;
    const std::string value = has_value ? std::string(body.substr(eq + 1)) : std::string();

    if (key == "verify_only" || key == "help") {
      if (has_value) throw std::invalid_argument("--sssp_" + std::string(key) + " takes no value");
      (key == "help" ? c.help : c.verify_only) = true;
      continue;
    }
    if (!has_value) {
      throw std::invalid_argument("--sssp_" + std::string(key) + " needs a value (=...)");
    }
    if (key == "graphs") {
      c.graph_types = split_list(key, value);
      for (const auto& t : c.graph_types) {
        if (!is_generator(t)) throw std::invalid_argument("unknown graph type '" + t + "'");
      }
    } else if (key == "sizes") {
      c.sizes.clear();
      for (const auto& s : split_list(key, value)) c.sizes.push_back(parse_positive_size(key, s));
    } else if (key == "avg_degree") {
      c.avg_degree = parse_positive_double(key, value);
    } else if (key == "weights") {
      c.weight_dists = split_list(key, value);
      for (const auto& w : c.weight_dists) {
        (void)weight_type_of(w);  // throws for unknown names
        if (w == "builtin") {
          throw std::invalid_argument(
              "--sssp_weights: 'builtin' is implied for road and adversarial_long_path and is not "
              "a distribution of the other families");
        }
      }
    } else if (key == "seed") {
      c.seed = parse_u64(key, value);
    } else if (key == "sources") {
      c.num_sources = parse_positive_size(key, value);
    } else if (key == "algorithms") {
      c.algorithms = split_list(key, value);
      for (const auto& name : c.algorithms) {
        if (find_algorithm(name) == nullptr) {
          throw std::invalid_argument("'" + name +
                                      "' is not a registered algorithm (sssp_bench_runner --list)");
        }
      }
    } else {
      throw std::invalid_argument("unknown option '--sssp_" + std::string(key) + "'");
    }
  }
  return c;
}

// ------------------------------------------------------------------ checking one case

std::optional<std::string> refusal_reason(const AlgorithmEntry& entry, std::size_t n,
                                          std::size_t m) {
  if (entry.max_vertices != 0 && n > entry.max_vertices) {
    return std::string(entry.name) + " is limited to n <= " + std::to_string(entry.max_vertices) +
           " (n = " + std::to_string(n) + ")";
  }
  const std::uint64_t work = static_cast<std::uint64_t>(n) * static_cast<std::uint64_t>(m);
  if (entry.max_work != 0 && work > entry.max_work) {
    return std::string(entry.name) + " is limited to n * m <= " + std::to_string(entry.max_work) +
           " (n * m = " + std::to_string(work) + ")";
  }
  return std::nullopt;
}

template <class W>
CaseCheck check_case(const AlgorithmEntry& entry, const CsrGraph<W>& g, Vertex source,
                     const std::vector<W>& reference) {
  CaseCheck out;
  const SolveFn<W> solve = entry.solver<W>();
  if (solve == nullptr) {
    out.message = std::string(entry.name) + " has no solver for this weight type";
    return out;
  }
  std::optional<SsspResult<W>> result;
  try {
    PhaseTimer timer;
    result = solve(g, source, Json::object(), timer);
  } catch (const std::overflow_error& e) {
    out.message = std::string("integer distance overflow: ") + e.what();
    return out;
  } catch (const std::exception& e) {
    out.message = std::string("solver threw: ") + e.what();
    return out;
  }
  for (const FallbackRecord& f : result->diag.fallbacks) {
    ++out.fallback_records;
    out.fallback_vertices += f.vertices;
    out.fallback_edges += f.edges;
  }
  // Same rule as sssp_bench_runner: relaxation baselines agree bit for bit even for double
  // (numerics.md §2); other algorithms get the two-rounded-results bound (numerics.md §3.2).
  const std::size_t n = g.num_vertices();
  const Tolerance tol = (!std::is_same_v<W, double> || entry.exact_relaxation)
                            ? Tolerance::exact()
                            : Tolerance::for_path_hops(n == 0 ? 0 : 2 * (n - 1));
  const VerificationReport rep = verify_sssp(g, source, *result, &reference, tol);
  if (!rep.ok) {
    out.message = "verification failed: " + first_lines(rep.summary(), 400);
    return out;
  }
  out.status = CaseCheck::Status::kOk;
  out.message = "matches dijkstra-binary; verify_sssp ok";
  return out;
}

template CaseCheck check_case<std::uint64_t>(const AlgorithmEntry&, const CsrGraph<std::uint64_t>&,
                                             Vertex, const std::vector<std::uint64_t>&);
template CaseCheck check_case<double>(const AlgorithmEntry&, const CsrGraph<double>&, Vertex,
                                      const std::vector<double>&);

// ------------------------------------------------------------------ the suite

template <>
std::vector<GraphInput<std::uint64_t>>& Suite::inputs<std::uint64_t>() noexcept {
  return u64_;
}
template <>
std::vector<GraphInput<double>>& Suite::inputs<double>() noexcept {
  return f64_;
}

Suite::Suite(const DriverConfig& config) {
  if (config.algorithms.empty()) {
    for (const AlgorithmEntry& e : algorithm_registry()) algorithms_.push_back(&e);
  } else {
    for (const auto& name : config.algorithms) {
      const AlgorithmEntry* e = find_algorithm(name);
      if (e == nullptr) throw std::invalid_argument("'" + name + "' is not a registered algorithm");
      algorithms_.push_back(e);
    }
  }
  for (const auto& type : config.graph_types) {
    for (const std::size_t n : config.sizes) {
      const std::vector<std::string> dists =
          has_builtin_weights(type) ? std::vector<std::string>{"builtin"} : config.weight_dists;
      for (const auto& dist : dists) {
        if (weight_type_of(dist) == "double") {
          add_input<double>(config, type, n, dist);
        } else {
          add_input<std::uint64_t>(config, type, n, dist);
        }
      }
    }
  }
  checks_.resize(cases_.size());
}

template <class W>
void Suite::add_input(const DriverConfig& config, const std::string& type, std::size_t n,
                      const std::string& dist) {
  gen::GraphSpec spec;
  spec.type = type;
  spec.n = n;
  spec.avg_degree = config.avg_degree;
  spec.degree = static_cast<std::size_t>(std::llround(config.avg_degree));
  spec.attach =
      std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(config.avg_degree / 2)));
  spec.options.seed = config.seed;
  if (dist != "builtin") spec.options.weights.dist = gen::parse_weight_dist(dist);

  EdgeList<W> el;
  try {
    el = gen::generate<W>(spec);
  } catch (const std::invalid_argument& e) {
    throw std::invalid_argument("cannot generate " + type + " with n = " + std::to_string(n) +
                                ": " + e.what());
  }
  GraphInput<W> in;
  in.graph_type = type;
  in.weight_dist = dist;
  in.graph = CsrGraph<W>::from_edges(el.num_vertices, std::move(el.edges));
  const std::size_t actual_n = in.graph.num_vertices();
  if (actual_n < config.num_sources) {
    throw std::invalid_argument(type + " with n = " + std::to_string(n) + " has " +
                                std::to_string(actual_n) + " vertices, fewer than --sssp_sources " +
                                std::to_string(config.num_sources));
  }
  in.sources = select_sources(actual_n, config.num_sources, config.seed);
  in.reference.resize(in.sources.size());
  in.reference_error.resize(in.sources.size());
  const auto degree = nominal_degree(spec);
  in.label = type + "/n:" + std::to_string(actual_n) +
             (degree ? "/d:" + format_degree(*degree) : std::string()) + "/" + dist + "/" +
             std::string(weight_type_of(dist));

  auto& list = inputs<W>();
  const std::size_t index = list.size();
  list.push_back(std::move(in));
  const GraphInput<W>& stored = list.back();
  for (std::size_t si = 0; si < stored.sources.size(); ++si) {
    for (const AlgorithmEntry* e : algorithms_) {
      if (e->solver<W>() == nullptr) continue;  // e.g. radix heap x double: not a case
      cases_.push_back(
          Case{e, weight_type_v<W>, index, si,
               std::string(e->name) + "/" + stored.label + "/src:" + std::to_string(si)});
    }
  }
}

template <class W>
CaseCheck Suite::check_impl(const Case& c) {
  GraphInput<W>& in = inputs<W>()[c.input];
  if (auto reason = refusal_reason(*c.algorithm, in.graph.num_vertices(), in.graph.num_edges())) {
    CaseCheck r;
    r.status = CaseCheck::Status::kRefused;
    r.message = std::move(*reason);
    return r;
  }
  const std::size_t si = c.source_index;
  const Vertex s = in.sources[si];
  if (!in.reference[si] && !in.reference_error[si]) {
    try {
      in.reference[si] = baseline::dijkstra_binary_heap(in.graph, s).dist;
    } catch (const std::overflow_error& e) {
      in.reference_error[si] =
          std::string("reference failed: integer distance overflow: ") + e.what();
    }
  }
  if (in.reference_error[si]) {
    CaseCheck r;
    r.message = *in.reference_error[si];
    return r;
  }
  return check_case<W>(*c.algorithm, in.graph, s, *in.reference[si]);
}

const CaseCheck& Suite::check(std::size_t i) {
  if (!checks_.at(i)) {
    const Case& c = cases_[i];
    checks_[i] =
        c.weight_type == WeightType::kU64 ? check_impl<std::uint64_t>(c) : check_impl<double>(c);
  }
  return *checks_[i];
}

std::vector<std::string> Suite::register_benchmarks() {
  std::vector<std::string> names;
  names.reserve(cases_.size());
  for (std::size_t i = 0; i < cases_.size(); ++i) {
    benchmark::RegisterBenchmark(cases_[i].name, [this, i](benchmark::State& state) {
      run_case(state, *this, i);
    })->UseManualTime();
    names.push_back(cases_[i].name);
  }
  return names;
}

bool Suite::verify_all(std::ostream& out) {
  std::size_t ok = 0;
  std::size_t refused = 0;
  std::size_t failed = 0;
  for (std::size_t i = 0; i < cases_.size(); ++i) {
    const CaseCheck& chk = check(i);
    out << cases_[i].name << ": ";
    switch (chk.status) {
      case CaseCheck::Status::kOk:
        ++ok;
        out << "ok";
        break;
      case CaseCheck::Status::kRefused:
        ++refused;
        out << "refused: " << chk.message;
        break;
      case CaseCheck::Status::kFailed:
        ++failed;
        out << "FAILED: " << chk.message;
        break;
    }
    out << '\n';
  }
  out << cases_.size() << " cases: " << ok << " ok, " << refused << " refused, " << failed
      << " failed\n";
  return failed == 0;
}

bool Suite::any_failed() const noexcept {
  return std::any_of(checks_.begin(), checks_.end(), [](const std::optional<CaseCheck>& c) {
    return c && c->status == CaseCheck::Status::kFailed;
  });
}

}  // namespace sssp::bench::gbench
