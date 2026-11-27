// Argument parsing for sssp_bench_runner (docs/design/benchmark_runner_cli.md §2).
// Graph options use the names, meaning and defaults of tools/sssp_gen.cpp.

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "sssp/bench_runner.hpp"
#include "sssp/generators.hpp"

namespace sssp::bench {

namespace {

constexpr std::string_view kUsage =
    R"(usage: sssp_bench_runner --type TYPE --algorithms A[,B...] --env-json E --run-id ID [options]
  graph (same names as sssp_gen):
    --type TYPE  --n N  --avg-degree D  --degree D  --attach M  --rows R --cols C
    --exponent G  --multi keep|drop|switch  --connected
    --jitter J --keep-prob P --shortcut-prob P --shortcut-radius R --int-scale S   (road)
    --weights DIST|builtin  --lo X --hi X --mean X --zero-prob P --wide-log2 K
    --weight-type uint64|double  --seed S  --graph-file PATH (edge list, or DIMACS *.gr)
  measurement:
    --sources K (default 1) | --source-list a,b,...   --trials T (5)  --warmup W (1)
    --order-seed S (0)  --timeout-sec X  --verify none|reference|bellman-ford (reference)
    --param KEY=VALUE (repeatable; non-baseline algorithms)  --isolate
  output:
    --env-json PATH  --run-id ID  --out PATH (append JSONL; default stdout)
    --list  --version  --help
Exit codes: 0 ok, 1 error, 2 bad arguments, 3 unsupported combination, 5 out of memory.
)";

/// Families whose weights are fixed by the construction (they ignore --weights).
bool builtin_weight_family(std::string_view type) noexcept {
  return type == "road" || type == "adversarial_long_path" || type == "adversarial_decrease_key";
}

[[noreturn]] void usage_error(const std::string& message) {
  throw RunnerError(kExitUsage, message);
}

template <class T>
T parse_number(std::string_view flag, std::string_view text) {
  T value{};
  const char* first = text.data();  // [first, end) is explicitly sized; no NUL is needed
  const char* end = first + text.size();
  const auto [ptr, ec] = std::from_chars(first, end, value);
  if (text.empty() || ec != std::errc{} || ptr != end) {
    usage_error("invalid value for " + std::string(flag) + ": '" + std::string(text) + "'");
  }
  return value;
}

std::vector<std::string_view> split_commas(std::string_view flag, std::string_view text) {
  std::vector<std::string_view> parts;
  std::size_t start = 0;
  for (;;) {
    const std::size_t comma = text.find(',', start);
    const std::string_view part = text.substr(start, comma - start);
    if (part.empty())
      usage_error("empty item in " + std::string(flag) + " '" + std::string(text) + "'");
    parts.push_back(part);
    if (comma == std::string_view::npos) return parts;
    start = comma + 1;
  }
}

/// --param value: JSON when it parses (4 -> 4, true -> true), otherwise the text itself.
Json param_value(std::string_view raw) {
  try {
    return Json::parse(raw);
  } catch (const JsonError&) {
    return Json::string(std::string(raw));
  }
}

bool known_type(std::string_view type) {
  const auto catalog = gen::generator_catalog();
  return std::any_of(catalog.begin(), catalog.end(),
                     [&](const gen::GeneratorInfo& g) { return g.name == type; });
}

void validate(Options& o, bool type_given, bool weights_given, bool sources_given) {
  if (o.list || o.version || o.help) return;
  if (o.algorithms.empty()) usage_error("--algorithms is required");
  for (std::size_t i = 0; i < o.algorithms.size(); ++i) {
    for (std::size_t j = i + 1; j < o.algorithms.size(); ++j) {
      if (o.algorithms[i] == o.algorithms[j]) {
        usage_error("--algorithms lists '" + o.algorithms[i] + "' twice");
      }
    }
  }
  if (o.isolate && o.algorithms.size() != 1) {
    usage_error("--isolate runs exactly one algorithm per process; got " +
                std::to_string(o.algorithms.size()));
  }
  if (o.env_json.empty()) usage_error("--env-json is required (scripts/bench/envinfo.py)");
  if (o.run_id.empty()) usage_error("--run-id is required");
  if (o.trials == 0) usage_error("--trials must be at least 1");
  if (sources_given && !o.source_list.empty()) {
    usage_error("--sources and --source-list are mutually exclusive");
  }
  if (o.num_sources == 0) usage_error("--sources must be at least 1");
  if (o.timeout_sec && !(*o.timeout_sec > 0.0)) usage_error("--timeout-sec must be positive");

  if (o.graph_file.empty()) {
    if (!type_given) usage_error("--type (or --graph-file) is required");
    if (!known_type(o.spec.type)) usage_error("unknown --type '" + o.spec.type + "'");
  } else if (!type_given) {
    o.spec.type = "file";
  }
  const bool builtin_family = builtin_weight_family(o.spec.type);
  if (!weights_given) {
    o.weight_dist = builtin_family ? "builtin" : o.weight_dist;
  } else if (builtin_family && o.weight_dist != "builtin") {
    usage_error("--type " + o.spec.type + " has built-in weights; pass --weights builtin");
  } else if (!builtin_family && o.weight_dist == "builtin" && o.graph_file.empty()) {
    usage_error("--weights builtin is only valid for road / adversarial_* graphs");
  }
  if (o.weight_dist != "builtin") {
    try {
      gen::validate_weight_spec(o.spec.options.weights);
    } catch (const std::invalid_argument& e) {
      usage_error(std::string("invalid weight parameters: ") + e.what());
    }
  }

  // Baselines take no parameters; reject --param when nothing selected can use it.
  if (!o.params.items().empty()) {
    const bool any_taker = std::any_of(o.algorithms.begin(), o.algorithms.end(), [](const auto& a) {
      const AlgorithmEntry* e = find_algorithm(a);
      return e == nullptr || e->takes_params;  // unknown names are reported later (exit 3)
    });
    if (!any_taker) usage_error("--param given, but the selected algorithms take no parameters");
  }
}

}  // namespace

std::string_view usage_text() noexcept {
  return kUsage;
}

Options parse_args(std::span<const std::string> args) {
  Options o;
  auto& s = o.spec;
  auto& w = s.options.weights;
  bool type_given = false;
  bool weights_given = false;
  bool sources_given = false;
  for (std::size_t i = 0; i < args.size(); ++i) {
    const std::string_view a = args[i];
    auto value = [&]() -> std::string_view {
      if (i + 1 >= args.size()) usage_error("missing value for " + std::string(a));
      return args[++i];
    };
    if (a == "--help" || a == "-h") {
      o.help = true;
    } else if (a == "--list") {
      o.list = true;
    } else if (a == "--version") {
      o.version = true;
    } else if (a == "--type") {
      s.type = std::string(value());
      type_given = true;
    } else if (a == "--n") {
      s.n = parse_number<std::size_t>(a, value());
    } else if (a == "--avg-degree") {
      s.avg_degree = parse_number<double>(a, value());
      o.avg_degree_given = true;
    } else if (a == "--degree") {
      s.degree = parse_number<std::size_t>(a, value());
      o.degree_given = true;
    } else if (a == "--attach") {
      s.attach = parse_number<std::size_t>(a, value());
      o.attach_given = true;
    } else if (a == "--multi") {
      try {
        s.multi = gen::parse_multi_edge_policy(value());
      } catch (const std::invalid_argument& e) {
        usage_error(e.what());
      }
    } else if (a == "--rows") {
      s.rows = parse_number<std::size_t>(a, value());
    } else if (a == "--cols") {
      s.cols = parse_number<std::size_t>(a, value());
    } else if (a == "--exponent") {
      s.exponent = parse_number<double>(a, value());
    } else if (a == "--jitter") {
      s.road.jitter = parse_number<double>(a, value());
    } else if (a == "--keep-prob") {
      s.road.keep_probability = parse_number<double>(a, value());
    } else if (a == "--shortcut-prob") {
      s.road.shortcut_probability = parse_number<double>(a, value());
    } else if (a == "--shortcut-radius") {
      s.road.shortcut_radius = parse_number<int>(a, value());
    } else if (a == "--int-scale") {
      s.road.integer_scale = parse_number<double>(a, value());
    } else if (a == "--connected") {
      s.options.connected = true;
    } else if (a == "--weights") {
      const std::string_view d = value();
      weights_given = true;
      if (d == "builtin") {
        o.weight_dist = "builtin";
      } else {
        try {
          w.dist = gen::parse_weight_dist(d);
        } catch (const std::invalid_argument&) {
          usage_error("unknown --weights '" + std::string(d) + "'");
        }
        o.weight_dist = std::string(gen::to_string(w.dist));
      }
    } else if (a == "--lo") {
      w.lo = parse_number<double>(a, value());
    } else if (a == "--hi") {
      w.hi = parse_number<double>(a, value());
    } else if (a == "--mean") {
      w.mean = parse_number<double>(a, value());
    } else if (a == "--zero-prob") {
      w.zero_probability = parse_number<double>(a, value());
    } else if (a == "--wide-log2") {
      w.wide_log2 = parse_number<int>(a, value());
    } else if (a == "--weight-type") {
      const std::string_view t = value();
      if (t != "uint64" && t != "double") usage_error("--weight-type must be uint64|double");
      o.use_double = (t == "double");
    } else if (a == "--seed") {
      s.options.seed = parse_number<std::uint64_t>(a, value());
    } else if (a == "--graph-file") {
      o.graph_file = std::string(value());
      if (o.graph_file.empty()) usage_error("--graph-file needs a path");
    } else if (a == "--sources") {
      o.num_sources = parse_number<std::size_t>(a, value());
      sources_given = true;
    } else if (a == "--source-list") {  // a repeated list option replaces the earlier one
      o.source_list.clear();
      for (const auto part : split_commas(a, value())) {
        o.source_list.push_back(parse_number<Vertex>(a, part));
      }
    } else if (a == "--algorithms") {
      o.algorithms.clear();
      for (const auto part : split_commas(a, value())) o.algorithms.emplace_back(part);
    } else if (a == "--param") {
      const std::string_view kv = value();
      const std::size_t eq = kv.find('=');
      if (eq == std::string_view::npos || eq == 0) {
        usage_error("--param expects KEY=VALUE, got '" + std::string(kv) + "'");
      }
      o.params.set(kv.substr(0, eq), param_value(kv.substr(eq + 1)));
    } else if (a == "--trials") {
      o.trials = parse_number<std::size_t>(a, value());
    } else if (a == "--warmup") {
      o.warmup = parse_number<std::size_t>(a, value());
    } else if (a == "--order-seed") {
      o.order_seed = parse_number<std::uint64_t>(a, value());
    } else if (a == "--timeout-sec") {
      o.timeout_sec = parse_number<double>(a, value());
    } else if (a == "--verify") {
      const std::string_view v = value();
      if (v == "none") {
        o.verify = VerifyMode::kNone;
      } else if (v == "reference") {
        o.verify = VerifyMode::kReference;
      } else if (v == "bellman-ford") {
        o.verify = VerifyMode::kBellmanFord;
      } else {
        usage_error("--verify must be none|reference|bellman-ford");
      }
    } else if (a == "--isolate") {
      o.isolate = true;
    } else if (a == "--env-json") {
      o.env_json = std::string(value());
    } else if (a == "--run-id") {
      o.run_id = std::string(value());
    } else if (a == "--out") {
      o.out = std::string(value());
    } else {
      usage_error("unknown option: " + std::string(a));
    }
  }
  validate(o, type_given, weights_given, sources_given);
  return o;
}

void check_supported(const Options& o) {
  for (const auto& name : o.algorithms) {
    const AlgorithmEntry* e = find_algorithm(name);
    if (e == nullptr) {
      if (schema_mode(name)) {
        throw RunnerError(kExitUnsupported,
                          "algorithm '" + name + "' is not implemented yet in this build");
      }
      throw RunnerError(kExitUnsupported, "unknown algorithm '" + name + "' (see --list)");
    }
    if ((o.use_double ? e->solve_f64 == nullptr : e->solve_u64 == nullptr)) {
      throw RunnerError(kExitUnsupported, "algorithm '" + name +
                                              "' does not support --weight-type " +
                                              (o.use_double ? "double" : "uint64"));
    }
  }
}

std::optional<double> nominal_avg_degree(const Options& o) {
  const std::string& t = o.spec.type;
  // A generated graph uses the generator default when the option is omitted; a graph file is
  // described only by what the caller states explicitly.
  const bool generated = o.graph_file.empty();
  if ((t == "er" || t == "chung_lu") && (generated || o.avg_degree_given)) {
    return o.spec.avg_degree;
  }
  if (t == "regular" && (generated || o.degree_given)) return static_cast<double>(o.spec.degree);
  if (t == "ba" && (generated || o.attach_given)) return 2.0 * static_cast<double>(o.spec.attach);
  return std::nullopt;
}

}  // namespace sssp::bench
