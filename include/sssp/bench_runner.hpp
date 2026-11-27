#pragma once
// Library behind tools/sssp_bench_runner (one isolated benchmark job per process).
//
//   CLI:      docs/design/benchmark_runner_cli.md
//   Records:  docs/design/benchmark_schema.md (JSONL, schema_version "1")
//   Protocol: docs/design/benchmark_protocol.md
//
// The pieces are separate functions so that argument parsing, source / order selection and
// record construction are unit-testable without running a benchmark.

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::bench {

// ------------------------------------------------------------------ exit codes (CLI §4)

inline constexpr int kExitOk = 0;
inline constexpr int kExitError = 1;
inline constexpr int kExitUsage = 2;
inline constexpr int kExitUnsupported = 3;
inline constexpr int kExitOutOfMemory = 5;

/// Error that ends the runner with a specific exit code.
class RunnerError : public std::runtime_error {
 public:
  RunnerError(int code, const std::string& message) : std::runtime_error(message), code_(code) {}
  [[nodiscard]] int code() const noexcept { return code_; }

 private:
  int code_;
};

/// True when the library was compiled with SSSP_DIAGNOSTICS=1 (operation counters present).
inline constexpr bool kDiagnosticsCompiled = SSSP_DIAGNOSTICS != 0;

// ------------------------------------------------------------------ JSON

class JsonError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

/// Small JSON value. Objects keep insertion order; numbers keep their source token, so values
/// copied from --env-json are written back unchanged (1.0 stays 1.0, integers stay integers).
class Json {
 public:
  enum class Kind : std::uint8_t { kNull, kBool, kNumber, kString, kArray, kObject };
  using Array = std::vector<Json>;
  using Object = std::vector<std::pair<std::string, Json>>;

  Json() = default;
  static Json null() { return {}; }
  static Json boolean(bool b);
  static Json integer(std::uint64_t v);
  static Json integer(std::int64_t v);
  /// Shortest round-trip representation; throws JsonError for NaN / infinity.
  static Json real(double v);
  /// Number from its JSON token (kept verbatim); throws JsonError if it is not a JSON number.
  static Json raw_number(std::string_view token);
  static Json string(std::string s);
  static Json array();
  static Json object();
  /// Strict RFC 8259 parser (one value, optional surrounding whitespace); throws JsonError.
  static Json parse(std::string_view text);

  [[nodiscard]] Kind kind() const noexcept { return kind_; }
  [[nodiscard]] bool is_null() const noexcept { return kind_ == Kind::kNull; }
  /// Number written without fraction or exponent.
  [[nodiscard]] bool is_integer() const noexcept;
  [[nodiscard]] bool as_bool() const;
  [[nodiscard]] const std::string& as_string() const;
  [[nodiscard]] const std::string& number_text() const;
  [[nodiscard]] const Array& elements() const;
  [[nodiscard]] const Object& items() const;

  /// Object member or nullptr.
  [[nodiscard]] const Json* find(std::string_view key) const;
  [[nodiscard]] Json* find_mut(std::string_view key);
  /// Sets an object member (replaces an existing key in place, else appends).
  void set(std::string_view key, Json value);
  void push_back(Json value);

  /// Compact single-line serialisation.
  [[nodiscard]] std::string dump() const;
  void dump_to(std::string& out) const;

 private:
  Kind kind_ = Kind::kNull;
  bool bool_ = false;
  std::string text_;  ///< string value or number token
  Array array_;
  Object object_;
};

// ------------------------------------------------------------------ options (CLI §2)

enum class VerifyMode : std::uint8_t { kNone, kReference, kBellmanFord };

struct Options {
  gen::GraphSpec spec{};  ///< graph options, same names and defaults as sssp_gen
  bool use_double = false;
  std::string weight_dist = "uniform_int";  ///< record name; "builtin" for fixed-weight families
  bool avg_degree_given = false;
  bool degree_given = false;
  bool attach_given = false;
  std::string graph_file;  ///< --graph-file (edge list, or DIMACS when it ends in .gr)

  std::size_t num_sources = 1;
  std::vector<Vertex> source_list;
  std::vector<std::string> algorithms;
  Json params = Json::object();  ///< --param KEY=VALUE (non-baseline algorithms only)
  std::size_t trials = 5;
  std::size_t warmup = 1;
  std::uint64_t order_seed = 0;
  std::optional<double> timeout_sec;
  VerifyMode verify = VerifyMode::kReference;
  bool isolate = false;

  std::string env_json;
  std::string run_id;
  std::string out;  ///< empty: records go to standard output

  bool list = false;
  bool version = false;
  bool help = false;
};

/// Parses the arguments (without argv[0]). Throws RunnerError(kExitUsage) for any argument error,
/// including missing required options and inconsistent combinations.
Options parse_args(std::span<const std::string> args);

/// Throws RunnerError(kExitUnsupported) for algorithms that are unknown or not implemented yet and
/// for algorithm x weight-type combinations an algorithm cannot run (radix heap x double).
void check_supported(const Options& options);

/// Nominal average degree written to records: --avg-degree (er, chung_lu), --degree (regular),
/// 2 * --attach (ba); nullopt for families whose degree is fixed by construction.
std::optional<double> nominal_avg_degree(const Options& options);

/// Usage text printed by --help.
std::string_view usage_text() noexcept;

// ------------------------------------------------------------------ algorithm registry

/// Phase timer handed to an algorithm adapter. Each adapter wraps exactly the region it can
/// separate (preprocess / sssp / tree); phases it does not time stay unmeasured (null).
class PhaseTimer {
 public:
  using Clock = std::chrono::steady_clock;
  enum class Phase : std::uint8_t { kPreprocess = 0, kSssp = 1, kTree = 2 };

  template <class F>
  decltype(auto) time(Phase p, F&& f) {
    const auto t0 = Clock::now();
    if constexpr (std::is_void_v<decltype(f())>) {
      f();
      add(p, Clock::now() - t0);
    } else {
      auto result = f();
      add(p, Clock::now() - t0);
      return result;
    }
  }
  [[nodiscard]] std::optional<std::uint64_t> ns(Phase p) const noexcept {
    return ns_[static_cast<std::size_t>(p)];
  }
  void reset() noexcept { ns_ = {}; }

  /// Clock reading in ns, at least 1: a region shorter than one clock tick is reported as 1 ns
  /// (the schema reserves 0 / null for "not measured").
  static std::uint64_t to_ns(Clock::duration d) noexcept;

 private:
  void add(Phase p, Clock::duration d) noexcept;
  std::array<std::optional<std::uint64_t>, 3> ns_{};
};

template <class W>
using SolveFn = SsspResult<W> (*)(const CsrGraph<W>& g, Vertex source, const Json& params,
                                  PhaseTimer& timer);

struct AlgorithmEntry {
  std::string_view name;  ///< schema name (bench.schema.ALGORITHM_MODES)
  std::string_view mode;  ///< baseline / paper-faithful / paper-experimental / practical-optimized
  SolveFn<std::uint64_t> solve_u64 = nullptr;  ///< nullptr: weight type not supported
  SolveFn<double> solve_f64 = nullptr;
  /// Relaxation algorithm: double distances are bit-identical to Dijkstra's
  /// (docs/design/numerics.md §2), so the reference comparison is exact. Otherwise the
  /// comparison uses Tolerance::for_path_hops(2(n-1)) (numerics.md §3.2).
  bool exact_relaxation = true;
  bool takes_params = false;   ///< accepts --param
  std::uint64_t max_work = 0;  ///< limit on n * m (0: none); larger graphs exit with code 3

  template <class W>
  [[nodiscard]] SolveFn<W> solver() const noexcept {
    if constexpr (std::is_same_v<W, double>) {
      return solve_f64;
    } else {
      return solve_u64;
    }
  }
};

/// Bellman-Ford is an O(nm) oracle: refused (exit 3) when n * m exceeds this.
inline constexpr std::uint64_t kBellmanFordMaxWork = std::uint64_t{1} << 30;

/// Every algorithm the runner can execute. New implementations are added here only.
std::span<const AlgorithmEntry> algorithm_registry() noexcept;
const AlgorithmEntry* find_algorithm(std::string_view name) noexcept;
/// Mode of a name known to the schema (including ones not implemented yet), else nullopt.
std::optional<std::string_view> schema_mode(std::string_view name) noexcept;

// ------------------------------------------------------------------ sources and order

/// k distinct sources in [0, n), a deterministic function of (n, seed); a prefix of the k + 1
/// selection. Throws std::invalid_argument if k > n.
std::vector<Vertex> select_sources(std::size_t n, std::size_t k, std::uint64_t seed);

/// Permutation of [0, count) used as the algorithm order of one round (one warm-up or trial)
/// for one source, a deterministic function of (order_seed, source_index, round).
std::vector<std::size_t> round_order(std::size_t count, std::uint64_t order_seed,
                                     std::uint64_t source_index, std::uint64_t round);

// ------------------------------------------------------------------ memory

/// ru_maxrss in bytes (macOS reports bytes, Linux KiB); nullopt when not positive.
std::optional<std::uint64_t> maxrss_to_bytes(long ru_maxrss, bool unit_is_bytes) noexcept;
/// getrusage(RUSAGE_SELF).ru_maxrss of this process in bytes, nullopt if unavailable.
std::optional<std::uint64_t> peak_rss_bytes() noexcept;

// ------------------------------------------------------------------ records

/// Names of bench.schema.ENV_FIELDS, in record order.
inline constexpr std::array<std::string_view, 13> kEnvFields = {
    "paper_version",       "git_sha", "git_dirty",   "os",         "cpu",
    "ram_bytes",           "arch",    "compiler",    "build_type", "build_flags",
    "diagnostics_enabled", "power",   "load_average"};

/// Parses --env-json text and returns an object with exactly kEnvFields (other keys are
/// ignored). Throws RunnerError(kExitUsage) if it is not a JSON object, a field is missing, the
/// paper version differs from this build's, or diagnostics_enabled disagrees with
/// `diagnostics_compiled` (an env file captured for another build).
Json load_env(std::string_view text, bool diagnostics_compiled);

struct GraphInfo {
  std::string graph_type;
  std::size_t n = 0;
  std::size_t m = 0;
  std::optional<double> avg_degree;
  std::string weight_dist;
  std::string weight_type;
  std::uint64_t seed = 0;
};

struct TimesNs {
  std::optional<std::uint64_t> graph_gen, preprocess, sssp, tree, end_to_end, verify;
};

struct Correctness {
  bool verified = false;
  std::optional<bool> ok;
  std::string method;
};

/// One SSSP run (or a run that failed with an error).
struct TrialRecord {
  std::string algorithm;
  std::string mode;
  Json params = Json::object();
  std::optional<Vertex> source;
  std::optional<std::size_t> trial;
  bool warmup = false;
  std::optional<std::uint64_t> order_seed;
  bool isolated = false;
  std::string status = "ok";
  std::optional<std::string> error;
  TimesNs times;
  std::optional<std::uint64_t> reachable;
  std::optional<std::uint64_t> peak_rss;
  std::optional<std::uint64_t> rss_after_graph;
  const Diagnostics* diag = nullptr;  ///< nullptr: no diagnostics (fallback = null)
  Correctness correctness;
};

/// Builds one schema record. ops / counters are written only if `diagnostics_compiled`;
/// peak_rss_bytes only if isolated; a failed verification forces status "error".
Json make_record(const Json& env, std::string_view run_id, std::string_view timestamp,
                 const GraphInfo& graph, const TrialRecord& trial, bool diagnostics_compiled);

/// "YYYY-MM-DDThh:mm:ss+00:00".
std::string utc_timestamp(std::chrono::system_clock::time_point t);

// ------------------------------------------------------------------ driver

/// --list output: graph types, weight distributions and the registry (name, mode, weight types).
Json list_json();

/// The whole runner: parses `args` (without argv[0]), runs every source x round x algorithm and
/// writes the records (CLI §3). --list / --version / --help and records without --out go to
/// `out`; diagnostics go to `err`. Returns the exit code of CLI §4.
int run_main(std::span<const std::string> args, std::ostream& out, std::ostream& err);

}  // namespace sssp::bench
