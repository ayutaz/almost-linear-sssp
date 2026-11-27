#pragma once
// Google Benchmark driver for every algorithm of the bench_runner registry (P2-BENCH-CPP).
//
// Fairness rules (docs/design/benchmark_protocol.md, .claude/skills/benchmark):
//   * every algorithm runs on the same CsrGraph instance and the same sources: one graph per
//     (graph type, n, degree, weight distribution), built once before any benchmark runs, and
//     sources from bench::select_sources(n, k, seed), the selection sssp_bench_runner uses;
//   * graph generation, CSR construction, the reference run and the verification are outside the
//     timed region. The benchmarks use manual time: one iteration is the solver call alone
//     (the runner's end_to_end region, §4); destroying the result is not timed;
//   * a case is timed only after its result was checked against binary-heap Dijkstra and by
//     verify_sssp (shortest-path-tree certificate). A case that fails is reported as an error,
//     never timed, and makes the process exit with 1;
//   * algorithm x weight-type pairs without a solver (radix heap x double) are not registered;
//     pairs the registry refuses by size (max_vertices / max_work) are registered and skipped
//     with a message, so they are visible in the output instead of silently missing.
// Repetitions, warm-up and interleaving are Google Benchmark flags and therefore identical for
// every algorithm of one invocation.

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "sssp/bench_runner.hpp"
#include "sssp/graph.hpp"

namespace sssp::bench::gbench {

// ------------------------------------------------------------------ configuration

/// Driver options. They are passed as --sssp_<name>=<value> after the Google Benchmark flags
/// (benchmark::Initialize removes its own flags first); lists are comma separated.
struct DriverConfig {
  std::vector<std::string> graph_types{"er", "grid"};  ///< --sssp_graphs (sssp_gen --type names)
  std::vector<std::size_t> sizes{1000, 10000};         ///< --sssp_sizes (requested n)
  double avg_degree = 4.0;                             ///< --sssp_avg_degree
  /// --sssp_weights; the weight type follows from the distribution (weight_type_of). Families
  /// with built-in weights (road, adversarial_long_path) always use "builtin" (uint64).
  std::vector<std::string> weight_dists{"uniform_int", "uniform_real"};
  std::uint64_t seed = 1;               ///< --sssp_seed: graph seed and source-selection seed
  std::size_t num_sources = 1;          ///< --sssp_sources: sources per graph
  std::vector<std::string> algorithms;  ///< --sssp_algorithms; empty = whole registry
  bool verify_only = false;             ///< --sssp_verify_only: check every case, time nothing
  bool help = false;                    ///< --sssp_help
};

/// Parses the arguments left after benchmark::Initialize (without argv[0]). Throws
/// std::invalid_argument for an unknown option, a malformed value, an unknown algorithm or
/// weight distribution, or an empty list.
DriverConfig parse_driver_args(std::span<const std::string> args);

/// Text printed by --sssp_help.
std::string_view driver_usage() noexcept;

/// "uint64" or "double" for a weight distribution name (benchmark_protocol.md §12, the same
/// mapping as scripts/bench/plan.py WEIGHT_TYPES). Throws std::invalid_argument otherwise.
std::string_view weight_type_of(std::string_view weight_dist);

/// Graph families whose weights are fixed by construction (weight distribution "builtin").
bool has_builtin_weights(std::string_view graph_type) noexcept;

// ------------------------------------------------------------------ checking one case

struct CaseCheck {
  enum class Status : std::uint8_t {
    kOk,       ///< matches the reference and passes verify_sssp: may be timed
    kRefused,  ///< the registry refuses this size (max_vertices / max_work); skipped
    kFailed,   ///< wrong result, or the solver threw (e.g. integer overflow); error
  };
  Status status = Status::kFailed;
  std::string message;
  /// FallbackRecords of the checked run (Diagnostics.fallbacks), summed. Recorded in every
  /// build: the paper entries return their base-case record even with diagnostics off.
  std::uint64_t fallback_records = 0;
  std::uint64_t fallback_vertices = 0;
  std::uint64_t fallback_edges = 0;
};

/// Reason why `entry` refuses a graph with n vertices and m edges, or nullopt.
std::optional<std::string> refusal_reason(const AlgorithmEntry& entry, std::size_t n,
                                          std::size_t m);

/// Runs `entry` once on (g, source) and compares the result with `reference` (binary-heap
/// Dijkstra distances) through verify_sssp. Integer weights and relaxation baselines are compared
/// exactly; other double results with Tolerance::for_path_hops(2(n-1)) (numerics.md §3.2), the
/// rule sssp_bench_runner applies. Instantiated for std::uint64_t and double.
template <class W>
CaseCheck check_case(const AlgorithmEntry& entry, const CsrGraph<W>& g, Vertex source,
                     const std::vector<W>& reference);

// ------------------------------------------------------------------ the suite

enum class WeightType : std::uint8_t { kU64, kF64 };

/// One shared input: the graph every algorithm of the suite runs on, and its sources.
template <class W>
struct GraphInput {
  std::string label;  ///< "<type>/n:<n>[/d:<d>]/<weight_dist>/<weight_type>"
  std::string graph_type;
  std::string weight_dist;
  CsrGraph<W> graph;
  std::vector<Vertex> sources;
  /// Reference distances per source (binary-heap Dijkstra), computed on first use.
  std::vector<std::optional<std::vector<W>>> reference;
  std::vector<std::optional<std::string>> reference_error;
};

class Suite {
 public:
  struct Case {
    const AlgorithmEntry* algorithm = nullptr;
    WeightType weight_type = WeightType::kU64;
    std::size_t input = 0;         ///< index into inputs<W>()
    std::size_t source_index = 0;  ///< index into the input's sources
    std::string name;              ///< "<algorithm>/<input label>/src:<i>"
  };

  /// Builds every input (outside any timed region) and enumerates the cases: input x source x
  /// algorithm, skipping algorithms without a solver for the input's weight type. Throws
  /// std::invalid_argument if a graph cannot be generated or has fewer vertices than sources.
  explicit Suite(const DriverConfig& config);

  [[nodiscard]] const std::vector<Case>& cases() const noexcept { return cases_; }
  [[nodiscard]] const std::vector<GraphInput<std::uint64_t>>& inputs_u64() const noexcept {
    return u64_;
  }
  [[nodiscard]] const std::vector<GraphInput<double>>& inputs_f64() const noexcept { return f64_; }

  /// check_case for case i (cached: every case is run and verified at most once).
  const CaseCheck& check(std::size_t i);

  /// Registers one Google Benchmark per case (manual time) and returns the names in
  /// registration order. The suite must outlive benchmark::RunSpecifiedBenchmarks.
  std::vector<std::string> register_benchmarks();

  /// Checks every case and writes one line per case; true iff no case failed.
  bool verify_all(std::ostream& out);

  /// True iff a case checked so far failed.
  [[nodiscard]] bool any_failed() const noexcept;

 private:
  template <class W>
  std::vector<GraphInput<W>>& inputs() noexcept;
  template <class W>
  CaseCheck check_impl(const Case& c);
  template <class W>
  void add_input(const DriverConfig& config, const std::string& type, std::size_t n,
                 const std::string& dist);

  std::vector<const AlgorithmEntry*> algorithms_;
  std::vector<GraphInput<std::uint64_t>> u64_;
  std::vector<GraphInput<double>> f64_;
  std::vector<Case> cases_;
  std::vector<std::optional<CaseCheck>> checks_;
};

}  // namespace sssp::bench::gbench
