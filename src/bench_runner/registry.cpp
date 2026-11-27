// Algorithm registry, phase timer, and seeded source / order selection for the runner.

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/bench_runner.hpp"
#include "sssp/generators.hpp"
#include "sssp/paper.hpp"

namespace sssp::bench {

// ------------------------------------------------------------------ phase timer

std::uint64_t PhaseTimer::to_ns(Clock::duration d) noexcept {
  const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
  return ns < 1 ? std::uint64_t{1} : static_cast<std::uint64_t>(ns);
}

void PhaseTimer::add(Phase p, Clock::duration d) noexcept {
  auto& slot = ns_[static_cast<std::size_t>(p)];
  slot = slot.value_or(0) + to_ns(d);
}

// ------------------------------------------------------------------ registry

namespace {

// Baselines have no separable preprocessing, and build the parent tree while they settle
// vertices: only the sssp phase is timed (preprocess and tree stay null, see
// benchmark_protocol.md §4).
template <class W, SsspResult<W> (*Fn)(const CsrGraph<W>&, Vertex)>
SsspResult<W> run_baseline(const CsrGraph<W>& g, Vertex source, const Json& /*params*/,
                           PhaseTimer& timer) {
  return timer.time(PhaseTimer::Phase::kSssp, [&] { return Fn(g, source); });
}

// Mode A entry (Alg 8, sssp/paper.hpp). paper_sssp is one call: zero-edge contraction (l.2), the
// base case (l.5) and the tree restoration of Lemma 8.5 (l.6) are not separable from outside, so
// the whole call is the sssp phase and preprocess / tree stay null (benchmark_protocol.md §4).
// The FallbackRecord of the base case is returned in result.diag and written to the record's
// fallback block unchanged (record.cpp), in every build.
template <class W, paper::PaperVariant V>
SsspResult<W> run_paper(const CsrGraph<W>& g, Vertex source, const Json& /*params*/,
                        PhaseTimer& timer) {
  return timer.time(PhaseTimer::Phase::kSssp,
                    [&] { return paper::paper_sssp(g, source, paper::PaperOptions{V}); });
}

using U64 = std::uint64_t;
using paper::PaperVariant;

// exact_relaxation = false for the Mode A entry: it is not one of the relaxation baselines whose
// double distances are bit-identical to Dijkstra's (numerics.md §2); a non-baseline result is
// compared with a relaxation reference under Tolerance::for_path_hops(2(n-1)) (numerics.md §3.2,
// §3.4). Integer weights are compared exactly regardless (runner.cpp).
constexpr std::array<AlgorithmEntry, 7> kRegistry{{
    {"bellman-ford", "baseline", &run_baseline<U64, &baseline::bellman_ford<U64>>,
     &run_baseline<double, &baseline::bellman_ford<double>>, true, false, kBellmanFordMaxWork},
    {"dijkstra-binary", "baseline", &run_baseline<U64, &baseline::dijkstra_binary_heap<U64>>,
     &run_baseline<double, &baseline::dijkstra_binary_heap<double>>, true, false, 0},
    {"dijkstra-4ary", "baseline", &run_baseline<U64, &baseline::dijkstra_4ary_heap<U64>>,
     &run_baseline<double, &baseline::dijkstra_4ary_heap<double>>, true, false, 0},
    {"dijkstra-fibonacci", "baseline", &run_baseline<U64, &baseline::dijkstra_fibonacci_heap<U64>>,
     &run_baseline<double, &baseline::dijkstra_fibonacci_heap<double>>, true, false, 0},
    // Radix heap: monotone integer keys only (no double solver -> exit 3 for double weights).
    {"dijkstra-radix", "baseline", &run_baseline<U64, &baseline::dijkstra_radix_heap>, nullptr,
     true, false, 0},
    {paper::algorithm_name(PaperVariant::kFaithful), "paper-faithful",
     &run_paper<U64, PaperVariant::kFaithful>, &run_paper<double, PaperVariant::kFaithful>, false,
     false, 0, kScanBaseCaseMaxVertices},
    {paper::algorithm_name(PaperVariant::kFaithfulD1), "paper-faithful",
     &run_paper<U64, PaperVariant::kFaithfulD1>, &run_paper<double, PaperVariant::kFaithfulD1>,
     false, false, 0, 0},
}};

// bench.schema.ALGORITHM_MODES (docs/design/benchmark_schema.md §2.2).
constexpr std::array<std::pair<std::string_view, std::string_view>, 10> kSchemaModes{{
    {"dijkstra-binary", "baseline"},
    {"dijkstra-4ary", "baseline"},
    {"dijkstra-fibonacci", "baseline"},
    {"dijkstra-radix", "baseline"},
    {"bellman-ford", "baseline"},
    {"paper", "paper-faithful"},  // reserved for the outer construction (not implemented)
    {paper::algorithm_name(PaperVariant::kFaithful), "paper-faithful"},
    {paper::algorithm_name(PaperVariant::kFaithfulD1), "paper-faithful"},
    {"paper-experimental", "paper-experimental"},
    {"paper-optimized", "practical-optimized"},
}};

}  // namespace

std::span<const AlgorithmEntry> algorithm_registry() noexcept {
  return kRegistry;
}

const AlgorithmEntry* find_algorithm(std::string_view name) noexcept {
  for (const auto& e : kRegistry) {
    if (e.name == name) return &e;
  }
  return nullptr;
}

std::optional<std::string_view> schema_mode(std::string_view name) noexcept {
  for (const auto& [n, mode] : kSchemaModes) {
    if (n == name) return mode;
  }
  return std::nullopt;
}

// ------------------------------------------------------------------ sources and order

namespace {
// Sub-stream ids (the generators use their own ids on the same seed; these are distinct).
constexpr std::uint64_t kSourceStream = 0x534f55524345ULL;  // "SOURCE"
constexpr std::uint64_t kOrderStream = 0x4f52444552ULL;     // "ORDER"
}  // namespace

std::vector<Vertex> select_sources(std::size_t n, std::size_t k, std::uint64_t seed) {
  if (k > n) {
    throw std::invalid_argument("cannot choose " + std::to_string(k) + " distinct sources from " +
                                std::to_string(n) + " vertices");
  }
  // Sequential draws without replacement: the i-th source does not depend on k.
  gen::Xoshiro256ss rng(gen::substream_seed(seed, kSourceStream));
  std::vector<bool> taken(n, false);
  std::vector<Vertex> out;
  out.reserve(k);
  while (out.size() < k) {
    const auto v = static_cast<Vertex>(rng.uniform_below(n));
    if (taken[v]) continue;
    taken[v] = true;
    out.push_back(v);
  }
  return out;
}

std::vector<std::size_t> round_order(std::size_t count, std::uint64_t order_seed,
                                     std::uint64_t source_index, std::uint64_t round) {
  std::vector<std::size_t> order(count);
  std::iota(order.begin(), order.end(), std::size_t{0});
  const std::uint64_t seed = gen::substream_seed(
      gen::substream_seed(gen::substream_seed(order_seed, kOrderStream), source_index), round);
  gen::Xoshiro256ss rng(seed);
  for (std::size_t i = count; i > 1; --i) {  // Fisher-Yates
    const auto j = static_cast<std::size_t>(rng.uniform_below(i));
    std::swap(order[i - 1], order[j]);
  }
  return order;
}

}  // namespace sssp::bench
