#pragma once
// Reproducible graph generators.
//
// Reproducibility contract: every generator is a pure function of (parameters, seed) and produces
// bit-identical output on every platform and standard library. To guarantee that we
//   * use our own PRNGs (splitmix64 for seeding, xoshiro256** for streams) instead of <random>
//     engines + std::*_distribution, whose algorithms are implementation-defined;
//   * derive integers with exact integer arithmetic (rejection sampling, no modulo bias);
//   * derive reals only from IEEE-754 basic operations (+ - * / sqrt, all correctly rounded),
//     with log/exp implemented in this library (det_log / det_exp) because libm's
//     transcendental functions are not correctly rounded and differ between platforms;
//   * never iterate over unordered containers when producing output.
// FMA contraction is disabled in the implementation file (it would change rounding).

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>

#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"

namespace sssp::gen {

// ------------------------------------------------------------------ PRNG

/// splitmix64 (Steele, Lea, Flood 2014; reference code by S. Vigna). Used to expand a 64-bit
/// seed into xoshiro state and to derive independent sub-stream seeds.
class SplitMix64 {
 public:
  explicit constexpr SplitMix64(std::uint64_t seed) noexcept : state_(seed) {}
  constexpr std::uint64_t next() noexcept {
    std::uint64_t z = (state_ += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }

 private:
  std::uint64_t state_;
};

/// xoshiro256** 1.0 (Blackman, Vigna 2018).
class Xoshiro256ss {
 public:
  using State = std::array<std::uint64_t, 4>;

  /// Seeds the 256-bit state with four splitmix64 outputs (the reference recommendation).
  explicit constexpr Xoshiro256ss(std::uint64_t seed) noexcept : s_{} {
    SplitMix64 sm(seed);
    for (auto& w : s_) w = sm.next();
  }
  /// Raw state (must not be all zero); used for known-answer tests.
  static constexpr Xoshiro256ss from_state(State s) noexcept {
    Xoshiro256ss r(0);
    r.s_ = s;
    return r;
  }

  constexpr std::uint64_t next() noexcept {
    const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }

  /// Uniform integer in [0, bound), bound >= 1. Unbiased: rejects the low (2^64 mod bound)
  /// values so that every residue class is equally likely.
  constexpr std::uint64_t uniform_below(std::uint64_t bound) noexcept {
    const std::uint64_t threshold = (0 - bound) % bound;  // == 2^64 mod bound
    for (;;) {
      const std::uint64_t r = next();
      if (r >= threshold) return r % bound;
    }
  }

  /// Uniform integer in [lo, hi] (inclusive), lo <= hi.
  constexpr std::uint64_t uniform_range(std::uint64_t lo, std::uint64_t hi) noexcept {
    const std::uint64_t span = hi - lo;
    if (span == ~std::uint64_t{0}) return next();
    return lo + uniform_below(span + 1);
  }

  /// Uniform double in [0, 1) with 53 random bits (exact conversion, no rounding).
  constexpr double uniform01() noexcept { return static_cast<double>(next() >> 11) * 0x1.0p-53; }

  /// Bernoulli(p).
  constexpr bool bernoulli(double p) noexcept { return uniform01() < p; }

 private:
  static constexpr std::uint64_t rotl(std::uint64_t x, int k) noexcept {
    return (x << k) | (x >> (64 - k));
  }
  State s_;
};

/// Derives the seed of an independent sub-stream (e.g. topology vs. weights) from a user seed.
constexpr std::uint64_t substream_seed(std::uint64_t seed, std::uint64_t stream) noexcept {
  SplitMix64 sm(seed ^ (0xD1B54A32D192ED03ULL * (stream + 1)));
  return sm.next();
}

// ------------------------------------------------------------------ deterministic math

/// Natural logarithm built only from correctly rounded IEEE operations (x > 0, finite).
/// Accurate to a few ulp; identical on every IEEE-754 platform.
double det_log(double x) noexcept;
/// Exponential built only from correctly rounded IEEE operations (|x| <= ~709).
double det_exp(double x) noexcept;

// ------------------------------------------------------------------ weight distributions

enum class WeightDist {
  kUniformInt,   ///< integer uniform in [ceil(lo), floor(hi)] (0 <= lo <= hi <= 2^53)
  kUniformReal,  ///< lo + (hi - lo) * U, U in [0,1)
  kLogUniform,   ///< exp(uniform(log lo, log hi)), 0 < lo <= hi
  kExponential,  ///< exponential with the given mean
  kConstant,     ///< always lo
  kZeroMixed,    ///< 0 with probability zero_probability, otherwise uniform in [lo, hi]
  kWideRange,    ///< binary exponent e uniform in [0, wide_log2], then uniform in that binade,
                 ///< capped at 2^wide_log2: covers [1, 2^wide_log2] (e.g. 1..2^50)
};

/// Parameters of a weight distribution. For integer W, real-valued draws are rounded to the
/// nearest integer (ties up); kUniformInt and kZeroMixed draw integers directly. For double W,
/// kUniformInt yields integral doubles and kZeroMixed draws reals in [lo, hi).
struct WeightSpec {
  WeightDist dist = WeightDist::kUniformInt;
  double lo = 1.0;
  double hi = 100.0;
  double mean = 10.0;              ///< kExponential
  double zero_probability = 0.25;  ///< kZeroMixed
  int wide_log2 = 50;              ///< kWideRange, 0..62
};

/// Throws std::invalid_argument when the parameters are inconsistent for spec.dist.
void validate_weight_spec(const WeightSpec& spec);

/// One draw from the distribution (spec must be valid). Explicitly instantiated for
/// std::uint64_t and double.
template <class W>
W sample_weight(const WeightSpec& spec, Xoshiro256ss& rng);

/// Canonical names used by the CLI: uniform_int, uniform_real, log_uniform, exponential,
/// constant, zero_mixed, wide_range.
std::string_view to_string(WeightDist d) noexcept;
/// Inverse of to_string; throws std::invalid_argument for an unknown name.
WeightDist parse_weight_dist(std::string_view name);

// ------------------------------------------------------------------ graph generators
//
// All generators return an EdgeList<W> (W = std::uint64_t or double) and are deterministic in
// (parameters, options). Topology, weights and connecting edges use three independent PRNG
// sub-streams derived from options.seed, so the same seed gives the same topology for every
// weight distribution and for both weight types.

struct GenOptions {
  WeightSpec weights{};
  std::uint64_t seed = 1;
  /// If true, components are joined afterwards by connect_components (random tree over the
  /// components, weights from `weights`). Generators that are connected by construction ignore it.
  bool connected = false;
};

/// Erdős–Rényi G(n, m) with m = round(n * avg_degree / 2): m distinct pairs drawn uniformly
/// (no self-loops, no multi-edges). Throws std::invalid_argument if m > n(n-1)/2.
template <class W>
EdgeList<W> erdos_renyi(std::size_t n, double avg_degree, const GenOptions& opt);
template <class W>
EdgeList<W> erdos_renyi_gnm(std::size_t n, std::size_t m, const GenOptions& opt);

/// What the configuration model does with self-loops and multi-edges.
enum class MultiEdgePolicy {
  kKeep,            ///< keep them (degrees exact, self-loop counts 2, graph not simple)
  kDrop,            ///< delete them (simple, degrees may fall below d)
  kSwitchToSimple,  ///< remove them by degree-preserving double-edge swaps (simple, exactly d)
};

/// Random d-regular graph via the configuration model (random perfect matching of n*d stubs).
/// Requires n*d even; kSwitchToSimple also requires d < n.
template <class W>
EdgeList<W> random_regular(std::size_t n, std::size_t d, MultiEdgePolicy policy,
                           const GenOptions& opt);

/// Barabási–Albert preferential attachment: a clique on m+1 vertices, then every new vertex
/// attaches to m distinct existing vertices chosen with probability proportional to degree.
/// Simple and connected; m * (n - m - 1) + m(m+1)/2 edges. Requires 1 <= m < n.
template <class W>
EdgeList<W> barabasi_albert(std::size_t n, std::size_t m, const GenOptions& opt);

/// rows x cols 4-neighbour grid; vertex (r, c) has id r * cols + c.
template <class W>
EdgeList<W> grid_2d(std::size_t rows, std::size_t cols, const GenOptions& opt);

struct RoadLikeOptions {
  double jitter = 0.3;            ///< each grid point moves by U[-jitter, jitter] per axis (< 0.5)
  double keep_probability = 0.9;  ///< each grid edge survives with this probability
  double shortcut_probability = 0.15;  ///< per vertex: add one local shortcut
  int shortcut_radius = 2;             ///< shortcut target within this Chebyshev radius
  double integer_scale = 1000.0;       ///< integer weights: round(length * scale), at least 1
};

/// Road-network-like graph: a perturbed grid with randomly removed streets and local
/// shortcuts, weighted by Euclidean length of the perturbed points (options.weights is ignored).
/// Always connected (removed streets are restored as needed). Simple.
template <class W>
EdgeList<W> road_like(std::size_t rows, std::size_t cols, const RoadLikeOptions& road,
                      const GenOptions& opt);

/// Uniformly random labelled tree on n vertices (random Prüfer sequence).
template <class W>
EdgeList<W> random_tree(std::size_t n, const GenOptions& opt);

/// Chung–Lu style scale-free graph: expected-degree weights w_i ∝ (i + 1)^(-1/(exponent-1))
/// (power-law degree tail with the given exponent > 2), m = round(n * avg_degree / 2) distinct
/// pairs, endpoints drawn independently ∝ w (self-loops and duplicates are re-drawn).
template <class W>
EdgeList<W> chung_lu(std::size_t n, double avg_degree, double exponent, const GenOptions& opt);

/// Decrease-key adversary for Dijkstra on k+1 vertices (source 0, v_i = i):
///   {0, 1} weight 1, {0, i} weight 2k for i > 1, {i, j} weight 2(j - i) - 1 for 1 <= i < j.
/// Rationale: dist(v_j) = j (chain of weight-1 edges; every edge {i,j} costs >= |j - i|), so
/// vertices are settled in order 1, 2, ..., k, and settling v_i lowers the tentative distance of
/// every v_j (j > i) from j + (j - i) to j + (j - i - 1). Hence a heap-based Dijkstra performs
/// exactly k(k-1)/2 = Theta(m) decrease-keys — the worst case for binary heaps (O(m log n)) and
/// the case Fibonacci heaps are designed for.
template <class W>
EdgeList<W> adversarial_decrease_key(std::size_t k);

/// Long-shortest-path adversary on n vertices: path {i, i+1} with weight 1 plus decoy skip edges
/// {i, i+2} with weight 3 (> 2). The shortest path from 0 to n-1 has n-1 hops, so any method
/// that bounds hop count, uses few rounds of Bellman–Ford, or recursion depth proportional to
/// path length is stressed, while every decoy is relaxed and rejected.
template <class W>
EdgeList<W> adversarial_long_path(std::size_t n);

/// Joins all connected components: components are ordered by smallest vertex, and component
/// i > 0 is linked to a uniformly random earlier component by an edge between uniformly random
/// members (a random tree over components). Appends exactly (components - 1) edges.
template <class W>
void connect_components(EdgeList<W>& g, const WeightSpec& weights, std::uint64_t seed);

/// Number of connected components (isolated vertices count).
template <class W>
std::size_t count_components(const EdgeList<W>& g);

// ------------------------------------------------------------------ generic entry point

/// Parameters for generate(); each type reads the fields listed in generator_catalog().
struct GraphSpec {
  std::string type = "er";
  std::size_t n = 1000;     ///< vertices (grid/road: rows = cols = floor(sqrt(n)) if unset)
  double avg_degree = 4.0;  ///< er, chung_lu
  std::size_t degree = 3;   ///< regular
  std::size_t attach = 2;   ///< ba
  std::size_t rows = 0;     ///< grid, road
  std::size_t cols = 0;     ///< grid, road
  double exponent = 2.5;    ///< chung_lu
  MultiEdgePolicy multi = MultiEdgePolicy::kSwitchToSimple;  ///< regular
  RoadLikeOptions road{};
  GenOptions options{};
};

struct GeneratorInfo {
  std::string_view name;
  std::string_view description;
  std::string_view parameters;
};

std::span<const GeneratorInfo> generator_catalog() noexcept;

/// Dispatches on spec.type; throws std::invalid_argument for an unknown type.
/// adversarial_decrease_key uses k = n - 1.
template <class W>
EdgeList<W> generate(const GraphSpec& spec);

MultiEdgePolicy parse_multi_edge_policy(std::string_view name);

}  // namespace sssp::gen
