#pragma once
// Support for the Las Vegas checks of task P4-RANDOMIZED (arXiv 2611.34911v1, Sec 6): instance
// families, checkers and statistics shared by the Mode A tests
// (tests/property/paper/test_randomized_las_vegas_property.cpp) and the Mode B tests
// (tests/property/paper_experimental/test_randomized_mode_b_property.cpp).
//
// Independence. The checkers in this file are written from Definition 6.1, Lemma 6.7(1) and
// Lemma 6.8(1); they do not call the library's check_route_names / check_estimate or its
// power-of-two arithmetic. Reference distances come from Floyd-Warshall (tests/support) and are
// cross-checked against the baseline Bellman-Ford on every instance.
//
// Exactness. uint64 weights: P * lambda >= d is evaluated as lambda >= ceil(d / 2^k) with shifts,
// so nothing can overflow. double weights: every weight is k / 8 with k < 2^20, all path sums
// stay below 2^53 / 8 and are therefore exact; scaling by a power of two (ldexp) is exact.
//
// Statistics. All thresholds below are derived for a false-positive probability of at most
// kDelta = 1e-9 per test under the lemma being checked, treating distinct WordRng seeds as
// independent randomness. The derivations are next to each function. Measured values are
// printed by the tests (Catch2 WARN); none is assumed.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numbers>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/pivots/estimate.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

namespace sssp::test::randomized {

namespace pv = sssp::paper::pivots;
using sssp::paper::model::WordRng;

/// False-positive budget of every statistical check in the P4-RANDOMIZED tests.
inline constexpr double kDelta = 1e-9;

// ------------------------------------------------------------------------------- weights

/// Integer k as a weight: k for integers, the dyadic rational k / 8 for double.
template <class W>
W weight_from(std::uint64_t k) {
  if constexpr (std::is_floating_point_v<W>) {
    return static_cast<W>(k) / W{8};
  } else {
    return static_cast<W>(k);
  }
}

template <class W>
const char* weight_name() {
  return std::is_floating_point_v<W> ? "double(dyadic k/8)" : "uint64";
}

/// lambda * 2^k >= d, exactly (see the file comment). Requires lambda > 0.
inline bool pow2_times_at_least(std::uint64_t lambda, unsigned k, std::uint64_t d) {
  if (k >= 64) return true;  // lambda >= 1 and d < 2^64
  const std::uint64_t q = d >> k;
  const bool remainder = (q << k) != d;
  return lambda >= q + (remainder ? 1U : 0U);  // q + 1 cannot wrap: remainder implies k >= 1
}
inline bool pow2_times_at_least(double lambda, unsigned k, double d) {
  return std::ldexp(lambda, static_cast<int>(k)) >= d;
}

// ------------------------------------------------------------------- independent constants

inline unsigned ceil_log2_u64(std::uint64_t x) {
  unsigned k = 0;
  while (k < 64 && (std::uint64_t{1} << k) < x) ++k;
  return k;
}

/// Alg 5 l.2-3: u = ceil(log2(t + 2)), g = 4u + 3, P = least power of two >= 8 * 256 * g * u^2.
inline unsigned paper_p_log2(std::uint64_t t) {
  const std::uint64_t u = ceil_log2_u64(t + 2);
  const std::uint64_t g = 4 * u + 3;
  return ceil_log2_u64(8 * 256 * g * u * u);
}

/// Sec 6.2 l.1491: H = least power of two >= (16t)^5. Valid for t <= 445 ((16t)^5 < 2^64).
inline unsigned paper_h_log2(std::uint64_t t) {
  REQUIRE(t <= 445);
  const std::uint64_t b = 16 * t;
  return ceil_log2_u64(b * b * b * b * b);
}

/// Deviation D-9 (Q-42): H = least power of two >= 4 t T^2 with T = least power of two >= 4t.
/// Since T^2 is a power of two, this is T^2 * (least power of two >= 4t) = T^3.
inline unsigned proof_minimum_h_log2(std::uint64_t t) {
  return 3 * ceil_log2_u64(4 * t);
}

inline double harmonic(std::uint64_t n) {
  double h = 0;
  for (std::uint64_t i = 1; i <= n; ++i) h += 1.0 / static_cast<double>(i);
  return h;
}

// ------------------------------------------------------------------------------ instances

template <class W>
using Instance = pv::PivotRoutingInstance<W>;

enum class Family : std::uint8_t {
  kComplete,            ///< K_t, weights uniform in 1..1000
  kCompleteEqual,       ///< K_t, every weight equal
  kCompleteFewWeights,  ///< K_t, weights in {1, 2, 3} (many ties)
  kCompleteWideRange,   ///< K_t, weights (1..9) * 23^(0..3) (several bands)
  kDense,               ///< each pair present with probability 4/5, weights 1..1000
  kMultiEdge,           ///< K_t with 1..3 parallel copies per pair, some of equal weight
  kTwoCliques,          ///< two disjoint cliques on floor(t/2) and ceil(t/2) pivots
  kSparseRandom,        ///< about 1.5 * C(t,2) random edges with repetitions, weights in decades
};

inline const char* family_name(Family f) {
  switch (f) {
    case Family::kComplete:
      return "complete";
    case Family::kCompleteEqual:
      return "complete-equal-weights";
    case Family::kCompleteFewWeights:
      return "complete-few-weights";
    case Family::kCompleteWideRange:
      return "complete-wide-range";
    case Family::kDense:
      return "dense-0.8";
    case Family::kMultiEdge:
      return "multi-edge";
    case Family::kTwoCliques:
      return "two-cliques";
    case Family::kSparseRandom:
      return "random-multigraph";
  }
  return "?";
}

template <class W>
Instance<W> make_instance(Family f, std::size_t t, std::uint64_t seed) {
  SplitMix64 g(seed ^ 0x5eed0fa417ULL);
  Instance<W> inst;
  inst.t = t;
  const auto add = [&](std::size_t x, std::size_t y, std::uint64_t k) {
    // Random orientation: the algorithm must treat (x, y) and (y, x) as the same pair.
    if (g.below(2) == 1) std::swap(x, y);
    inst.edges.push_back({static_cast<Vertex>(x), static_cast<Vertex>(y), weight_from<W>(k)});
  };
  const auto decade = [&] {
    std::uint64_t k = 1 + g.below(9);
    for (std::uint64_t i = g.below(4); i > 0; --i) k *= 23;
    return k;
  };
  const std::size_t half = t / 2;
  for (std::size_t x = 0; x < t; ++x) {
    for (std::size_t y = x + 1; y < t; ++y) {
      switch (f) {
        case Family::kComplete:
          add(x, y, 1 + g.below(1000));
          break;
        case Family::kCompleteEqual:
          add(x, y, 7);
          break;
        case Family::kCompleteFewWeights:
          add(x, y, 1 + g.below(3));
          break;
        case Family::kCompleteWideRange:
          add(x, y, decade());
          break;
        case Family::kDense:
          if (g.below(5) != 0) add(x, y, 1 + g.below(1000));
          break;
        case Family::kMultiEdge: {
          const std::uint64_t first = 1 + g.below(1000);
          add(x, y, first);
          for (std::uint64_t c = g.below(3); c > 0; --c) {
            add(x, y, g.below(3) == 0 ? first : 1 + g.below(1000));  // equal-weight copies too
          }
          break;
        }
        case Family::kTwoCliques:
          if ((x < half) == (y < half)) add(x, y, 1 + g.below(1000));
          break;
        case Family::kSparseRandom:
          break;
      }
    }
  }
  if (f == Family::kSparseRandom && t >= 2) {
    const std::size_t m = 3 * (t * (t - 1) / 2) / 2 + static_cast<std::size_t>(g.below(t + 1));
    for (std::size_t i = 0; i < m; ++i) {
      const std::size_t x = static_cast<std::size_t>(g.below(t));
      std::size_t y = static_cast<std::size_t>(g.below(t));
      if (x == y) y = (y + 1) % t;
      add(x, y, decade());
    }
  }
  return inst;
}

/// A family at a fixed number of pivots.
struct FamilySpec {
  Family family;
  std::size_t t;
};

/// Families on which the paper's loop (Alg 5 l.5, |U| > 32t) runs: more than 32t distinct pairs,
/// which needs t >= 66 (C(66, 2) = 2145 > 2112), and at most 64t, so that one accepted round
/// ends the loop (Alg 5 l.13-14: an accepted round removes at least half of U).
inline constexpr FamilySpec kLoopFamilies[] = {
    {Family::kComplete, 66},
    {Family::kCompleteEqual, 68},
    {Family::kCompleteFewWeights, 70},
    {Family::kCompleteWideRange, 72},
    {Family::kDense, 90},
    {Family::kMultiEdge, 67},
    {Family::kTwoCliques, 132},
};

/// Seeds per family of the [slow] sweeps.
inline constexpr std::uint64_t kSweepSeeds = 1000;
/// Seeds per family of the [slow] sweeps whose single run is expensive: NameRoutes with uint64
/// weights at t >= 66 (Estimate then works on multiword integers; measured at about 0.1 s per
/// run optimised and 0.7 s unoptimised). Optimised builds (NDEBUG: release, ci-release and the
/// sanitizer preset) run the full 1000 seeds; unoptimised builds run 125 so that the cases stay
/// inside the ctest timeouts. The count used is printed with the results.
#ifdef NDEBUG
inline constexpr std::uint64_t kHeavySweepSeeds = 1000;
#else
inline constexpr std::uint64_t kHeavySweepSeeds = 125;
#endif

/// Seeds of a [slow] NameRoutes sweep at t >= 66 for weight type W.
template <class W>
constexpr std::uint64_t loop_sweep_seeds() {
  return std::is_floating_point_v<W> ? kSweepSeeds : kHeavySweepSeeds;
}

/// Number of distinct endpoint pairs (the representatives of Alg 5 l.4), counted independently.
template <class W>
std::uint64_t distinct_pairs(const Instance<W>& inst) {
  std::vector<std::uint64_t> keys;
  keys.reserve(inst.edges.size());
  for (const auto& e : inst.edges) {
    keys.push_back(std::uint64_t{std::min(e.u, e.v)} * inst.t + std::max(e.u, e.v));
  }
  std::sort(keys.begin(), keys.end());
  return static_cast<std::uint64_t>(std::unique(keys.begin(), keys.end()) - keys.begin());
}

/// Reference all-pairs distances: Floyd-Warshall, required to equal the baseline Bellman-Ford
/// run from every source (two independent oracles; exact for both weight types, see the file
/// comment).
template <class W>
std::vector<std::vector<W>> reference_distances(std::size_t t, const std::vector<Edge<W>>& edges) {
  auto fw = floyd_warshall<W>(t, std::span<const Edge<W>>(edges));
  const auto g = CsrGraph<W>::from_edges(t, edges);
  for (Vertex s = 0; s < t; ++s) {
    const auto bf = sssp::baseline::bellman_ford(g, s).dist;
    REQUIRE(bf == fw[s]);
  }
  return fw;
}

// ------------------------------------------------------------- certificates (Definition 6.1)

/// Definition 6.1 / Lemma 6.8(1) against reference distances of G_P: one name per edge, inside
/// the dictionary; the dictionary strictly increasing and made of input weights;
/// lambda_e <= a_e; dist(x, y) <= 2^p_log2 * lambda_e. Returns the violations (empty: valid).
template <class W>
std::vector<std::string> certificate_violations(const Instance<W>& inst,
                                                const pv::RouteNames<W>& names,
                                                const std::vector<std::vector<W>>& dist) {
  std::vector<std::string> out;
  const auto fail = [&](const std::string& what, std::size_t i) {
    if (out.size() < 8) out.push_back(what + " at index " + std::to_string(i));
  };
  if (names.name.size() != inst.edges.size()) {
    out.emplace_back("name count differs from the edge count");
    return out;
  }
  std::vector<W> weights;
  weights.reserve(inst.edges.size());
  for (const auto& e : inst.edges) weights.push_back(e.weight);
  std::sort(weights.begin(), weights.end());
  for (std::size_t i = 0; i < names.dictionary.size(); ++i) {
    if (i > 0 && !(names.dictionary[i - 1] < names.dictionary[i])) {
      fail("dictionary not strictly increasing", i);
    }
    if (!std::binary_search(weights.begin(), weights.end(), names.dictionary[i])) {
      fail("dictionary entry is not an input weight", i);
    }
  }
  for (std::size_t i = 0; i < inst.edges.size(); ++i) {
    const auto& e = inst.edges[i];
    if (names.name[i] >= names.dictionary.size()) {
      fail("name outside the dictionary", i);
      continue;
    }
    const W lambda = names.dictionary[names.name[i]];
    if (!(lambda <= e.weight)) fail("lambda_e > a_e", i);
    const W d = dist[e.u][e.v];
    if (d == WeightTraits<W>::infinity()) {
      fail("reference distance between the endpoints of an edge is infinite", i);
    } else if (!pow2_times_at_least(lambda, names.p_log2, d)) {
      fail("dist(x, y) > P * lambda_e", i);
    }
  }
  return out;
}

inline std::uint64_t counter(const Diagnostics& d, const std::string& key) {
  const auto it = d.counters.find(key);
  return it == d.counters.end() ? 0 : it->second;
}

using Histogram = std::map<std::uint64_t, std::uint64_t>;

inline std::string describe(const Histogram& h) {
  std::ostringstream s;
  bool first = true;
  for (const auto& [k, v] : h) {
    s << (first ? "" : " ") << k << ":" << v;
    first = false;
  }
  return first ? "(empty)" : s.str();
}

/// What one NameRoutes run was asked to do (taken from the parameters, not from its output).
struct RunSetting {
  std::uint64_t residual_limit = 0;        ///< loop while |U| > this
  std::optional<std::uint64_t> retry_cap;  ///< D-10
  unsigned p_log2 = 0;                     ///< expected log2 P
};

/// Measured trial statistics of a sweep. Histograms are per run.
struct Tally {
  std::uint64_t runs = 0;
  std::uint64_t runs_with_trials = 0;
  std::uint64_t trials = 0;
  std::uint64_t accepted = 0;
  std::uint64_t rejected = 0;
  std::uint64_t runs_without_rejection = 0;  ///< among runs with trials
  std::uint64_t cap_hits = 0;
  std::uint64_t fallback_edges = 0;
  std::uint64_t edges = 0;
  std::uint64_t dictionary_total = 0;
  std::uint64_t dictionary_max = 0;
  Histogram trials_per_run;
  Histogram accepted_rounds_per_run;
  Histogram max_trials_per_round;
  Histogram max_consecutive_rejections;

  void merge(const Tally& o) {
    runs += o.runs;
    runs_with_trials += o.runs_with_trials;
    trials += o.trials;
    accepted += o.accepted;
    rejected += o.rejected;
    runs_without_rejection += o.runs_without_rejection;
    cap_hits += o.cap_hits;
    fallback_edges += o.fallback_edges;
    edges += o.edges;
    dictionary_total += o.dictionary_total;
    dictionary_max = std::max(dictionary_max, o.dictionary_max);
    for (const auto& [k, v] : o.trials_per_run) trials_per_run[k] += v;
    for (const auto& [k, v] : o.accepted_rounds_per_run) accepted_rounds_per_run[k] += v;
    for (const auto& [k, v] : o.max_trials_per_round) max_trials_per_round[k] += v;
    for (const auto& [k, v] : o.max_consecutive_rejections) max_consecutive_rejections[k] += v;
  }
};

/// Checks one run (certificates against `dist`, the loop bookkeeping of Alg 5 and of D-10) and
/// adds its trial counts to `tally`.
template <class W>
void check_run(const Instance<W>& inst, const std::vector<std::vector<W>>& dist,
               const pv::RouteNames<W>& names, const RunSetting& setting, Tally& tally) {
  const auto violations = certificate_violations(inst, names, dist);
  INFO("t " << inst.t << " M " << inst.edges.size() << " P 2^" << names.p_log2
            << (violations.empty() ? "" : ": " + violations.front()));
  REQUIRE(violations.empty());
  CHECK(names.p_log2 == setting.p_log2);

  const auto& d = names.diag;
  const std::uint64_t reps = distinct_pairs(inst);
  const std::uint64_t trials = counter(d, "name_routes.trials");
  const std::uint64_t accepted = counter(d, "name_routes.accepted_rounds");
  const std::uint64_t rejected = counter(d, "name_routes.rejected_trials");
  const std::uint64_t terminal = counter(d, "name_routes.terminal_edges");
  const std::uint64_t consecutive = counter(d, "name_routes.max_consecutive_rejections");
  const std::uint64_t round_max = counter(d, "name_routes.max_trials_per_round");
  CHECK(trials == accepted + rejected);
  CHECK(terminal <= reps);
  // Alg 5 l.5: a trial runs iff more than `limit` representatives remain.
  CHECK((trials > 0) == (reps > setting.residual_limit));
  // Alg 5 l.13-14: an accepted round removes at least half of U, and it ran with |U| > limit, so
  // the k-th accepted round saw limit < |U| <= reps / 2^(k-1).
  if (accepted >= 1) {
    REQUIRE(accepted < 64);
    CHECK((std::uint64_t{1} << (accepted - 1)) * (setting.residual_limit + 1) <= reps);
  }
  const bool capped = !d.fallbacks.empty();
  if (capped) {
    // D-10: recorded once, with the residual it left; only a cap can leave |U| above the limit.
    REQUIRE(setting.retry_cap.has_value());
    REQUIRE(d.fallbacks.size() == 1);
    CHECK(d.fallbacks[0].vertices == inst.t);
    CHECK(d.fallbacks[0].edges == terminal);
    CHECK(terminal > setting.residual_limit);
    CHECK(consecutive == *setting.retry_cap);
    CHECK(counter(d, "name_routes.retry_cap_hits") == 1);
  } else {
    CHECK(terminal <= setting.residual_limit);
    CHECK(counter(d, "name_routes.retry_cap_hits") == 0);
    if (setting.retry_cap.has_value()) CHECK(consecutive < *setting.retry_cap);
  }

  ++tally.runs;
  tally.edges += inst.edges.size();
  tally.trials += trials;
  tally.accepted += accepted;
  tally.rejected += rejected;
  tally.dictionary_total += names.dictionary.size();
  tally.dictionary_max = std::max<std::uint64_t>(tally.dictionary_max, names.dictionary.size());
  ++tally.trials_per_run[trials];
  ++tally.accepted_rounds_per_run[accepted];
  if (trials > 0) {
    ++tally.runs_with_trials;
    if (rejected == 0) ++tally.runs_without_rejection;
    ++tally.max_trials_per_round[round_max];
    ++tally.max_consecutive_rejections[consecutive];
  }
  if (capped) {
    ++tally.cap_hits;
    tally.fallback_edges += terminal;
  }
}

/// The WordRng seed of run `seed` of a sweep; `salt` separates sweeps.
inline std::uint64_t rng_seed(std::uint64_t seed, std::uint64_t salt) {
  return seed * 0x9E3779B97F4A7C15ULL + salt;
}

// ------------------------------------------------------------------ acceptance statistics

/// Test of Lemma 6.8(2), "every trial is accepted with probability >= 1/2".
///
/// Let X_1, X_2, ... be the accept indicators of all trials of a sweep in execution order. The
/// lemma bounds the conditional probability of each trial given everything before it (the
/// residual U is arbitrary in its proof), so under the lemma S_n = sum_{i<=n} (X_i - 1/2) is a
/// submartingale with increments in an interval of length 1, and Azuma-Hoeffding gives
///   Pr[A_n <= n/2 - s] <= exp(-2 s^2 / n)        (A_n = accepted among the first n).
/// The number of trials of a sweep is itself random (a run stops when its rounds are done), so
/// the bound is applied to every n at once with delta_n = 6 delta / (pi^2 n^2), whose sum over
/// n >= 1 is delta:
///   radius(n) = sqrt(ln(pi^2 n^2 / (6 delta)) / (2 n)).
/// `rejects_half` is true iff the measured rate is below 1/2 - radius(n). Under the lemma this
/// happens with probability at most delta = kDelta, whatever the stopping rule. `lower` is the
/// matching lower confidence bound on the average conditional acceptance probability.
struct AcceptanceTest {
  double rate = 1.0;
  double radius = 1.0;
  double lower = 0.0;
  bool rejects_half = false;
};

inline AcceptanceTest acceptance_test(std::uint64_t accepted, std::uint64_t trials) {
  AcceptanceTest r;
  if (trials == 0) return r;
  const auto n = static_cast<double>(trials);
  r.rate = static_cast<double>(accepted) / n;
  r.radius =
      std::sqrt(std::log(std::numbers::pi * std::numbers::pi * n * n / (6.0 * kDelta)) / (2.0 * n));
  r.lower = r.rate - r.radius;
  r.rejects_half = r.rate < 0.5 - r.radius;
  return r;
}

/// Fixed-n Hoeffding radius for the mean of n independent variables with range 1:
/// Pr[mean <= E - r] <= exp(-2 n r^2) = delta for r = sqrt(ln(1 / delta) / (2 n)).
inline double hoeffding_radius(std::uint64_t n, double delta = kDelta) {
  return n == 0 ? 1.0 : std::sqrt(std::log(1.0 / delta) / (2.0 * static_cast<double>(n)));
}

inline std::string describe(const std::string& label, const Tally& t) {
  const auto a = acceptance_test(t.accepted, t.trials);
  std::ostringstream s;
  s << label << ": runs " << t.runs << " (with trials " << t.runs_with_trials << "), trials "
    << t.trials << ", accepted " << t.accepted << ", rejected " << t.rejected
    << ", acceptance rate " << a.rate << " (anytime radius " << a.radius << ", lower bound "
    << a.lower << "), runs without rejection " << t.runs_without_rejection << ", cap hits "
    << t.cap_hits << ", fallback edges " << t.fallback_edges << " of " << t.edges
    << ", dictionary mean "
    << (t.runs == 0 ? 0.0 : static_cast<double>(t.dictionary_total) / static_cast<double>(t.runs))
    << " max " << t.dictionary_max << "\n    trials per run {" << describe(t.trials_per_run)
    << "}\n    accepted rounds per run {" << describe(t.accepted_rounds_per_run)
    << "}\n    max trials in one round, per run {" << describe(t.max_trials_per_round)
    << "}\n    max consecutive rejections, per run {" << describe(t.max_consecutive_rejections)
    << "}";
  return s.str();
}

// ------------------------------------------------------------------- Estimate (Lemma 6.7)

enum class SampleFamily : std::uint8_t {
  kSparseDecades,  ///< t in 2..31, t..3t-1 random edges, weights (1..9) * 37^(0..3)
  kSparseEqual,    ///< the same topology, every weight equal
  kMultiEdge,      ///< t in 2..16, each chosen pair 1..3 times, weights 1..40 or repeated
  kComplete,       ///< K_t for t in 2..14, weights 1..50
  kDenseTies,      ///< t in 4..20, each pair with probability 1/2, weights in {1, 2}
};

inline const char* sample_family_name(SampleFamily f) {
  switch (f) {
    case SampleFamily::kSparseDecades:
      return "sparse-decades";
    case SampleFamily::kSparseEqual:
      return "sparse-equal-weights";
    case SampleFamily::kMultiEdge:
      return "multi-edge";
    case SampleFamily::kComplete:
      return "complete";
    case SampleFamily::kDenseTies:
      return "dense-ties";
  }
  return "?";
}

template <class W>
struct Sample {
  std::size_t t = 0;
  std::vector<Edge<W>> edges;
};

template <class W>
Sample<W> make_sample(SampleFamily f, std::uint64_t seed) {
  SplitMix64 g(seed ^ 0x5a3b1e5ULL);
  Sample<W> s;
  const auto add = [&](std::size_t x, std::size_t y, std::uint64_t k) {
    if (g.below(2) == 1) std::swap(x, y);
    s.edges.push_back({static_cast<Vertex>(x), static_cast<Vertex>(y), weight_from<W>(k)});
  };
  const auto random_pair = [&](std::size_t& x, std::size_t& y) {
    x = static_cast<std::size_t>(g.below(s.t));
    y = static_cast<std::size_t>(g.below(s.t));
    if (x == y) y = (y + 1) % s.t;
  };
  switch (f) {
    case SampleFamily::kSparseDecades:
    case SampleFamily::kSparseEqual: {
      s.t = 2 + static_cast<std::size_t>(g.below(30));
      const std::size_t m = s.t + static_cast<std::size_t>(g.below(2 * s.t));
      for (std::size_t i = 0; i < m; ++i) {
        std::size_t x = 0;
        std::size_t y = 0;
        random_pair(x, y);
        std::uint64_t k = 1 + g.below(9);
        for (std::uint64_t d = g.below(4); d > 0; --d) k *= 37;
        add(x, y, f == SampleFamily::kSparseEqual ? 5 : k);
      }
      break;
    }
    case SampleFamily::kMultiEdge: {
      s.t = 2 + static_cast<std::size_t>(g.below(15));
      const std::size_t pairs = 1 + static_cast<std::size_t>(g.below(2 * s.t));
      for (std::size_t i = 0; i < pairs; ++i) {
        std::size_t x = 0;
        std::size_t y = 0;
        random_pair(x, y);
        const std::uint64_t first = 1 + g.below(40);
        add(x, y, first);
        for (std::uint64_t c = g.below(3); c > 0; --c) {
          add(x, y, g.below(2) == 0 ? first : 1 + g.below(40));
        }
      }
      break;
    }
    case SampleFamily::kComplete: {
      s.t = 2 + static_cast<std::size_t>(g.below(13));
      for (std::size_t x = 0; x < s.t; ++x) {
        for (std::size_t y = x + 1; y < s.t; ++y) add(x, y, 1 + g.below(50));
      }
      break;
    }
    case SampleFamily::kDenseTies: {
      s.t = 4 + static_cast<std::size_t>(g.below(17));
      for (std::size_t x = 0; x < s.t; ++x) {
        for (std::size_t y = x + 1; y < s.t; ++y) {
          if (g.below(2) == 0) add(x, y, 1 + g.below(2));
        }
      }
      break;
    }
  }
  return s;
}

inline std::vector<pv::Demand> all_pairs(std::size_t t) {
  std::vector<pv::Demand> q;
  for (std::size_t x = 0; x < t; ++x) {
    for (std::size_t y = x + 1; y < t; ++y) {
      q.push_back({static_cast<Vertex>(x), static_cast<Vertex>(y)});
    }
  }
  return q;
}

/// Lemma 6.7 against reference sample distances: one index per demand; a demand is flagged iff
/// its endpoints are disconnected; otherwise its index is inside the table, the anchor is a
/// sample weight and d_S(x, y) <= Z = 4 * 2^scale * anchor + anchor. The last test is evaluated
/// as d <= anchor or 2^(scale+2) * anchor >= d - anchor, which is exact (file comment).
template <class W>
std::vector<std::string> estimate_violations(const Sample<W>& s,
                                             std::span<const pv::Demand> demands,
                                             const pv::EstimateResult<W>& r,
                                             const std::vector<std::vector<W>>& dist) {
  std::vector<std::string> out;
  const auto fail = [&](const std::string& what, std::size_t i) {
    if (out.size() < 8) out.push_back(what + " at demand " + std::to_string(i));
  };
  if (r.index.size() != demands.size()) {
    out.emplace_back("index count differs from the demand count");
    return out;
  }
  std::vector<W> weights;
  for (const auto& e : s.edges) weights.push_back(e.weight);
  std::sort(weights.begin(), weights.end());
  for (std::size_t i = 0; i < demands.size(); ++i) {
    const W d = dist[demands[i].x][demands[i].y];
    const bool connected = d != WeightTraits<W>::infinity();
    if (!connected) {
      if (r.index[i] != pv::kDisconnected) fail("disconnected pair not flagged", i);
      continue;
    }
    if (r.index[i] == pv::kDisconnected) {
      fail("connected pair flagged as disconnected", i);
      continue;
    }
    if (r.index[i] >= r.table.size()) {
      fail("index outside the table", i);
      continue;
    }
    const auto& z = r.table[r.index[i]];
    if (!std::binary_search(weights.begin(), weights.end(), z.anchor)) {
      fail("table anchor is not a sample weight", i);
    }
    if (!(d <= z.anchor) &&
        !pow2_times_at_least(z.anchor, z.scale + 2, static_cast<W>(d - z.anchor))) {
      fail("Z < d_S", i);
    }
  }
  return out;
}

/// Measured quantities of an Estimate sweep.
struct EstimateTally {
  std::uint64_t runs = 0;
  std::uint64_t demands = 0;
  std::uint64_t connected = 0;
  std::uint64_t disconnected = 0;
  std::uint64_t table_max = 0;
  double normalised_ratio_sum = 0;  ///< sum over connected demands of Z / (256 u^2 d_S)
  double ratio_sum = 0;             ///< sum over connected demands of Z / d_S
  double ratio_max = 0;
  std::uint64_t scans = 0;        ///< estimate.adjacency_scans
  double scan_allowed = 0;        ///< sum of H_t * estimate.vertex_scales
  double scan_range_squares = 0;  ///< sum of (t * vertex_scales)^2, for the Hoeffding slack
  std::uint64_t vertex_scales = 0;

  [[nodiscard]] double mean_normalised_ratio() const {
    return connected == 0 ? 0.0 : normalised_ratio_sum / static_cast<double>(connected);
  }
  [[nodiscard]] double mean_ratio() const {
    return connected == 0 ? 0.0 : ratio_sum / static_cast<double>(connected);
  }
  /// Lemma 6.6: E[scans of run i] <= H_t * vertex_scales_i, and 0 <= scans_i <= t *
  /// vertex_scales_i (a vertex is scanned at most once per center and scale). The runs are
  /// independent, so by Hoeffding
  ///   Pr[sum scans > sum allowed + s] <= exp(-2 s^2 / sum (t * vertex_scales_i)^2),
  /// which is kDelta for s = sqrt(ln(1 / kDelta) * sum (t * vertex_scales_i)^2 / 2). The ranges
  /// depend only on the samples, not on the estimator's randomness.
  [[nodiscard]] double scan_slack() const {
    return std::sqrt(std::log(1.0 / kDelta) * scan_range_squares / 2.0);
  }
};

inline std::string describe(const std::string& label, const EstimateTally& t) {
  std::ostringstream s;
  s << label << ": runs " << t.runs << ", demands " << t.demands << " (connected " << t.connected
    << ", disconnected " << t.disconnected << "), largest table " << t.table_max << ", mean Z/d_S "
    << t.mean_ratio() << " max " << t.ratio_max << ", mean Z/(256 u^2 d_S) "
    << t.mean_normalised_ratio() << ", adjacency scans " << t.scans << " vs H_t-weighted bound "
    << t.scan_allowed << " (ratio "
    << (t.scan_allowed > 0 ? static_cast<double>(t.scans) / t.scan_allowed : 0.0)
    << ", Hoeffding slack " << t.scan_slack() << ")";
  return s.str();
}

/// Runs Estimate on `seeds` samples of a family with all pairs as demands and fresh randomness,
/// checks Lemma 6.7(1) on every run and accumulates the measured quantities. `params_of(t)`
/// yields the EstimateParams; `eta_of(t)` the expected log2 H (computed independently).
template <class W, class ParamsOf, class EtaOf>
EstimateTally sweep_estimate(SampleFamily family, std::uint64_t seeds, std::uint64_t salt,
                             ParamsOf params_of, EtaOf eta_of) {
  EstimateTally tally;
  for (std::uint64_t seed = 1; seed <= seeds; ++seed) {
    const auto s = make_sample<W>(family, seed);
    const auto q = all_pairs(s.t);
    const auto dist = reference_distances<W>(s.t, s.edges);
    const pv::EstimateParams params = params_of(s.t);
    WordRng rng(rng_seed(seed, salt));
    const auto r = pv::estimate<W>(s.t, s.edges, q, rng, params);
    const auto violations = estimate_violations<W>(s, q, r, dist);
    INFO(sample_family_name(family)
         << " " << weight_name<W>() << " seed " << seed << " t " << s.t << " edges "
         << s.edges.size() << (violations.empty() ? "" : ": " + violations.front()));
    REQUIRE(violations.empty());
    CHECK(pv::estimate_constants(s.t, params).eta == eta_of(s.t));

    const auto u = static_cast<double>(ceil_log2_u64(s.t + 2));
    ++tally.runs;
    tally.demands += q.size();
    tally.table_max = std::max<std::uint64_t>(tally.table_max, r.table.size());
    for (std::size_t k = 0; k < q.size(); ++k) {
      if (r.index[k] == pv::kDisconnected) {
        ++tally.disconnected;
        continue;
      }
      ++tally.connected;
      const auto& z = r.table[r.index[k]];
      const double zv = static_cast<double>(z.anchor) +
                        std::ldexp(static_cast<double>(z.anchor), static_cast<int>(z.scale) + 2);
      const double ratio = zv / static_cast<double>(dist[q[k].x][q[k].y]);
      tally.ratio_sum += ratio;
      tally.ratio_max = std::max(tally.ratio_max, ratio);
      tally.normalised_ratio_sum += ratio / (256.0 * u * u);
    }
    const std::uint64_t vs = counter(r.diag, "estimate.vertex_scales");
    tally.scans += counter(r.diag, "estimate.adjacency_scans");
    tally.vertex_scales += vs;
    tally.scan_allowed += harmonic(s.t) * static_cast<double>(vs);
    const double range = static_cast<double>(s.t) * static_cast<double>(vs);
    tally.scan_range_squares += range * range;
  }
  return tally;
}

/// The assertions on an Estimate sweep that are statistical (Lemma 6.6, Lemma 6.7(2)).
inline void check_estimate_tally(const EstimateTally& t) {
  REQUIRE(t.connected > 0);
  // The searches ran and were counted (the counters are written with SSSP_DIAGNOSTICS = 0 too).
  REQUIRE(t.scans > 0);
  REQUIRE(t.vertex_scales > 0);
  // Lemma 6.6, with a false-positive probability of at most kDelta (EstimateTally::scan_slack).
  CHECK(static_cast<double>(t.scans) <= t.scan_allowed + t.scan_slack());
  // The same bound without slack, in the form of the existing test of P3-ESTIMATE
  // (test_estimate_property.cpp). It has no a-priori false-positive bound; it is kept because
  // Q-64 asks for exactly this check, and the measured ratio is printed next to it.
  CHECK(static_cast<double>(t.scans) <= t.scan_allowed);
  // Lemma 6.7(2): E[Z | S] <= 256 u^2 d_S, as the mean of Z / (256 u^2 d_S) over all connected
  // demands. Z / d_S has no useful a-priori range (the top scale is t T^2 times the anchor), so
  // no concentration bound is claimed for this check; the measured mean is printed.
  CHECK(t.mean_normalised_ratio() <= 1.0);
}

// ------------------------------------------------------------- Lemma 6.5 (one-scale cuts)

/// Measured result of the separation test.
struct SeparationResult {
  std::uint64_t hypotheses = 0;  ///< (band, scale, pair) triples tested
  std::uint64_t nontrivial = 0;  ///< those whose bound is below 1
  std::uint64_t violations = 0;  ///< frequency > bound + radius
  double radius = 0;             ///< Hoeffding radius per hypothesis (union bound)
  double worst_excess = -1.0;    ///< max of frequency - bound
  double worst_ratio = 0;        ///< max of frequency / bound over nontrivial bounds
  std::string worst;             ///< description of the worst-excess triple
};

/// Lemma 6.5: two quotient vertices at quotient distance d receive different centers at scale
/// Delta = 2^i a with probability <= H_t (d / Delta + 1 / H). For one fixed sample the event is
/// counted over `seeds` independent runs for every band, scale and connected pair of quotient
/// vertices (the quotients depend only on the sample and the demands). With K triples, each
/// frequency is the mean of `seeds` independent indicators, so by Hoeffding and a union bound
///   Pr[some frequency > bound + r] <= K exp(-2 seeds r^2) = kDelta
/// for r = sqrt(ln(K / kDelta) / (2 seeds)). Quotient distances are computed here by
/// Floyd-Warshall on the traced quotient edges.
template <class ParamsOf>
SeparationResult separation_test(const Sample<std::uint64_t>& s, std::uint64_t seeds,
                                 std::uint64_t salt, ParamsOf params_of) {
  using W = std::uint64_t;
  const auto q = all_pairs(s.t);
  const pv::EstimateParams params = params_of(s.t);
  struct BandRef {
    W anchor = 0;
    std::size_t nq = 0;
    std::vector<std::vector<W>> qd;
    std::vector<std::uint32_t> scale;                   // per scale slot
    std::vector<std::vector<std::uint64_t>> separated;  // [slot][x * nq + y], x < y
  };
  std::vector<BandRef> bands;
  unsigned eta = 0;
  for (std::uint64_t seed = 1; seed <= seeds; ++seed) {
    WordRng rng(rng_seed(seed, salt));
    pv::EstimateTrace<W> tr;
    (void)pv::estimate<W>(s.t, s.edges, q, rng, params, &tr);
    if (seed == 1) {
      eta = tr.constants.eta;
      for (const auto& b : tr.bands) {
        BandRef ref;
        ref.anchor = b.anchor;
        ref.nq = b.num_quotient_vertices;
        ref.qd = floyd_warshall<W>(ref.nq, b.quotient_edges);
        for (const auto& sc : b.scales) {
          ref.scale.push_back(sc.scale);
          ref.separated.emplace_back(ref.nq * ref.nq, 0);
        }
        bands.push_back(std::move(ref));
      }
    }
    REQUIRE(tr.bands.size() == bands.size());
    for (std::size_t bi = 0; bi < bands.size(); ++bi) {
      auto& ref = bands[bi];
      const auto& b = tr.bands[bi];
      REQUIRE(b.anchor == ref.anchor);
      REQUIRE(b.num_quotient_vertices == ref.nq);
      REQUIRE(b.scales.size() == ref.scale.size());
      for (std::size_t si = 0; si < b.scales.size(); ++si) {
        const auto& sc = b.scales[si];
        REQUIRE(sc.scale == ref.scale[si]);
        REQUIRE(sc.center.size() == ref.nq);
        for (std::size_t x = 0; x < ref.nq; ++x) {
          for (std::size_t y = x + 1; y < ref.nq; ++y) {
            if (sc.center[x] != sc.center[y]) ++ref.separated[si][x * ref.nq + y];
          }
        }
      }
    }
  }

  SeparationResult out;
  for (const auto& ref : bands) {
    for (std::size_t x = 0; x < ref.nq; ++x) {
      for (std::size_t y = x + 1; y < ref.nq; ++y) {
        if (ref.qd[x][y] != WeightTraits<W>::infinity()) out.hypotheses += ref.scale.size();
      }
    }
  }
  if (out.hypotheses == 0) return out;
  out.radius = std::sqrt(std::log(static_cast<double>(out.hypotheses) / kDelta) /
                         (2.0 * static_cast<double>(seeds)));
  const double ht = harmonic(s.t);
  const double inv_h = std::ldexp(1.0, -static_cast<int>(eta));
  for (std::size_t bi = 0; bi < bands.size(); ++bi) {
    const auto& ref = bands[bi];
    for (std::size_t si = 0; si < ref.scale.size(); ++si) {
      const double delta_scale =
          std::ldexp(static_cast<double>(ref.anchor), static_cast<int>(ref.scale[si]));
      for (std::size_t x = 0; x < ref.nq; ++x) {
        for (std::size_t y = x + 1; y < ref.nq; ++y) {
          const W d = ref.qd[x][y];
          if (d == WeightTraits<W>::infinity()) continue;
          const double bound = ht * (static_cast<double>(d) / delta_scale + inv_h);
          const double freq =
              static_cast<double>(ref.separated[si][x * ref.nq + y]) / static_cast<double>(seeds);
          if (bound < 1.0) {
            ++out.nontrivial;
            out.worst_ratio = std::max(out.worst_ratio, freq / bound);
          }
          if (freq - bound > out.worst_excess) {
            out.worst_excess = freq - bound;
            std::ostringstream w;
            w << "band " << bi << " (anchor " << ref.anchor << ") scale " << ref.scale[si]
              << " pair (" << x << ", " << y << ") d " << d << ": frequency " << freq << ", bound "
              << bound;
            out.worst = w.str();
          }
          if (freq > bound + out.radius) ++out.violations;
        }
      }
    }
  }
  return out;
}

inline std::string describe(const std::string& label, const SeparationResult& r) {
  std::ostringstream s;
  s << label << ": " << r.hypotheses << " (band, scale, pair) triples, " << r.nontrivial
    << " with a bound below 1, Hoeffding radius " << r.radius << ", violations " << r.violations
    << ", largest frequency - bound " << r.worst_excess << ", largest frequency / bound "
    << r.worst_ratio << "\n    worst: " << r.worst;
  return s.str();
}

}  // namespace sssp::test::randomized
