#pragma once
// Sec 6 "Short Pivot Routes from Random Samples": Definition 6.1, Sec 6.1 "Sampling and
// Certificate Assignment", Algorithm 5 NameRoutes(G_P) and Lemma 6.8 (arXiv 2611.34911v1,
// txt l.1367-1635), Mode A.
//
// Definition 6.1: a pivot-routing instance is a loopless undirected multigraph on t pivots; edge
// e = xy carries the positive length a_e of a physical x-y walk. For a common factor P >= 1 a
// route certificate for e is an input weight lambda_e <= a_e such that a physical x-y walk of
// length <= P lambda_e exists. The dictionary lists the distinct certificate weights, sorted.
// Here every certificate walk is a concatenation of instance edges, so it is checked as
// dist_{G_P}(x, y) <= P lambda_e (check_route_names).
//
// Lemma 6.8 (proved deterministically for the output, Las Vegas for the time): every edge gets a
// certificate on every seed; a trial is accepted with probability >= 1/2 (measured in the tests);
// there is no retry limit in the paper (almost-sure termination).
//
// Parameters (NameRoutesParams). Mode A uses the defaults, which are the paper's values: sample
// 16t, residual threshold 32t, P = least power of two >= 8 * 256 * g * u^2, no retry cap. These
// only affect time and acceptance (Lemma 6.8(2)); certificates stay valid for any value. The
// correctness constants of Estimate (T, top scale, window) are not parameters (estimate.hpp).
//   * sample_factor / residual_factor / p_log2 other than the defaults: Mode B.
//   * retry_cap: deviation D-10 (research note Q-45), NOT in the paper and off by default. After
//     `retry_cap` consecutive rejected trials the loop stops and the residual edges certify their
//     own weight (as in Alg 5 l.15); this is recorded in Diagnostics.fallbacks.
//
// Randomness: only the caller's WordRng (uniform subsets by partial Fisher-Yates with the
// rejection method of Sec 2, and Estimate's draws). A fixed seed is reproducible everywhere.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/pivots/estimate.hpp"
#include "paper/primitives/search.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::pivots {

/// Definition 6.1: t pivots and M edges (x, y, a_e), loopless, a_e > 0.
template <class W>
struct PivotRoutingInstance {
  std::size_t t = 0;
  std::vector<Edge<W>> edges;
};

struct NameRoutesParams {
  std::uint64_t sample_factor = 16;    ///< |S| = sample_factor * t (Alg 5 l.6)
  std::uint64_t residual_factor = 32;  ///< loop while |U| > residual_factor * t (Alg 5 l.5)
  /// log2 P. nullopt: the paper's P (Alg 5 l.3). Other values: Mode B.
  std::optional<unsigned> p_log2;
  /// Deviation D-10 (not in the paper): stop after this many consecutive rejected trials.
  /// nullopt (default, Mode A): no cap.
  std::optional<std::uint64_t> retry_cap;
  EstimateParams estimate{};
};

/// Alg 5 l.2-3: u = ceil(log2(t + 2)), g = 4u + 3, P = 2^p_log2 = least power of two
/// >= 8 * 256 * g * u^2.
struct NameRoutesConstants {
  unsigned u = 0;
  unsigned g = 0;
  unsigned p_log2 = 0;
};
NameRoutesConstants name_routes_constants(std::size_t t);

template <class W>
struct RouteNames {
  std::vector<W> dictionary;        ///< distinct certificate weights, strictly increasing
  std::vector<std::uint32_t> name;  ///< per input edge: index into `dictionary`
  unsigned p_log2 = 0;              ///< P = 2^p_log2
  /// Counters "name_routes.*" (representatives, trials, accepted_rounds, rejected_trials,
  /// max_trials_per_round, max_consecutive_rejections, terminal_edges, dictionary_records,
  /// dictionary_size, retry_cap_hits) plus the summed "estimate.*" counters of every trial.
  Diagnostics diag;
};

/// Algorithm 5. Throws std::out_of_range for an endpoint >= t, std::invalid_argument for a loop
/// (discard loops before forming the instance, Sec 6), a nonpositive or non-finite weight, or
/// invalid parameters (sample_factor = 0, residual_factor < sample_factor, retry_cap = 0).
template <class W>
RouteNames<W> name_routes(const PivotRoutingInstance<W>& instance, model::WordRng& rng,
                          const NameRoutesParams& params = {});

/// Sec 6.1 median pruning: the least key k with feasible(k), for a predicate that is monotone
/// (false below some key, true from it on), or nullopt. Test a median (Lemma 3.1); on success
/// keep it as candidate and continue with the strictly smaller keys, otherwise with the strictly
/// larger keys (ties of the median are discarded either way). Each test at least halves the
/// keys, so at most floor(log2 n) + 1 tests run (counted into *tests if given).
template <class T, class Pred>
std::optional<T> median_prune(std::span<const T> keys, Pred feasible, OpCounters* ops = nullptr,
                              std::uint64_t* tests = nullptr) {
  std::vector<T> cur(keys.begin(), keys.end());
  std::optional<T> best;
  while (!cur.empty()) {
    const T m = select_bfprt(cur, (cur.size() - 1) / 2, std::ranges::less{}, ops);
    if (tests != nullptr) ++*tests;
    const bool ok = feasible(m);
    std::vector<T> next;
    for (const T& k : cur) {
      if (ops != nullptr) SSSP_COUNT(++ops->comparisons);
      if (ok ? k < m : m < k) next.push_back(k);
    }
    if (ok) best = m;
    cur = std::move(next);
  }
  return best;
}

/// Alg 5 l.6: moves a uniformly random k-element subset of `items` (in uniformly random order)
/// to items[0..k) by a partial Fisher-Yates shuffle: position i receives a uniform element of
/// i..n-1 (uniform_below, Sec 2). Throws std::invalid_argument if k > items.size().
template <class T>
void partial_shuffle(std::span<T> items, std::size_t k, model::WordRng& rng) {
  if (k > items.size()) throw std::invalid_argument("partial_shuffle: k exceeds the item count");
  for (std::size_t i = 0; i < k; ++i) {
    const auto j = i + static_cast<std::size_t>(model::uniform_below(rng, items.size() - i));
    using std::swap;
    swap(items[i], items[j]);
  }
}

/// Output invariants of Definition 6.1 / Lemma 6.8 against reference distances `dist` of G_P
/// (t x t, WeightTraits<W>::infinity() for unreachable pairs): one name per edge inside the
/// dictionary; the dictionary strictly increasing and made of input weights; lambda_e <= a_e;
/// dist(x, y) <= P lambda_e. Returns one message per violation (empty: all hold).
template <class W>
std::vector<std::string> check_route_names(const PivotRoutingInstance<W>& instance,
                                           const RouteNames<W>& names,
                                           const std::vector<std::vector<W>>& dist);

}  // namespace sssp::paper::pivots
