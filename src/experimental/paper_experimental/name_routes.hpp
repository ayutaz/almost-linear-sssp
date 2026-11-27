#pragma once
// Mode B (paper-experimental) entry to NameRoutes, task P5-EXP-NAMEROUTES (arXiv 2611.34911v1,
// Sec 6.1 Algorithm 5 l.2-17, Definition 6.1, Lemma 6.8; Sec 6.2 l.1491 for H).
//
// The algorithm is the shared Mode A component paper::pivots::name_routes, called with the
// parameters that PaperParams derives for the instance's t (ADR 0003, Q-04: one implementation,
// parameterised). Nothing of Alg 5 is re-implemented here. The deviations it can carry:
//   D-9  H = least power of two >= 4tT^2 instead of >= (16t)^5 (Sec 6.2 l.1491; Q-42). Z >= d_S
//        holds for every H >= 1 (Lemma 6.7(1)); the proof of Lemma 6.7(2), and so of
//        Lemma 6.8(2), needs only H >= 4tT^2.
//   D-10 retry cap: after r consecutive rejected trials the loop stops and the residual edges
//        certify their own weight (as Alg 5 l.15). Not in the paper (Q-45). The shared
//        component records the FallbackRecord; this entry reports it as `retry_cap_reached`
//        and `fallback_edges` (the fallback share).
//   D-17 |S| = sample_factor * t, loop while |U| > residual_factor * t, P = 2^p_log2 instead of
//        16t, 32t and the P of Alg 5 l.3. Certificates stay valid (Lemma 6.8(1)); Lemma 6.8(2)
//        no longer applies.
// Under the paper's values the output, the counters, the counted operations, the fallbacks and
// the random words drawn are those of Mode A.
//
// The result carries the constants used next to the paper's values for the same t
// (NameRoutesSetting), so that dictionary sizes and trial counts of a run can be read against
// the paper's constants. This entry does not call PaperParams::record(); logging the parameter
// set of a run is the job of the Mode B driver (P5-EXP-THRESHOLDS).
//
// Correctness: Lemma 6.8(1) holds deterministically for every parameter value, so no input is
// unsupported here; invalid instances throw as in paper::pivots::name_routes.

#include <cstdint>
#include <optional>

#include "experimental/paper_experimental/params.hpp"
#include "paper/model/random.hpp"
#include "paper/pivots/name_routes.hpp"

namespace sssp::paper_experimental {

/// The NameRoutes constants for t pivots: the value used under the parameters and the paper's
/// value for the same t. For t < 2 the instance has no edges and Estimate never runs (Alg 5,
/// l.1422); H is then evaluated at t = 1 (it is only reported).
struct NameRoutesSetting {
  std::uint64_t t = 0;
  unsigned p_log2 = 0;                     ///< P used (Alg 5 l.3, or D-17)
  unsigned paper_p_log2 = 0;               ///< P of Alg 5 l.3
  unsigned h_log2 = 0;                     ///< H used (Sec 6.2 l.1491, or D-9)
  unsigned paper_h_log2 = 0;               ///< H of Sec 6.2 l.1491
  std::uint64_t sample_size = 0;           ///< |S| used (Alg 5 l.6, or D-17), capped at 2^64 - 1
  std::uint64_t paper_sample_size = 0;     ///< 16t
  std::uint64_t residual_limit = 0;        ///< loop while |U| > this (Alg 5 l.5, or D-17)
  std::uint64_t paper_residual_limit = 0;  ///< 32t
  std::optional<std::uint64_t> retry_cap;  ///< D-10; nullopt in the paper

  /// Every used value equals the paper's and there is no cap.
  [[nodiscard]] bool is_paper() const;
};

/// The setting of a run on t pivots under `params`. Throws std::length_error for t >= 2^32.
NameRoutesSetting name_routes_setting(const PaperParams& params, std::uint64_t t);

template <class W>
struct ExperimentalRouteNames {
  paper::pivots::RouteNames<W> names;  ///< output of paper::pivots::name_routes, unchanged
  NameRoutesSetting setting;
  /// D-10 stopped the loop; names.diag.fallbacks holds its FallbackRecord.
  bool retry_cap_reached = false;
  /// Edges certified by their own weight because of the cap (0 unless retry_cap_reached).
  std::uint64_t fallback_edges = 0;
};

/// Algorithm 5 under the Mode B parameters `params` (see the file comment). Throws as
/// paper::pivots::name_routes does for an invalid instance.
template <class W>
ExperimentalRouteNames<W> name_routes(const paper::pivots::PivotRoutingInstance<W>& instance,
                                      paper::model::WordRng& rng, const PaperParams& params);

}  // namespace sssp::paper_experimental
