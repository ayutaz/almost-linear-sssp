#pragma once
// Mode B threshold settings, task P5-EXP-THRESHOLDS (arXiv 2611.34911v1; deviation D-7 with the
// limits decided in Q-03).
//
// The paper's size tests (Sec 2 l.286) are N > L_N^64 for the outer construction (Alg 8 l.4) and
// k > 2 L_k^4 for the local one (Alg 4 l.2). Below 2^64 the outer test never holds and the local
// one first holds at k = 209953, so neither construction runs on a feasible input in Mode A.
// Mode B lowers them (PaperParams, D-7). This file adds what is needed to choose lowered values
// safely and to state where they take effect:
//
//   * least_applying_size: the least k >= 1 at which a size test holds. The tests are not upward
//     closed (k > L_k^2 holds on 26..30, not on 31..36), so the value is found per range of
//     constant L_k: on {k : L_k = l} = [2^(l-1) - 1, 2^l - 2] the right side c * l^e is fixed
//     and the test is monotone in k.
//   * least_radial_local_size / least_unsupported_local_size: the least k at which Alg 4 takes
//     the radial branch, and the least k at which the size test selects it but the assembly's
//     preconditions fail (PaperParams::local_branch == kUnsupported; Sec 5.4 l.1206: k >= 32 and
//     h = b + 1 <= k for the paper's formulas with the grid, k >= 2 for D-8 / D-11). On each
//     range of constant L_k the separator b is fixed (a constant or L_k^e), so "supported" is
//     upward closed there as well; both values are exact for every k < 2^64.
//   * conservative_local_params(b): the Q-03 conservative variant, the paper's formulas with
//     exactly the sizes k >= max(32, b + 1) radial (k > max(31, b) * L_k^0), b constant.
//   * q03_threshold_settings: the named settings used for the Mode B differential checks
//     (P4-DIFF-MODE-B); each changes only D-7 entries and never meets an unsupported size.
//
// Model of computation: integer word arithmetic on sizes only, with products capped as in
// paper/model/params.hpp (Sec 2 l.304); no floating-point log or pow, no length is touched.

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "experimental/paper_experimental/params.hpp"

namespace sssp::paper_experimental {

/// The least k in [1, 2^64 - 1] with t.applies(k); nullopt if there is none (e.g. the paper's
/// outer test 1 * L^64).
[[nodiscard]] std::optional<std::uint64_t> least_applying_size(const SizeThreshold& t);

/// The least k >= 1 with p.local_branch(k) == LocalBranch::kRadial; nullopt if none below 2^64.
[[nodiscard]] std::optional<std::uint64_t> least_radial_local_size(const PaperParams& p);

/// The least k >= 1 with p.local_branch(k) == LocalBranch::kUnsupported; nullopt if Alg 4 can
/// run on every size under p. Mode B rejects the outer branch for a p with a value here.
[[nodiscard]] std::optional<std::uint64_t> least_unsupported_local_size(const PaperParams& p);

/// Q-03 conservative local setting: threshold max(31, b) * L^0 and separator b (constant), the
/// paper's alignment and grid; radial exactly for k >= max(32, b + 1). Throws
/// std::invalid_argument for b = 0.
[[nodiscard]] LocalParams conservative_local_params(std::uint64_t b);

/// A named Mode B setting with the sizes at which its constructions start.
struct ThresholdSetting {
  std::string name;
  PaperParams params;
  std::uint64_t least_outer_size = 0;                    ///< least_applying_size(outer threshold)
  std::optional<std::uint64_t> least_local_radial_size;  ///< least_radial_local_size(params)
};

/// The settings of Q-03 (all but the thresholds and separators are the paper's):
///   q03-outer-forced-b64        outer 0*L^0, b = 64; local as in the paper (atoms take Alg 4's
///                               base case: the minimal forcing set of note 04 Sec 3.4)
///   q03-both-forced-b64-l4      outer 0*L^0, b = 64; local conservative_local_params(4)
///   q03-both-forced-b1024-l31   outer 0*L^0, b = 1024; local conservative_local_params(31)
///   q03-scaled-L3-L2            outer 1*L^3 (from N = 1001; not on 1023..1331), b = L_N^2;
///                               local conservative_local_params(8)
[[nodiscard]] std::vector<ThresholdSetting> q03_threshold_settings();

}  // namespace sssp::paper_experimental
