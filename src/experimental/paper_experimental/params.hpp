#pragma once
// Mode B (paper-experimental) parameter set, task P5-EXP-PARAMS (arXiv 2611.34911v1).
//
// Mode B runs the paper's constructions with the time-only constants turned into parameters so
// that the outer and local constructions trigger on feasible sizes (research note Q-03,
// docs/research/thresholds.md §6-7). Mode A (src/paper/) is unchanged and keeps the constants of
// src/paper/model/params.hpp; the shared components are called with the values derived here
// (ADR 0003, Q-04: one implementation, parameterised).
//
// Tunable (each value other than the paper's is a listed deviation):
//   local.threshold   k > c * L_k^e, paper c = 2, e = 4 (Sec 2 l.286, Alg 4 l.2)      D-7
//   local.separator   b, paper L_k^4 (Sec 5.1 l.967)                                   D-7
//   local.alignment   Sec 5.4 c = 10 / N^16 (paper) or the Lambda formulas             D-8
//   local.assembly    Sec 5.4 grid (paper) or the exact-height assembly                D-11
//   outer.threshold   N > c * L_N^e, paper c = 1, e = 64 (Sec 2 l.286, Alg 8 l.4)      D-7
//   outer.separator   b, paper L_N^64 (Sec 7.1 l.1644)                                 D-7
//   outer.sample_factor / residual_factor / p_log2   16t, 32t, P (Alg 5 l.3-6)         D-17
//   outer.h_rule      H = pow2 >= (16t)^5 (Sec 6.2 l.1491) or pow2 >= 4tT^2            D-9
//   outer.retry_cap   none in the paper; NameRoutes stops after r rejections           D-10
//
// Fixed (PAPER_MAPPING §2 "Cannot be changed even in Mode B"; Sec 5.4 c and N^16 only as a
// whole through D-8). They are fields so that a caller cannot silently assume another value:
// PaperParams rejects any value other than the paper's with std::invalid_argument.
//
// Derived constants are recomputed from b (and t) on every query, never stored: h = b + 1, the
// local K = least power of two >= 8h (Sec 5.2 l.1029), u, g, P (Alg 5 l.2-3), T, the top scale,
// H (Sec 6.2 l.1473-1491) and F = least power of two > 2P(2h+1)(8K^2h)^2 with the outer K = 16
// (Sec 8.2 l.2080). The size arithmetic is exact integer arithmetic (multiword where b = L^64
// does not fit a word); no floating-point log or pow is used.
//
// Logging: parameter_log() lists every tunable parameter with its value, the paper's value and
// the deviation ID if they differ; record() writes it into Diagnostics.counters as
// "params.<key>=<value>" = 1, "params.deviation.<D-id>.<key>" = 1 and "params.deviations" = the
// number of deviating entries. Like FallbackRecords (ADR 0006) this is part of the result and is
// written also when SSSP_DIAGNOSTICS = 0, so that timings from the release build carry the
// parameter set they were measured with. A parameter choice is not a fallback: record() pushes
// no FallbackRecord.

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "paper/local/local_hierarchy.hpp"
#include "paper/model/params.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/result.hpp"

namespace sssp::paper_experimental {

/// Size test "N > coeff * L_N^exponent" (Sec 2 l.284-286), evaluated with products capped at N.
/// coeff = 0 makes the construction apply to every N >= 1.
struct SizeThreshold {
  std::uint64_t coeff = 0;
  unsigned exponent = 0;

  [[nodiscard]] bool applies(std::uint64_t n) const {
    return paper::model::exceeds_scaled_pow(n, coeff, paper::model::L(n), exponent);
  }
  [[nodiscard]] std::string describe() const;  ///< "<coeff>*L^<exponent>"
  friend bool operator==(const SizeThreshold&, const SizeThreshold&) = default;
};

/// b evaluated for one instance size; `saturated` iff the exact b is >= 2^64 (value = 2^64 - 1).
struct SeparatorValue {
  std::uint64_t value = 0;
  bool saturated = false;
};

/// Separator size b of the owner decomposition (Def 5.1): b = L_N^e or a constant b >= 1.
class SeparatorRule {
 public:
  enum class Kind : std::uint8_t { kPowerOfL, kConstant };

  /// b = L_N^e (paper: e = 4 local, e = 64 outer). e <= 64 is required by PaperParams.
  static constexpr SeparatorRule power_of_l(unsigned e) { return {Kind::kPowerOfL, e}; }
  /// b fixed for every N (thresholds.md §7.2: constants for tests). b >= 1 is required.
  static constexpr SeparatorRule constant(std::uint64_t b) { return {Kind::kConstant, b}; }

  [[nodiscard]] constexpr Kind kind() const { return kind_; }
  [[nodiscard]] constexpr std::uint64_t parameter() const { return parameter_; }
  [[nodiscard]] SeparatorValue evaluate(std::uint64_t n) const;
  [[nodiscard]] std::string describe() const;  ///< "L^<e>" or "<b>"
  friend constexpr bool operator==(const SeparatorRule&, const SeparatorRule&) = default;

 private:
  constexpr SeparatorRule(Kind kind, std::uint64_t parameter)
      : kind_(kind), parameter_(parameter) {}
  Kind kind_;
  std::uint64_t parameter_;
};

/// Sec 5.4 alignment test and support: the paper's s_U >= sigma_U / N^10 with support
/// [., N^16 sigma] (valid for k >= max(32, b + 1), Q-03 conservative variant), or deviation D-8:
/// Lambda = 2N * 3hK^2 and support end 2(3hK^2)^2 N^3 sigma (claimed for N >= 2, b >= 1).
enum class LocalAlignment : std::uint8_t { kPaper, kLambda };

/// Sec 5.4 grid assembly (paper), or the exact-height assembly of deviation D-11.
enum class LocalAssembly : std::uint8_t { kGrid, kExactHeight };

/// H of Sec 6.2 l.1491: the paper's least power of two >= (16t)^5, or D-9's least power of two
/// >= 4tT^2 (all the proof of Lemma 6.7 needs, Q-42).
enum class HRule : std::uint8_t { kPaper, kProofMinimum };

/// Local correctness constants (fixed).
struct LocalFixed {
  unsigned shell_ratio_factor = 8;   ///< Sec 5.2 l.1029: K = least power of two >= 8h
  unsigned near_ratio_exponent = 2;  ///< Sec 5.2: near means ratio <= K^2
  unsigned macro_exponent = 3;       ///< Sec 5.2 l.1021: macro iff w <= N^3 sigma
  unsigned align_exponent = 10;      ///< Sec 5.4 l.1206: c in s_U >= sigma_U / N^c
  unsigned support_exponent = 16;    ///< Sec 5.4 l.1216: support end N^16 sigma
  friend bool operator==(const LocalFixed&, const LocalFixed&) = default;
};

/// Outer correctness constants of Sec 8 (fixed).
struct OuterFixed {
  unsigned shell_ratio = 16;        ///< Sec 8.1 l.1977: K = 16
  unsigned near_ratio = 256;        ///< Sec 8.1 l.1977: near means ratio <= K^2 = 256
  unsigned scale_factor_coeff = 2;  ///< Sec 8.2 l.2080: F > 2 P (2h+1) (8K^2h)^2
  unsigned support_cube = 3;        ///< Sec 8.2: support [a, F^3 a]
  friend bool operator==(const OuterFixed&, const OuterFixed&) = default;
};

/// Estimate correctness constants of Sec 6.2 (fixed).
struct EstimateFixed {
  unsigned band_factor = 4;      ///< l.1473: T = least power of two >= 4t
  unsigned window_exponent = 2;  ///< l.1490: top scale 2^l >= tT^2, window [a/T^2, T^2 a]
  friend bool operator==(const EstimateFixed&, const EstimateFixed&) = default;
};

/// Calendar correctness constants of Sec 7.3 (fixed).
struct CalendarFixed {
  unsigned octaves_below = 7;     ///< l.1697: grid from 7 octaves below the origin
  unsigned octaves_above = 2;     ///< l.1697: to 2 octaves above the maximum R
  unsigned width_offset = 5;      ///< l.1698: 32q <= r < 64q (2^5)
  unsigned overlap_factor = 4;    ///< Lemma 7.6: piece overlap [r/4, 4R]
  unsigned crossing_factor = 16;  ///< Lemma 7.8: piece crossing r_v > 16 R_u
  unsigned atom_factor = 2;       ///< Lemma 7.5: atom edges w' < 2 max(q_u, q_v)
  friend bool operator==(const CalendarFixed&, const CalendarFixed&) = default;
};

/// Sec 4.1 normalization (shared by every hierarchy, fixed).
struct HierarchyFixed {
  unsigned normalization_link_factor = 4;  ///< Sec 4.1: 4 S(C) >= h(P)
  friend bool operator==(const HierarchyFixed&, const HierarchyFixed&) = default;
};

/// Local construction (Alg 4). Defaults are the paper's values.
struct LocalParams {
  SizeThreshold threshold{paper::model::kLocalCoeff, paper::model::kLocalExponent};
  SeparatorRule separator = SeparatorRule::power_of_l(paper::model::kLocalExponent);
  LocalAlignment alignment = LocalAlignment::kPaper;
  LocalAssembly assembly = LocalAssembly::kGrid;
  LocalFixed fixed{};
};

/// Outer construction (Alg 8 l.8-17) with NameRoutes (Alg 5) and Estimate (Alg 6). Defaults are
/// the paper's values.
struct OuterParams {
  SizeThreshold threshold{1, paper::model::kOuterExponent};
  SeparatorRule separator = SeparatorRule::power_of_l(paper::model::kOuterExponent);
  std::uint64_t sample_factor = 16;    ///< Alg 5 l.6: |S| = 16t
  std::uint64_t residual_factor = 32;  ///< Alg 5 l.5: loop while |U| > 32t
  std::optional<unsigned> p_log2;      ///< Alg 5 l.3: nullopt = the paper's P
  HRule h_rule = HRule::kPaper;
  std::optional<std::uint64_t> retry_cap;  ///< nullopt = no cap (the paper)
  OuterFixed fixed{};
  EstimateFixed estimate{};
  CalendarFixed calendars{};
};

/// Local derived constants for one instance of size k. `saturated` iff b or h = b + 1 does not
/// fit a word (the word fields then hold 2^64 - 1); shell_ratio_log2 is exact in every case.
struct LocalDerived {
  std::uint64_t separator = 0;  ///< b
  std::uint64_t hop_bound = 0;  ///< h = b + 1
  bool saturated = false;
  unsigned shell_ratio_log2 = 0;  ///< K = 2^shell_ratio_log2 >= 8h
};

/// Outer derived constants for N vertices and t pivots. Saturation as in LocalDerived; f_log2
/// is exact in every case.
struct OuterDerived {
  std::uint64_t separator = 0;  ///< b
  std::uint64_t hop_bound = 0;  ///< h = b + 1
  bool saturated = false;
  unsigned u = 0;          ///< ceil(log2(t + 2))        (Alg 5 l.2)
  unsigned g = 0;          ///< 4u + 3                   (Alg 5 l.2)
  unsigned p_log2 = 0;     ///< P actually used          (Alg 5 l.3 or D-17)
  unsigned t_log2 = 0;     ///< T = least pow2 >= 4t     (Sec 6.2 l.1473)
  unsigned top_scale = 0;  ///< least l, 2^l >= tT^2     (Sec 6.2 l.1490)
  unsigned h_log2 = 0;     ///< H                        (Sec 6.2 l.1491 or D-9)
  unsigned f_log2 = 0;     ///< F                        (Sec 8.2 l.2080)
};

/// Which branch Alg 4 takes on k vertices under these parameters. kUnsupported: the size test
/// selects the radial construction but its preconditions do not hold (Sec 5.4 l.1206 for the
/// paper's formulas, N >= 2 for D-8). The caller must treat such an input as a failure; it may
/// not substitute another construction.
enum class LocalBranch : std::uint8_t { kBaseCase, kRadial, kUnsupported };

/// One logged parameter. `deviation` is empty iff value == paper.
struct ParamEntry {
  std::string key;
  std::string value;
  std::string paper;
  std::string deviation;  ///< "D-7", ..., or empty
};

/// A validated Mode B parameter set. The constructor throws std::invalid_argument if a fixed
/// constant differs from the paper (the message names the field and the paper location) or if a
/// tunable value is outside its domain (b = 0, a separator exponent above 64, sample_factor = 0,
/// residual_factor < sample_factor, retry_cap = 0, p_log2 > 256).
class PaperParams {
 public:
  PaperParams();  ///< the paper's values
  PaperParams(const LocalParams& local, const OuterParams& outer,
              const HierarchyFixed& hierarchy = {});

  [[nodiscard]] const LocalParams& local() const { return local_; }
  [[nodiscard]] const OuterParams& outer() const { return outer_; }
  [[nodiscard]] const HierarchyFixed& hierarchy() const { return hierarchy_; }

  /// Alg 4 l.2 size test on k.
  [[nodiscard]] bool local_applies(std::uint64_t k) const { return local_.threshold.applies(k); }
  /// Alg 8 l.4 size test on N.
  [[nodiscard]] bool outer_applies(std::uint64_t n) const { return outer_.threshold.applies(n); }

  [[nodiscard]] LocalBranch local_branch(std::uint64_t k) const;
  [[nodiscard]] LocalDerived local_derived(std::uint64_t k) const;
  /// Parameters for the shared paper::local_hierarchy on k vertices (b = 0 keeps the paper's
  /// L_k^4 path). Throws std::logic_error for D-8 or D-11, which the shared component does not
  /// provide (P5-EXP-LOCAL).
  [[nodiscard]] paper::LocalHierarchyParams local_hierarchy_params(std::uint64_t k) const;

  /// Throws std::invalid_argument for t = 0 and std::length_error for t >= 2^32.
  [[nodiscard]] OuterDerived outer_derived(std::uint64_t n, std::uint64_t t) const;
  /// Parameters for paper::pivots::name_routes on t pivots (D-9 needs t for H). Throws
  /// std::invalid_argument for t = 0 under D-9.
  [[nodiscard]] paper::pivots::NameRoutesParams name_routes_params(std::uint64_t t) const;

  [[nodiscard]] std::vector<ParamEntry> parameter_log() const;
  [[nodiscard]] std::vector<ParamEntry> deviations() const;
  [[nodiscard]] bool is_paper() const { return deviations().empty(); }
  /// Writes the parameter log into diag.counters (see the file comment). Idempotent.
  void record(Diagnostics& diag) const;

 private:
  LocalParams local_;
  OuterParams outer_;
  HierarchyFixed hierarchy_;
};

}  // namespace sssp::paper_experimental
