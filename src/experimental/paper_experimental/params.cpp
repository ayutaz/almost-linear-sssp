// Mode B parameter set, task P5-EXP-PARAMS. See params.hpp for the parameters, the fixed
// constants and the deviation IDs (arXiv 2611.34911v1; docs/research/thresholds.md §6-7).

#include "experimental/paper_experimental/params.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

#include "paper/local/assemble.hpp"
#include "paper/pivots/estimate.hpp"

namespace sssp::paper_experimental {

namespace {

constexpr std::uint64_t kWordMax = std::numeric_limits<std::uint64_t>::max();
constexpr unsigned kMaxSeparatorExponent = 64;  // the paper's largest separator is L^64
constexpr unsigned kMaxPLog2 = 256;             // the bound paper::pivots::name_routes accepts

// The defaults in params.hpp are the values of the shared Mode A components.
static_assert(LocalFixed{}.align_exponent == paper::kLocalAlignExponent);
static_assert(LocalFixed{}.support_exponent == paper::kLocalSupportExponent);
static_assert(OuterParams{}.sample_factor == paper::pivots::NameRoutesParams{}.sample_factor);
static_assert(OuterParams{}.residual_factor == paper::pivots::NameRoutesParams{}.residual_factor);
static_assert(sizeof(std::size_t) == sizeof(std::uint64_t));

/// Unbounded natural number (little-endian 32-bit limbs) for the exact size products that do
/// not fit a word: b = L^64, h = b + 1, (2h+1) h^2 and the like. Only integer word operations
/// (multiply, add, shift, compare) are used.
class Nat {
 public:
  explicit Nat(std::uint64_t x) {
    while (x != 0) {
      limbs_.push_back(static_cast<std::uint32_t>(x));
      x >>= 32U;
    }
  }

  [[nodiscard]] Nat operator*(const Nat& o) const {
    Nat r(0);
    if (limbs_.empty() || o.limbs_.empty()) return r;
    r.limbs_.assign(limbs_.size() + o.limbs_.size(), 0);
    for (std::size_t i = 0; i < limbs_.size(); ++i) {
      std::uint64_t carry = 0;
      for (std::size_t j = 0; j < o.limbs_.size(); ++j) {
        // (2^32 - 1)^2 + 2 (2^32 - 1) = 2^64 - 1: no overflow.
        const std::uint64_t cur = std::uint64_t{limbs_[i]} * o.limbs_[j] + r.limbs_[i + j] + carry;
        r.limbs_[i + j] = static_cast<std::uint32_t>(cur);
        carry = cur >> 32U;
      }
      r.limbs_[i + o.limbs_.size()] = static_cast<std::uint32_t>(carry);
    }
    r.trim();
    return r;
  }

  [[nodiscard]] Nat plus(std::uint32_t x) const {
    Nat r = *this;
    std::uint64_t carry = x;
    for (std::size_t i = 0; i < r.limbs_.size() && carry != 0; ++i) {
      const std::uint64_t cur = std::uint64_t{r.limbs_[i]} + carry;
      r.limbs_[i] = static_cast<std::uint32_t>(cur);
      carry = cur >> 32U;
    }
    if (carry != 0) r.limbs_.push_back(static_cast<std::uint32_t>(carry));
    return r;
  }

  /// this * 2^k.
  [[nodiscard]] Nat shifted(unsigned k) const {
    Nat r(0);
    if (limbs_.empty()) return r;
    r.limbs_.assign(k / 32, 0);
    const unsigned s = k % 32;
    std::uint32_t carry = 0;
    for (const std::uint32_t limb : limbs_) {
      r.limbs_.push_back(s == 0 ? limb : static_cast<std::uint32_t>((limb << s) | carry));
      carry = s == 0 ? 0 : limb >> (32 - s);
    }
    if (carry != 0) r.limbs_.push_back(carry);
    return r;
  }

  /// Number of binary digits (0 for 0), by shifting the top limb.
  [[nodiscard]] unsigned bit_width() const {
    if (limbs_.empty()) return 0;
    unsigned w = 0;
    for (std::uint32_t top = limbs_.back(); top != 0; top >>= 1U) ++w;
    return static_cast<unsigned>(32 * (limbs_.size() - 1)) + w;
  }

  [[nodiscard]] bool is_power_of_two() const {
    if (limbs_.empty()) return false;
    for (std::size_t i = 0; i + 1 < limbs_.size(); ++i) {
      if (limbs_[i] != 0) return false;
    }
    const std::uint32_t top = limbs_.back();
    return (top & (top - 1)) == 0;
  }

  /// Least k with 2^k >= this; this >= 1.
  [[nodiscard]] unsigned ceil_log2() const { return bit_width() - (is_power_of_two() ? 1 : 0); }

  [[nodiscard]] bool fits_word() const { return limbs_.size() <= 2; }

  /// min(this, 2^64 - 1).
  [[nodiscard]] std::uint64_t saturated() const {
    if (!fits_word()) return kWordMax;
    std::uint64_t r = 0;
    for (std::size_t i = limbs_.size(); i-- > 0;) r = (r << 32U) | limbs_[i];
    return r;
  }

 private:
  void trim() {
    while (!limbs_.empty() && limbs_.back() == 0) limbs_.pop_back();
  }
  std::vector<std::uint32_t> limbs_;
};

Nat nat_pow(std::uint64_t base, std::uint64_t e) {
  Nat r(1);
  const Nat b(base);
  for (std::uint64_t i = 0; i < e; ++i) r = r * b;
  return r;
}

/// The exact b of `rule` on n vertices (the rule's exponent is <= 64 after validation).
Nat separator_exact(const SeparatorRule& rule, std::uint64_t n) {
  if (rule.kind() == SeparatorRule::Kind::kConstant) return Nat(rule.parameter());
  return nat_pow(paper::model::L(n), rule.parameter());  // Sec 5.1 l.967 / Sec 7.1 l.1644
}

void require_fixed(const char* field, std::uint64_t value, std::uint64_t paper, const char* where) {
  if (value == paper) return;
  throw std::invalid_argument(std::string("paper_experimental: ") + field + " = " +
                              std::to_string(value) + ", but " + where + " fixes it at " +
                              std::to_string(paper) +
                              "; this correctness constant cannot be changed even in Mode B "
                              "(PAPER_MAPPING §2)");
}

void require_separator(const char* field, const SeparatorRule& rule) {
  if (rule.kind() == SeparatorRule::Kind::kConstant && rule.parameter() == 0) {
    throw std::invalid_argument(std::string("paper_experimental: ") + field +
                                " = 0; Def 5.1 needs b >= 1");
  }
  if (rule.kind() == SeparatorRule::Kind::kPowerOfL && rule.parameter() > kMaxSeparatorExponent) {
    throw std::invalid_argument(std::string("paper_experimental: ") + field + " = L^" +
                                std::to_string(rule.parameter()) +
                                "; exponents above 64 (the paper's outer b) are not accepted");
  }
}

// The paper's values below repeat the default member initializers of params.hpp on purpose: an
// edited default no longer matches the paper and makes even PaperParams{} throw.
void validate(const LocalParams& l) {
  const LocalFixed& f = l.fixed;
  require_fixed("local.shell_ratio_factor", f.shell_ratio_factor, 8, "Sec 5.2 l.1029");
  require_fixed("local.near_ratio_exponent", f.near_ratio_exponent, 2, "Sec 5.2 l.1029");
  require_fixed("local.macro_exponent", f.macro_exponent, 3, "Sec 5.2 l.1021");
  require_fixed("local.align_exponent", f.align_exponent, 10,
                "Sec 5.4 l.1206 (replaceable only as a whole by D-8)");
  require_fixed("local.support_exponent", f.support_exponent, 16,
                "Sec 5.4 l.1216 (replaceable only as a whole by D-8)");
  require_separator("local.separator", l.separator);
}

void validate(const OuterParams& o) {
  require_fixed("outer.shell_ratio", o.fixed.shell_ratio, 16, "Sec 8.1 l.1977");
  require_fixed("outer.near_ratio", o.fixed.near_ratio, 256, "Sec 8.1 l.1977");
  require_fixed("outer.scale_factor_coeff", o.fixed.scale_factor_coeff, 2, "Sec 8.2 l.2080");
  require_fixed("outer.support_cube", o.fixed.support_cube, 3, "Sec 8.2 l.2080");
  require_fixed("outer.band_factor", o.estimate.band_factor, 4, "Sec 6.2 l.1473");
  require_fixed("outer.window_exponent", o.estimate.window_exponent, 2, "Sec 6.2 l.1490");
  const CalendarFixed& c = o.calendars;
  require_fixed("outer.grid_octaves_below", c.octaves_below, 7, "Sec 7.3 l.1697");
  require_fixed("outer.grid_octaves_above", c.octaves_above, 2, "Sec 7.3 l.1697");
  require_fixed("outer.width_offset", c.width_offset, 5, "Sec 7.3 l.1698 (32q <= r < 64q)");
  require_fixed("outer.overlap_factor", c.overlap_factor, 4, "Sec 7.3, Lemma 7.6");
  require_fixed("outer.crossing_factor", c.crossing_factor, 16, "Sec 7.3, Lemma 7.8");
  require_fixed("outer.atom_factor", c.atom_factor, 2, "Sec 7.3, Lemma 7.5");
  require_separator("outer.separator", o.separator);
  if (o.sample_factor == 0) {
    throw std::invalid_argument("paper_experimental: outer.sample_factor = 0");
  }
  if (o.residual_factor < o.sample_factor) {
    throw std::invalid_argument("paper_experimental: outer.residual_factor < sample_factor");
  }
  if (o.retry_cap.has_value() && *o.retry_cap == 0) {
    throw std::invalid_argument("paper_experimental: outer.retry_cap = 0");
  }
  if (o.p_log2.has_value() && *o.p_log2 > kMaxPLog2) {
    throw std::invalid_argument("paper_experimental: outer.p_log2 > 256");
  }
}

void validate(const HierarchyFixed& h) {
  require_fixed("hierarchy.normalization_link_factor", h.normalization_link_factor, 4, "Sec 4.1");
}

ParamEntry entry(std::string key, std::string value, std::string paper, const char* id) {
  ParamEntry e{std::move(key), std::move(value), std::move(paper), ""};
  if (e.value != e.paper) e.deviation = id;
  return e;
}

std::string optional_string(const std::optional<std::uint64_t>& v, const char* none) {
  return v.has_value() ? std::to_string(*v) : std::string(none);
}

}  // namespace

std::string SizeThreshold::describe() const {
  return std::to_string(coeff) + "*L^" + std::to_string(exponent);
}

SeparatorValue SeparatorRule::evaluate(std::uint64_t n) const {
  if (kind_ == Kind::kConstant) return {parameter_, false};
  // min(L^e, 2^64 - 1) with capped products. L^e = 2^64 - 1 exactly is impossible (2^64 - 1 is
  // squarefree, so not a perfect power, and L <= 65), hence the cap is reached iff L^e >= 2^64.
  // Exponents above 64 (rejected by PaperParams) are evaluated as 65: L >= 2 then saturates as
  // well, and L = 1 (n = 0) gives 1 for every exponent.
  const auto e = static_cast<unsigned>(parameter_ > kMaxSeparatorExponent ? 65 : parameter_);
  const std::uint64_t v = paper::model::capped_pow(paper::model::L(n), e, kWordMax);
  return {v, v == kWordMax};
}

std::string SeparatorRule::describe() const {
  return kind_ == Kind::kConstant ? std::to_string(parameter_) : "L^" + std::to_string(parameter_);
}

PaperParams::PaperParams() : PaperParams(LocalParams{}, OuterParams{}) {}

PaperParams::PaperParams(const LocalParams& local, const OuterParams& outer,
                         const HierarchyFixed& hierarchy)
    : local_(local), outer_(outer), hierarchy_(hierarchy) {
  validate(local_);
  validate(outer_);
  validate(hierarchy_);
}

LocalBranch PaperParams::local_branch(std::uint64_t k) const {
  if (!local_applies(k)) return LocalBranch::kBaseCase;  // Alg 4 l.2
  if (local_.alignment == LocalAlignment::kLambda) {
    // D-8: research note 02 §4.3 claims the Lambda formulas for N >= 2 and b >= 1.
    return k >= 2 ? LocalBranch::kRadial : LocalBranch::kUnsupported;
  }
  // Sec 5.4 l.1206 as written: N >= 32 and h = b + 1 <= N (Q-03 conservative variant).
  const SeparatorValue b = local_.separator.evaluate(k);
  return !b.saturated && paper::local_assembly_applies(k, b.value) ? LocalBranch::kRadial
                                                                   : LocalBranch::kUnsupported;
}

LocalDerived PaperParams::local_derived(std::uint64_t k) const {
  const Nat b = separator_exact(local_.separator, k);
  const Nat h = b.plus(1);  // h = b + 1 (Sec 5.1)
  LocalDerived d;
  d.separator = b.saturated();
  d.hop_bound = h.saturated();
  d.saturated = !h.fits_word();
  // Sec 5.2 l.1029: K = least power of two >= 8h.
  d.shell_ratio_log2 = (Nat(local_.fixed.shell_ratio_factor) * h).ceil_log2();
  return d;
}

paper::LocalHierarchyParams PaperParams::local_hierarchy_params(std::uint64_t k) const {
  if (local_.alignment != LocalAlignment::kPaper) {
    throw std::logic_error(
        "paper_experimental: the Lambda alignment (D-8) is not provided by "
        "paper::local_hierarchy (P5-EXP-LOCAL)");
  }
  if (local_.assembly != LocalAssembly::kGrid) {
    throw std::logic_error(
        "paper_experimental: the exact-height assembly (D-11) is not provided by "
        "paper::local_hierarchy (P5-EXP-LOCAL)");
  }
  paper::LocalHierarchyParams p;
  p.threshold_coeff = local_.threshold.coeff;
  p.threshold_exponent = local_.threshold.exponent;
  // separator = 0 keeps the shared component on the paper's b = L_k^4 (Sec 5.1 l.967).
  p.separator = local_.separator == SeparatorRule::power_of_l(paper::model::kLocalExponent)
                    ? 0
                    : static_cast<std::size_t>(local_.separator.evaluate(k).value);
  return p;
}

OuterDerived PaperParams::outer_derived(std::uint64_t n, std::uint64_t t) const {
  if (t == 0) throw std::invalid_argument("paper_experimental: outer_derived needs t >= 1");
  const auto nr = paper::pivots::name_routes_constants(static_cast<std::size_t>(t));
  const auto ec = paper::pivots::estimate_constants(static_cast<std::size_t>(t));
  const Nat b = separator_exact(outer_.separator, n);
  const Nat h = b.plus(1);  // h = b + 1 (Sec 7.1 l.1644)
  OuterDerived d;
  d.separator = b.saturated();
  d.hop_bound = h.saturated();
  d.saturated = !h.fits_word();
  d.u = nr.u;                                    // Alg 5 l.2
  d.g = nr.g;                                    // Alg 5 l.2
  d.p_log2 = outer_.p_log2.value_or(nr.p_log2);  // Alg 5 l.3, or D-17
  d.t_log2 = ec.tau;                             // Sec 6.2 l.1473
  d.top_scale = ec.ell;                          // Sec 6.2 l.1490
  // Sec 6.2 l.1491, or D-9: least power of two >= 4tT^2 = 2^(2 + 2 tau) t, i.e. 2 + ell.
  d.h_log2 = outer_.h_rule == HRule::kPaper ? ec.eta : ec.ell + 2;
  // Sec 8.2 l.2080: F = least power of two > 2P(2h+1)(8K^2h)^2 with K = 16; the least f with
  // 2^f > X is the bit width of X. P is the value actually used, so F follows D-17.
  const std::uint64_t k = outer_.fixed.shell_ratio;
  const Nat eight_k2_h = Nat(8 * k * k) * h;
  const Nat x = Nat(outer_.fixed.scale_factor_coeff).shifted(d.p_log2) * h.shifted(1).plus(1) *
                eight_k2_h * eight_k2_h;
  d.f_log2 = x.bit_width();
  return d;
}

paper::pivots::NameRoutesParams PaperParams::name_routes_params(std::uint64_t t) const {
  paper::pivots::NameRoutesParams p;
  p.sample_factor = outer_.sample_factor;      // Alg 5 l.6, or D-17
  p.residual_factor = outer_.residual_factor;  // Alg 5 l.5, or D-17
  p.p_log2 = outer_.p_log2;                    // Alg 5 l.3, or D-17
  p.retry_cap = outer_.retry_cap;              // D-10 (not in the paper)
  if (outer_.h_rule == HRule::kProofMinimum) {
    if (t == 0) throw std::invalid_argument("paper_experimental: D-9 needs t >= 1");
    // D-9: H = least power of two >= 4tT^2 (see outer_derived).
    p.estimate.h_log2 = paper::pivots::estimate_constants(static_cast<std::size_t>(t)).ell + 2;
  }
  return p;
}

std::vector<ParamEntry> PaperParams::parameter_log() const {
  const LocalParams paper_local;
  const OuterParams paper_outer;
  std::vector<ParamEntry> log;
  log.push_back(entry("local.threshold", local_.threshold.describe(),
                      paper_local.threshold.describe(), "D-7"));
  log.push_back(entry("local.separator", local_.separator.describe(),
                      paper_local.separator.describe(), "D-7"));
  log.push_back(entry("local.alignment",
                      local_.alignment == LocalAlignment::kPaper ? "paper" : "lambda", "paper",
                      "D-8"));
  log.push_back(entry("local.assembly",
                      local_.assembly == LocalAssembly::kGrid ? "grid" : "exact-height", "grid",
                      "D-11"));
  log.push_back(entry("outer.threshold", outer_.threshold.describe(),
                      paper_outer.threshold.describe(), "D-7"));
  log.push_back(entry("outer.separator", outer_.separator.describe(),
                      paper_outer.separator.describe(), "D-7"));
  log.push_back(entry("outer.sample_factor", std::to_string(outer_.sample_factor),
                      std::to_string(paper_outer.sample_factor), "D-17"));
  log.push_back(entry("outer.residual_factor", std::to_string(outer_.residual_factor),
                      std::to_string(paper_outer.residual_factor), "D-17"));
  log.push_back(entry("outer.p_log2",
                      outer_.p_log2.has_value() ? std::to_string(*outer_.p_log2) : "paper", "paper",
                      "D-17"));
  log.push_back(
      entry("outer.h_rule", outer_.h_rule == HRule::kPaper ? "paper" : "4tT^2", "paper", "D-9"));
  log.push_back(entry("outer.retry_cap", optional_string(outer_.retry_cap, "off"), "off", "D-10"));
  return log;
}

std::vector<ParamEntry> PaperParams::deviations() const {
  std::vector<ParamEntry> out;
  for (auto& e : parameter_log()) {
    if (!e.deviation.empty()) out.push_back(std::move(e));
  }
  return out;
}

void PaperParams::record(Diagnostics& diag) const {
  std::uint64_t deviating = 0;
  for (const ParamEntry& e : parameter_log()) {
    diag.counters["params." + e.key + "=" + e.value] = 1;
    if (!e.deviation.empty()) {
      diag.counters["params.deviation." + e.deviation + "." + e.key] = 1;
      ++deviating;
    }
  }
  diag.counters["params.deviations"] = deviating;
}

}  // namespace sssp::paper_experimental
