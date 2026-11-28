#include "experimental/paper_experimental/thresholds.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/model/params.hpp"

namespace sssp::paper_experimental {

namespace {

constexpr std::uint64_t kWordMax = std::numeric_limits<std::uint64_t>::max();

/// The sizes k >= 1 with L_k = l (L_k = floor(log2(k + 1)) + 1, paper/model/params.hpp).
struct SizeRange {
  std::uint64_t lo;
  std::uint64_t hi;
};

/// [2^(l-1) - 1, 2^l - 2] for l = 2..64, then {2^64 - 1} (L = 65); l = 1 holds only k = 0.
constexpr std::array<SizeRange, 64> kRanges = [] {
  std::array<SizeRange, 64> r{};
  for (unsigned l = 2; l <= 64; ++l) {
    const std::uint64_t lo = (std::uint64_t{1} << (l - 1)) - 1;
    const std::uint64_t hi = l == 64 ? kWordMax - 1 : (std::uint64_t{1} << l) - 2;
    r[l - 2] = {lo, hi};
  }
  r[63] = {kWordMax, kWordMax};
  return r;
}();

/// The least k in `range` with pred(k), for a pred that is upward closed on the range.
template <class Pred>
std::optional<std::uint64_t> least_in(const SizeRange& range, Pred pred) {
  if (!pred(range.hi)) return std::nullopt;
  std::uint64_t lo = range.lo;
  std::uint64_t hi = range.hi;  // pred(hi) holds
  while (lo < hi) {
    const std::uint64_t mid = lo + (hi - lo) / 2;
    if (pred(mid)) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  return lo;
}

/// The least k over all ranges with pred(k), for a pred upward closed on every range.
template <class Pred>
std::optional<std::uint64_t> least_over_ranges(Pred pred) {
  for (const SizeRange& range : kRanges) {
    if (const auto k = least_in(range, pred)) return k;
  }
  return std::nullopt;
}

}  // namespace

std::optional<std::uint64_t> least_applying_size(const SizeThreshold& t) {
  // On a range of constant L_k the test is k > c * l^e with a fixed right side: upward closed.
  return least_over_ranges([&](std::uint64_t k) { return t.applies(k); });
}

std::optional<std::uint64_t> least_radial_local_size(const PaperParams& p) {
  // Radial = the size test holds and the assembly's preconditions hold (local_branch). On a
  // range of constant L_k both are upward closed: the size test as above, and the
  // preconditions (k >= 32 and b < k with a fixed b, or k >= 2) as well.
  return least_over_ranges(
      [&](std::uint64_t k) { return p.local_branch(k) == LocalBranch::kRadial; });
}

std::optional<std::uint64_t> least_unsupported_local_size(const PaperParams& p) {
  // On a range of constant L_k the sizes where the test holds are [a, hi] and the unsupported
  // ones among them are downward closed in it (the preconditions are upward closed), so an
  // unsupported size exists in the range iff a itself is unsupported.
  for (const SizeRange& range : kRanges) {
    const auto a = least_in(range, [&](std::uint64_t k) { return p.local_applies(k); });
    if (a && p.local_branch(*a) == LocalBranch::kUnsupported) return a;
  }
  return std::nullopt;
}

LocalParams conservative_local_params(std::uint64_t b) {
  if (b == 0) {
    throw std::invalid_argument("paper_experimental: conservative_local_params needs b >= 1");
  }
  // Q-03 (decided): keep the formulas of Sec 5.4 as written and take the radial branch exactly
  // for k >= max(32, b + 1), i.e. k > max(31, b) * L_k^0 (Sec 5.4 l.1206 then holds).
  LocalParams l;
  l.threshold = SizeThreshold{b > 31 ? b : 31, 0};
  l.separator = SeparatorRule::constant(b);
  return l;
}

std::vector<ThresholdSetting> q03_threshold_settings() {
  const auto forced_outer = [](SeparatorRule b) {
    OuterParams o;
    o.threshold = SizeThreshold{0, 0};  // Alg 8 l.4 holds for every N >= 1 (D-7)
    o.separator = b;
    return o;
  };
  OuterParams scaled;
  scaled.threshold = SizeThreshold{1, 3};           // N > L_N^3 (D-7; thresholds.md §7.2)
  scaled.separator = SeparatorRule::power_of_l(2);  // b = L_N^2 (D-7)

  std::vector<ThresholdSetting> out;
  const auto add = [&](const char* name, const LocalParams& l, const OuterParams& o) {
    ThresholdSetting s{name, PaperParams{l, o}, 0, std::nullopt};
    // Every setting has an outer size test that holds somewhere below 2^64.
    s.least_outer_size = least_applying_size(o.threshold).value();
    s.least_local_radial_size = least_radial_local_size(s.params);
    out.push_back(std::move(s));
  };
  add("q03-outer-forced-b64", LocalParams{}, forced_outer(SeparatorRule::constant(64)));
  add("q03-both-forced-b64-l4", conservative_local_params(4),
      forced_outer(SeparatorRule::constant(64)));
  add("q03-both-forced-b1024-l31", conservative_local_params(31),
      forced_outer(SeparatorRule::constant(1024)));
  add("q03-scaled-L3-L2", conservative_local_params(8), scaled);
  return out;
}

}  // namespace sssp::paper_experimental
