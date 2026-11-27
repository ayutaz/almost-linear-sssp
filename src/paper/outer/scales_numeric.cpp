#include "paper/outer/scales_numeric.hpp"

#include <algorithm>
#include <array>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include "paper/local/saturating.hpp"
#include "sssp/graph.hpp"

namespace sssp::paper::outer {

namespace {

constexpr std::uint64_t kU64Max = std::numeric_limits<std::uint64_t>::max();

template <class W>
constexpr W kInf = WeightTraits<W>::infinity();

/// Positive and finite (false for NaN).
template <class W>
bool positive_finite(W x) {
  return x > W{0} && x != kInf<W>;
}

template <class W>
void require_positive_finite(W x, const char* what) {
  if (!positive_finite(x)) {
    throw std::invalid_argument(std::string("outer scales: ") + what +
                                " must be positive and finite");
  }
}

// ---- fixed-width nonnegative integers for h^2 (2h + 1) < 2^193: seven 32-bit limbs,
// little-endian. Integer word work on sizes (Sec 2 l.302-304), never on lengths.
constexpr std::size_t kLimbs = 7;
using Limbs = std::array<std::uint32_t, kLimbs>;

/// a * b; the caller guarantees that the product fits kLimbs limbs.
constexpr Limbs multiply(const Limbs& a, const Limbs& b) {
  Limbs r{};
  for (std::size_t i = 0; i < kLimbs; ++i) {
    std::uint64_t carry = 0;
    for (std::size_t j = 0; i + j < kLimbs; ++j) {
      // a[i] * b[j] + r[i + j] + carry <= (2^32 - 1)^2 + 2 (2^32 - 1) = 2^64 - 1: no wrap.
      const std::uint64_t cur = static_cast<std::uint64_t>(a[i]) * b[j] + r[i + j] + carry;
      r[i + j] = static_cast<std::uint32_t>(cur);
      carry = cur >> 32;
    }
  }
  return r;
}

/// Number of bits of a (0 for zero), by shifts of the top nonzero limb.
constexpr std::uint64_t bit_length(const Limbs& a) {
  for (std::size_t i = kLimbs; i-- > 0;) {
    if (a[i] == 0) continue;
    std::uint64_t bits = 32 * static_cast<std::uint64_t>(i);
    for (std::uint32_t top = a[i]; top != 0; top >>= 1) ++bits;
    return bits;
  }
  return 0;
}

/// min(a + b, 2^64 - 1).
constexpr std::uint64_t sat_add_u64(std::uint64_t a, std::uint64_t b) {
  return a > kU64Max - b ? kU64Max : a + b;
}

/// F^2 s and F a share this form: x * F^times.
template <class W>
W times_f_pow(const OuterScaleFactor& f, unsigned times, W x) {
  std::uint64_t k = 0;
  for (unsigned i = 0; i < times; ++i) k = sat_add_u64(k, f.log2);
  // A saturated doubling count is exact as well: any positive x saturates long before.
  return scale_pow2(x, k);
}

/// Index of the largest point <= w in points[lo, hi), or hi if there is none.
template <class W>
std::size_t floor_index(std::span<const W> points, std::size_t lo, std::size_t hi, W w) {
  const auto first = points.begin() + static_cast<std::ptrdiff_t>(lo);
  const auto last = points.begin() + static_cast<std::ptrdiff_t>(hi);
  const auto it = std::upper_bound(first, last, w);
  return it == first ? hi : static_cast<std::size_t>(it - points.begin()) - 1;
}

template <class W>
const OuterSupport<W>& support_of(const OuterGlobalGrid<W>& grid, std::uint32_t name) {
  if (name >= grid.supports.size()) throw std::invalid_argument("outer scales: unknown name");
  return grid.supports[name];
}

}  // namespace

OuterScaleFactor outer_scale_factor(unsigned log2_p, std::uint64_t h, unsigned log2_k) {
  if (h < 2) throw std::invalid_argument("outer_scale_factor: h must be at least 2");
  // l.2080: 2 P (2h+1) (8 K^2 h)^2 = 2^(7 + p + 4k) * h^2 (2h+1) for P = 2^p, K = 2^k. The least
  // power of two strictly above 2^c * M is 2^(c + bitlen(M)), also when M is a power of two.
  const Limbs hl{static_cast<std::uint32_t>(h), static_cast<std::uint32_t>(h >> 32)};
  // 2h + 1 < 2^65: three limbs.
  const Limbs twice_plus_one{static_cast<std::uint32_t>(h << 1) | 1U,
                             static_cast<std::uint32_t>(h >> 31),
                             static_cast<std::uint32_t>(h >> 63)};
  const Limbs m = multiply(multiply(hl, hl), twice_plus_one);  // < 2^193
  // 7 + p + 4k + bitlen < 2^3 + 2^32 + 2^34 + 2^8: no wrap.
  const std::uint64_t e = 7 + std::uint64_t{log2_p} + 4 * std::uint64_t{log2_k} + bit_length(m);
  return OuterScaleFactor{e, e < 64 ? std::uint64_t{1} << e : kU64Max};
}

template <class W>
W scale_pow2(W a, std::uint64_t k) {
  if (a == W{0}) return a;  // doubling zero never ends otherwise
  // Sec 2 l.297: k doublings; stops once the value has saturated.
  for (std::uint64_t i = 0; i < k && a != kInf<W>; ++i) a = sat::add(a, a);
  return a;
}

template <class W>
W outer_block_scale(std::span<const W> radii) {
  if (radii.empty()) throw std::invalid_argument("outer_block_scale: empty block");
  W s = radii.front();
  for (const W r : radii) {
    require_positive_finite(r, "a radius");
    if (r < s) s = r;
  }
  return s;
}

template <class W>
bool outer_block_aligned(const OuterScaleFactor& f, W s, W a) {
  require_positive_finite(s, "s_U");
  require_positive_finite(a, "a_U");
  return !(times_f_pow(f, 2, s) < a);  // l.2085: F^2 s >= a
}

template <class W>
bool outer_aligned_scale_in_range(const OuterScaleFactor& f, W s, W a) {
  return outer_block_aligned(f, s, a) && !(times_f_pow(f, 1, a) < s);  // l.2118
}

template <class W>
W outer_window_edge_bound(W s, std::uint64_t h, unsigned log2_k) {
  require_positive_finite(s, "s_U");
  if (h < 2) throw std::invalid_argument("outer_window_edge_bound: h must be at least 2");
  // l.2109: 4 K^2 h s = (h s) * 2^(2k + 2).
  return scale_pow2(sat::mul_int(s, h), 2 * static_cast<std::uint64_t>(log2_k) + 2);
}

template <class W>
std::size_t first_edge_not_below_half(std::span<const W> weights, W a) {
  require_positive_finite(a, "a_U");
  for (std::size_t i = 0; i < weights.size(); ++i) {
    if (!(sat::add(weights[i], weights[i]) < a)) return i;  // w < a / 2 as w + w < a
  }
  return weights.size();
}

template <class W>
bool outer_cutoff_valid(W s, W t) {
  return !(s < t) && s < sat::add(t, t);  // t <= s and s / 2 < t
}

template <class W>
OuterGlobalGrid<W> build_outer_global_grid(const OuterScaleFactor& f, std::span<const W> names,
                                           std::span<const OuterAlignedScale<W>> aligned) {
  OuterGlobalGrid<W> out;
  // l.2081: the support [a, F^3 a] of every name.
  out.supports.reserve(names.size());
  for (const W a : names) {
    require_positive_finite(a, "a name");
    out.supports.push_back(OuterSupport<W>{a, times_f_pow(f, 3, a), 0, 0});
  }
  // l.2086: an aligned block lowers its name's interval start to s_U.
  for (const auto& block : aligned) {
    if (block.name >= names.size()) {
      throw std::invalid_argument("build_outer_global_grid: block with an unknown name");
    }
    const W a = names[block.name];
    if (!outer_block_aligned(f, block.s, a)) {
      throw std::invalid_argument("build_outer_global_grid: block is not aligned");
    }
    if (times_f_pow(f, 1, a) < block.s) {
      throw std::logic_error("build_outer_global_grid: aligned block with s_U > F a_U");
    }
    auto& support = out.supports[block.name];
    if (block.s < support.lower) support.lower = block.s;
  }
  if (names.empty()) return out;

  // l.2088: sort the intervals by start and merge overlaps and gaps below a factor two.
  std::vector<std::uint32_t> order(names.size());
  std::iota(order.begin(), order.end(), std::uint32_t{0});
  std::sort(order.begin(), order.end(), [&](std::uint32_t x, std::uint32_t y) {
    if (out.supports[x].lower < out.supports[y].lower) return true;
    if (out.supports[y].lower < out.supports[x].lower) return false;
    return x < y;
  });
  std::vector<std::pair<W, W>> components;
  for (const std::uint32_t c : order) {
    const auto& support = out.supports[c];
    if (!components.empty() &&
        support.lower < sat::add(components.back().second, components.back().second)) {
      if (components.back().second < support.upper) components.back().second = support.upper;
    } else {
      components.emplace_back(support.lower, support.upper);
    }
  }
  // Lemma 8.4 item 2: a doubling grid per component, through the first point strictly above
  // its end ("one final doubling per merged grid"). A doubling that does not fit is the
  // saturated sentinel (D-13).
  for (const auto& [start, end] : components) {
    if (!out.points.empty() && out.points.back() == start) {
      // The previous sentinel is this start (a gap of exactly two): one shared point.
      out.component_begin.push_back(out.points.size() - 1);
    } else {
      out.component_begin.push_back(out.points.size());
      out.points.push_back(start);
    }
    for (W g = start; !(end < g) && g != kInf<W>;) {
      g = sat::add(g, g);
      out.points.push_back(g);
    }
  }
  // Pointers from each name into the table (l.2089).
  const std::span<const W> points(out.points);
  for (auto& support : out.supports) {
    support.lower_index = floor_index(points, 0, points.size(), support.lower);
    const auto above = std::upper_bound(points.begin(), points.end(), support.upper);
    // No point lies above a saturated end: the saturated sentinel (the last point) stands in.
    support.upper_index = above == points.end() ? points.size() - 1
                                                : static_cast<std::size_t>(above - points.begin());
  }
  return out;
}

template <class W>
std::size_t outer_grid_class(const OuterGlobalGrid<W>& grid, std::uint32_t name, W w) {
  const auto& support = support_of(grid, name);
  const std::span<const W> points(grid.points);
  if (!(w < points[support.upper_index])) return kNoGridClass;
  const std::size_t i = floor_index(points, support.lower_index, support.upper_index, w);
  return i == support.upper_index ? kNoGridClass : i;
}

template <class W>
W outer_aligned_cutoff(const OuterGlobalGrid<W>& grid, std::uint32_t name, W s) {
  const std::size_t i = outer_grid_class(grid, name, s);
  if (i == kNoGridClass) {
    throw std::logic_error("outer_aligned_cutoff: s_U outside the grid of its name");
  }
  const W t = grid.points[i];
  if (!outer_cutoff_valid(s, t)) {
    throw std::logic_error("outer_aligned_cutoff: s_U / 2 < t_U <= s_U fails");
  }
  return t;
}

template <class W>
W outer_separate_cutoff(std::span<const W> table, W s) {
  const std::size_t i = floor_index(table, 0, table.size(), s);
  if (i == table.size()) {
    throw std::logic_error("outer_separate_cutoff: the table starts above s_U");
  }
  if (!outer_cutoff_valid(s, table[i])) {
    throw std::logic_error("outer_separate_cutoff: s_U / 2 < t_U <= s_U fails");
  }
  return table[i];
}

// The macro argument W is a type in template argument lists (`OuterAlignedScale<W>>` reads as a
// shift to bugprone-macro-parentheses); a type cannot be parenthesised there.
// NOLINTBEGIN(bugprone-macro-parentheses)
#define SSSP_OUTER_SCALES_NUMERIC_INSTANTIATE(W)                                           \
  template W scale_pow2<W>(W, std::uint64_t);                                              \
  template W outer_block_scale<W>(std::span<const W>);                                     \
  template bool outer_block_aligned<W>(const OuterScaleFactor&, W, W);                     \
  template bool outer_aligned_scale_in_range<W>(const OuterScaleFactor&, W, W);            \
  template W outer_window_edge_bound<W>(W, std::uint64_t, unsigned);                       \
  template std::size_t first_edge_not_below_half<W>(std::span<const W>, W);                \
  template bool outer_cutoff_valid<W>(W, W);                                               \
  template OuterGlobalGrid<W> build_outer_global_grid<W>(                                  \
      const OuterScaleFactor&, std::span<const W>, std::span<const OuterAlignedScale<W>>); \
  template std::size_t outer_grid_class<W>(const OuterGlobalGrid<W>&, std::uint32_t, W);   \
  template W outer_aligned_cutoff<W>(const OuterGlobalGrid<W>&, std::uint32_t, W);         \
  template W outer_separate_cutoff<W>(std::span<const W>, W)
// NOLINTEND(bugprone-macro-parentheses)

SSSP_OUTER_SCALES_NUMERIC_INSTANTIATE(std::uint64_t);
SSSP_OUTER_SCALES_NUMERIC_INSTANTIATE(double);
#undef SSSP_OUTER_SCALES_NUMERIC_INSTANTIATE

}  // namespace sssp::paper::outer
