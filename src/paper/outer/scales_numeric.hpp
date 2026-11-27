#pragma once
// Numerical scales of the outer construction (arXiv 2611.34911v1, Sec 8.2, txt l.2076-2104, and
// Lemma 8.4, l.2105-2142), Mode A: the parts that depend only on numbers, not on the outer
// state representation of Sec 8.1. The caller (the Sec 8.2 assembly) supplies radii, names,
// block-MST weights and grid tables as plain sequences.
//
//   * F (l.2080): the least power of two strictly above 2P(2h+1)(8K^2h)^2, chosen once from
//     integer parameters. P = 2^p and K = 2^k are powers of two (Alg 5 l.3; Sec 8.1: K = 16),
//     so the product is 2^(7+p+4k) h^2 (2h+1) and F = 2^(7+p+4k+bitlen(h^2 (2h+1))). F is held
//     as (value, exponent) (Q-06): the exponent is exact, the 64-bit value saturates. h is a
//     64-bit word here; the paper's own h = L_N^64 + 1 never fits one (589 bits at the first N
//     where the outer construction runs), so Mode A parameters cannot reach this function.
//   * Block scale (Sec 8.1 l.1999, Q-35): s_U = min of the original radii r_v over the block.
//   * Alignment (l.2085): U is aligned iff F^2 s_U >= a_U. An unaligned block has every MST edge
//     below a_U / 2 (l.2085-2086; Lemma 8.4 item 1: w <= 4 K^2 h s_U < a_U / 2). An aligned
//     block has a_U / F^2 <= s_U <= F a_U (l.2118).
//   * Supports and global grid (l.2081-2089, Lemma 8.4 item 2): [a, F^3 a] per name, the start
//     lowered to the least aligned s_U of that name; intervals that overlap or leave a gap
//     below a factor two are merged; each merged component gets a doubling grid from its start
//     through the first point strictly above its end (the lookup sentinel). Each support stores
//     the index of the grid floor of its start and of the first point above its end.
//   * Cutoff (l.2090-2092, Q-36): t_U is the grid floor of s_U, on the name's global grid for an
//     aligned block and on a separate table otherwise; s_U / 2 < t_U <= s_U in both cases.
//
// Model of computation (Sec 2 l.292-300): lengths are only copied, compared and added;
// multiplication by 2^j is j doublings and by an integer its binary expansion. There is no
// halving: x < a / 2 is tested as x + x < a, and s / 2 < t as s < t + t.
//
// Numerics [D-13, adopted for Sec 5 under Q-38 and used here the same way]: F^2 s, F a, F^3 a,
// 4 K^2 h s and the grid doublings saturate at WeightTraits<W>::infinity()
// (paper/local/saturating.hpp). They are only compared with finite lengths, and a product that
// does not fit exceeds every finite length of the type, so each comparison is exact. A
// saturated support end merges the support with every later one; its grid ends with the
// saturated sentinel infinity(). Only grid points move; every class still has upper <= 2 lower.
//
// double (Q-55 option A): doublings are exact, so every decision here is exact for any finite
// positive double, except outer_window_edge_bound, whose integer multiple rounds once the
// product needs more than 53 significant bits. Correctness is promised for integer and dyadic
// weights.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace sssp::paper::outer {

/// F = 2^log2 (Sec 8.2 l.2080) as an exact exponent with a saturated 64-bit value.
struct OuterScaleFactor {
  std::uint64_t log2 = 0;   ///< exact: F = 2^log2
  std::uint64_t value = 0;  ///< min(F, 2^64 - 1)

  /// True iff `value` equals F (F < 2^64).
  [[nodiscard]] constexpr bool value_exact() const noexcept { return log2 < 64; }
};

/// The least power of two strictly above 2 P (2h+1) (8 K^2 h)^2 with P = 2^log2_p and
/// K = 2^log2_k (Sec 8.1: K = 16, log2_k = 4). Integer arithmetic only; h^2 (2h+1) is formed
/// exactly in 224 bits. h must be the exact hop bound b + 1: a value that does not fit 64 bits
/// cannot be passed, and a saturated stand-in gives a wrong F (the caller must reject such
/// parameters). Throws std::invalid_argument unless h >= 2 (b >= 1, Sec 7 l.1644).
OuterScaleFactor outer_scale_factor(unsigned log2_p, std::uint64_t h, unsigned log2_k);

/// a * 2^k by doublings, saturating at infinity. Zero stays zero for every k.
template <class W>
W scale_pow2(W a, std::uint64_t k);

/// s_U = min r_v (Sec 8.1 l.1999; the original radii, Q-35). Throws std::invalid_argument if the
/// list is empty or holds a value that is not positive and finite (a positive ordinary block).
template <class W>
W outer_block_scale(std::span<const W> radii);

/// l.2085: F^2 s_U >= a_U. Throws std::invalid_argument unless both are positive and finite.
template <class W>
bool outer_block_aligned(const OuterScaleFactor& f, W s, W a);

/// l.2118: a_U / F^2 <= s_U <= F a_U, i.e. aligned and not above F a_U (Lemma 8.2 gives the
/// upper part for every aligned block). Same preconditions as outer_block_aligned.
template <class W>
bool outer_aligned_scale_in_range(const OuterScaleFactor& f, W s, W a);

/// 4 K^2 h s with K = 2^log2_k (Lemma 8.4 item 1, l.2109: the bound on an edge inside a
/// two-shell window). Throws std::invalid_argument unless h >= 2 and s is positive and finite.
template <class W>
W outer_window_edge_bound(W s, std::uint64_t h, unsigned log2_k);

/// l.2085-2086: every MST edge of an unaligned block is below a_U / 2. Returns the index of the
/// first weight w with w + w >= a, or weights.size() if there is none. Throws
/// std::invalid_argument unless a is positive and finite.
template <class W>
std::size_t first_edge_not_below_half(std::span<const W> weights, W a);

/// l.2092: s_U / 2 < t_U <= s_U (half-open at s_U / 2, Q-36).
template <class W>
bool outer_cutoff_valid(W s, W t);

/// An aligned block as seen by the grid: its scale and the index of its name a_U.
template <class W>
struct OuterAlignedScale {
  W s{};
  std::uint32_t name = 0;
};

/// Support [lower, upper] of one name (l.2081, l.2086) with its grid indices (Lemma 8.4 item
/// 2): lower_index is the grid floor of lower; upper_index the first grid point strictly above
/// upper, or the saturated sentinel infinity() when upper itself saturated.
template <class W>
struct OuterSupport {
  W lower{};
  W upper{};
  std::size_t lower_index = 0;
  std::size_t upper_index = 0;
};

/// The global grid of Sec 8.2. `points` is strictly increasing. Component c owns the points
/// from component_begin[c] up to its sentinel; inside a component each point is the double of
/// the previous one (the last may be the saturated infinity()). When a component's sentinel
/// equals the next component's start, the point is stored once and shared. Class i is
/// [points[i], points[i + 1]) for every i that is not a sentinel.
template <class W>
struct OuterGlobalGrid {
  std::vector<W> points;
  std::vector<std::size_t> component_begin;
  std::vector<OuterSupport<W>> supports;  ///< per name, in input order
};

inline constexpr std::size_t kNoGridClass = std::numeric_limits<std::size_t>::max();

/// l.2081-2089 and Lemma 8.4 item 2. `names` holds the distinct finite join names and corridor
/// weights; `aligned` the aligned nonsingleton positive ordinary blocks.
/// Throws std::invalid_argument if a name or a scale is not positive and finite, if a block
/// refers to an unknown name, or if a block is not aligned (F^2 s < a). Throws std::logic_error
/// if an aligned block has s > F a (the relation of l.2118 fails).
template <class W>
OuterGlobalGrid<W> build_outer_global_grid(const OuterScaleFactor& f, std::span<const W> names,
                                           std::span<const OuterAlignedScale<W>> aligned);

/// Class of a value labelled with `name`: binary search over the support's index range
/// (O(log log F) comparisons). Returns i with points[i] <= w < points[i + 1], or kNoGridClass if
/// w is below the floor of the support start (a local weight, l.2122) or at or above the
/// support's upper sentinel. Throws std::invalid_argument for an unknown name.
template <class W>
std::size_t outer_grid_class(const OuterGlobalGrid<W>& grid, std::uint32_t name, W w);

/// l.2090: t_U of an aligned block, the floor of s_U on its name's global grid. Throws
/// std::invalid_argument for an unknown name and std::logic_error if s has no class there or
/// s / 2 < t <= s fails.
template <class W>
W outer_aligned_cutoff(const OuterGlobalGrid<W>& grid, std::uint32_t name, W s);

/// l.2091, l.2099-2103: t_U of a block on a separate grid, the floor of s_U on `table` (strictly
/// increasing doubling points, e.g. the extended piece table of Sec 7.3). Throws
/// std::logic_error if the table has no point <= s or if s / 2 < t <= s fails.
template <class W>
W outer_separate_cutoff(std::span<const W> table, W s);

}  // namespace sssp::paper::outer
