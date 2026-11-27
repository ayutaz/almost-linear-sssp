#pragma once
// Lengths in the comparison-addition model of Sec 2, "Operations on lengths and indices"
// (arXiv 2611.34911v1, txt l.292-300).
//
// - Len<W>: one opaque real register plus a separate infinity flag (l.293). Only copy, compare,
//   add and doubling are offered; there is deliberately no floor, log, bit access or use as an
//   index (l.300).
// - SignedPair<W>: a signed value p - q stored as two nonnegative sums (p, q) (l.294-296);
//   addition is componentwise, negation swaps, and p - q < r - s is decided by p + s < r + q.
// - EdgeKey<W>: the (length, edge index) order used to break ties between equal lengths.
// - LenOps<W>: the operations, counted into an optional sssp::OpCounters sink via SSSP_COUNT
//   (compiled out when SSSP_DIAGNOSTICS=0). Only real-register operations are counted: work on
//   the infinity flag or on edge indices is integer word work and is not a comparison/addition.
//
// Supported W: unsigned integers (exact; overflow throws std::overflow_error, and the value
// 2^w - 1 is reserved for WeightTraits<W>::infinity()) and floating point (IEEE-754 rounding;
// a sum that is not finite throws std::overflow_error).

#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::model {

template <class W>
concept LengthScalar = std::unsigned_integral<W> || std::floating_point<W>;

/// A nonnegative length or +infinity. Default-constructed value is the finite zero.
template <LengthScalar W>
class Len {
 public:
  using scalar_type = W;

  constexpr Len() noexcept = default;

  /// A finite length. Throws std::invalid_argument for NaN, negative or non-finite doubles and
  /// for the integer value reserved as the shared infinity sentinel. -0.0 becomes +0.0.
  static constexpr Len finite(W v) {
    if constexpr (std::floating_point<W>) {
      if (!std::isfinite(v) || v < W{0}) {
        throw std::invalid_argument("length must be finite and nonnegative");
      }
      if (v == W{0}) v = W{0};
    } else {
      if (v == WeightTraits<W>::infinity()) {
        throw std::invalid_argument("integer length collides with the infinity sentinel");
      }
    }
    return Len(v, false);
  }
  static constexpr Len infinity() noexcept { return Len(W{0}, true); }
  static constexpr Len zero() noexcept { return Len(W{0}, false); }

  /// Interoperation with CsrGraph weights and SsspResult::dist, where WeightTraits<W>::infinity()
  /// encodes "unreachable".
  static constexpr Len from_weight(W w) {
    return w == WeightTraits<W>::infinity() ? infinity() : finite(w);
  }
  [[nodiscard]] constexpr W to_weight() const noexcept {
    return inf_ ? WeightTraits<W>::infinity() : value_;
  }

  [[nodiscard]] constexpr bool is_infinite() const noexcept { return inf_; }
  [[nodiscard]] constexpr bool is_finite() const noexcept { return !inf_; }

  /// The register content; throws std::logic_error if the length is infinite.
  [[nodiscard]] constexpr W value() const {
    if (inf_) throw std::logic_error("value() of an infinite length");
    return value_;
  }

  /// Structural equality (uncounted; use LenOps::equal inside algorithms).
  friend constexpr bool operator==(const Len& a, const Len& b) noexcept {
    return a.inf_ == b.inf_ && (a.inf_ || a.value_ == b.value_);
  }
  /// Uncounted strict order for containers and tests; algorithms use LenOps::less.
  friend constexpr bool operator<(const Len& a, const Len& b) noexcept {
    if (a.inf_ || b.inf_) return !a.inf_ && b.inf_;
    return a.value_ < b.value_;
  }

 private:
  constexpr Len(W v, bool inf) noexcept : value_(v), inf_(inf) {}

  W value_{};
  bool inf_ = false;
};

/// The signed value positive() - negative() as two finite nonnegative sums (Sec 2, l.294-296).
/// Infinite components are rejected: inf - inf has no meaning in this representation.
template <LengthScalar W>
class SignedPair {
 public:
  constexpr SignedPair() noexcept = default;
  constexpr SignedPair(Len<W> positive, Len<W> negative) : pos_(positive), neg_(negative) {
    if (pos_.is_infinite() || neg_.is_infinite()) {
      throw std::invalid_argument("signed pair components must be finite");
    }
  }
  /// The value a, i.e. the pair (a, 0).
  static constexpr SignedPair from(Len<W> a) { return SignedPair(a, Len<W>::zero()); }

  [[nodiscard]] constexpr Len<W> positive() const noexcept { return pos_; }
  [[nodiscard]] constexpr Len<W> negative() const noexcept { return neg_; }
  /// Negation interchanges the components (Sec 2, l.295); no real operation is performed.
  [[nodiscard]] constexpr SignedPair neg() const noexcept { return SignedPair(neg_, pos_, 0); }

 private:
  constexpr SignedPair(Len<W> p, Len<W> q, int /*unchecked*/) noexcept : pos_(p), neg_(q) {}

  Len<W> pos_{};
  Len<W> neg_{};
};

/// (length, edge index): lexicographic strict total order that breaks ties between equal lengths
/// by the original edge index.
template <LengthScalar W>
struct EdgeKey {
  W w;
  EdgeId id;

  /// Uncounted comparisons for containers and tests; algorithms use LenOps::less.
  friend constexpr bool operator==(const EdgeKey& a, const EdgeKey& b) noexcept {
    return a.w == b.w && a.id == b.id;
  }
  friend constexpr bool operator<(const EdgeKey& a, const EdgeKey& b) noexcept {
    if (a.w < b.w) return true;
    if (b.w < a.w) return false;
    return a.id < b.id;
  }
};

/// Counted operations on Len, SignedPair and EdgeKey.
template <LengthScalar W>
class LenOps {
 public:
  constexpr LenOps() noexcept = default;
  explicit constexpr LenOps(OpCounters* sink) noexcept : sink_(sink) {}

  [[nodiscard]] constexpr OpCounters* counters() const noexcept { return sink_; }

  /// a + b; infinity absorbs (no real addition is performed then). One counted addition.
  [[nodiscard]] Len<W> add(Len<W> a, Len<W> b) const {
    if (a.is_infinite() || b.is_infinite()) return Len<W>::infinity();
    count_add(1);
    return Len<W>::finite(checked_sum(a.value(), b.value()));
  }

  /// a < b with every finite length below infinity. One counted comparison if both are finite.
  [[nodiscard]] bool less(Len<W> a, Len<W> b) const {
    if (a.is_infinite() || b.is_infinite()) return a.is_finite() && b.is_infinite();
    count_cmp(1);
    return a.value() < b.value();
  }
  [[nodiscard]] bool less_equal(Len<W> a, Len<W> b) const { return !less(b, a); }
  [[nodiscard]] bool equal(Len<W> a, Len<W> b) const {
    if (a.is_infinite() || b.is_infinite()) return a.is_infinite() && b.is_infinite();
    count_cmp(1);
    return a.value() == b.value();
  }

  /// a * 2^j by j doublings (Sec 2, l.297: "multiplication of a real by 2^j takes j
  /// doublings"); j counted additions. Exact for doubles unless it overflows (then throws).
  [[nodiscard]] Len<W> mul_pow2(Len<W> a, unsigned j) const {
    if (a.is_infinite()) return a;
    W v = a.value();
    for (unsigned i = 0; i < j; ++i) v = checked_sum(v, v);
    count_add(j);
    return Len<W>::finite(v);
  }

  /// Componentwise addition (Sec 2, l.295): two counted additions.
  [[nodiscard]] SignedPair<W> add(SignedPair<W> x, SignedPair<W> y) const {
    return SignedPair<W>(add(x.positive(), y.positive()), add(x.negative(), y.negative()));
  }
  /// x - y = x + (-y).
  [[nodiscard]] SignedPair<W> sub(SignedPair<W> x, SignedPair<W> y) const {
    return add(x, y.neg());
  }
  /// (a - b) < (c - d) decided as a + d < c + b (Sec 2, l.293-294): two additions, one comparison.
  [[nodiscard]] bool less(SignedPair<W> x, SignedPair<W> y) const {
    return less(add(x.positive(), y.negative()), add(y.positive(), x.negative()));
  }
  [[nodiscard]] bool equal(SignedPair<W> x, SignedPair<W> y) const {
    return equal(add(x.positive(), y.negative()), add(y.positive(), x.negative()));
  }

  /// (w, id) lexicographic order: one counted comparison when the lengths differ, two when they
  /// are equal (the index comparison is integer work and is not counted).
  [[nodiscard]] bool less(const EdgeKey<W>& a, const EdgeKey<W>& b) const {
    count_cmp(1);
    if (a.w < b.w) return true;
    count_cmp(1);
    if (b.w < a.w) return false;
    return a.id < b.id;
  }

 private:
  static W checked_sum(W a, W b) {
    if constexpr (std::floating_point<W>) {
      const W s = a + b;
      if (!std::isfinite(s)) throw std::overflow_error("floating-point length overflow");
      return s;
    } else {
      W s{};
      if (__builtin_add_overflow(a, b, &s) || s == WeightTraits<W>::infinity()) {
        throw std::overflow_error("integer length overflow");
      }
      return s;
    }
  }

  void count_add([[maybe_unused]] std::uint64_t k) const {
    SSSP_COUNT(if (sink_ != nullptr) sink_->additions += k);
  }
  void count_cmp([[maybe_unused]] std::uint64_t k) const {
    SSSP_COUNT(if (sink_ != nullptr) sink_->comparisons += k);
  }

  OpCounters* sink_ = nullptr;
};

extern template class Len<std::uint64_t>;
extern template class Len<double>;
extern template class SignedPair<std::uint64_t>;
extern template class SignedPair<double>;
extern template class LenOps<std::uint64_t>;
extern template class LenOps<double>;

}  // namespace sssp::paper::model
