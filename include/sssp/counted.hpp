#pragma once
// Counted<W>: a weight wrapper that counts operations of the comparison-addition model.
//
// The paper (Sec 2) measures cost in a model where real numbers may only be compared and
// added. Instantiating an algorithm with Counted<double> or Counted<std::uint64_t> therefore
// (1) counts its comparisons and additions on weights/distances exactly, and
// (2) proves at compile time that it uses no other arithmetic on them: Counted deliberately
//     has no subtraction, multiplication, division, rounding or conversion back to W
//     (value() is the explicit escape hatch for I/O and tests).
//
// Counters are thread_local: each thread counts its own operations. A parallel algorithm must
// sum per-thread deltas itself (e.g. CountedScope in every worker).
// Constructions, copies and value() are free; only + / += and the six comparisons count.
// Comparisons against the infinity sentinel count too, as they would in the model.

#include <cstdint>
#include <stdexcept>
#include <type_traits>

#include "sssp/graph.hpp"

namespace sssp {

struct CountedStats {
  std::uint64_t additions = 0;
  std::uint64_t comparisons = 0;

  friend constexpr CountedStats operator-(const CountedStats& a, const CountedStats& b) noexcept {
    return {a.additions - b.additions, a.comparisons - b.comparisons};
  }
  friend constexpr bool operator==(const CountedStats&, const CountedStats&) noexcept = default;
};

/// The calling thread's counters.
inline CountedStats& counted_stats() noexcept {
  thread_local CountedStats stats;
  return stats;
}

inline void reset_counted_stats() noexcept {
  counted_stats() = {};
}

/// Records the calling thread's counters at construction; delta() is the work done since.
class CountedScope {
 public:
  CountedScope() noexcept : start_(counted_stats()) {}
  [[nodiscard]] CountedStats delta() const noexcept { return counted_stats() - start_; }

 private:
  CountedStats start_;
};

template <class W>
class Counted {
  static_assert(std::is_arithmetic_v<W>, "Counted wraps an arithmetic weight type");

 public:
  using value_type = W;

  constexpr Counted() noexcept = default;
  constexpr Counted(W v) noexcept : v_(v) {}  // NOLINT(google-explicit-constructor): W -> Counted

  [[nodiscard]] constexpr W value() const noexcept { return v_; }

  friend Counted operator+(Counted a, Counted b) noexcept {
    ++counted_stats().additions;
    return Counted(static_cast<W>(a.v_ + b.v_));
  }
  Counted& operator+=(Counted b) noexcept {
    ++counted_stats().additions;
    v_ = static_cast<W>(v_ + b.v_);
    return *this;
  }

  friend bool operator<(Counted a, Counted b) noexcept { return count(a.v_ < b.v_); }
  friend bool operator<=(Counted a, Counted b) noexcept { return count(a.v_ <= b.v_); }
  friend bool operator>(Counted a, Counted b) noexcept { return count(a.v_ > b.v_); }
  friend bool operator>=(Counted a, Counted b) noexcept { return count(a.v_ >= b.v_); }
  friend bool operator==(Counted a, Counted b) noexcept { return count(a.v_ == b.v_); }
  friend bool operator!=(Counted a, Counted b) noexcept { return count(a.v_ != b.v_); }

 private:
  static bool count(bool r) noexcept {
    ++counted_stats().comparisons;
    return r;
  }

  W v_{};
};

/// Weight semantics of Counted<W> are those of W. checked_add counts one model addition; its
/// overflow guard is bookkeeping outside the model and is not counted.
template <class W>
struct WeightTraits<Counted<W>, void> {
  using Base = WeightTraits<W>;
  static constexpr bool is_exact = Base::is_exact;
  static constexpr Counted<W> infinity() noexcept { return Counted<W>(Base::infinity()); }
  static constexpr Counted<W> zero() noexcept { return Counted<W>(Base::zero()); }
  static bool is_valid_weight(Counted<W> w) noexcept { return Base::is_valid_weight(w.value()); }
  static Counted<W> checked_add(Counted<W> a, Counted<W> b) {
    const W sum = Base::checked_add(a.value(), b.value());
    ++counted_stats().additions;
    return Counted<W>(sum);
  }
};

}  // namespace sssp
