#pragma once
// Wide path lengths shared by the Mode A modules that sum 64-bit lengths along a path of a
// forest (child_order.hpp: link forest F, Sec 4.2; owner_decomposition.hpp: MST T, Sec 5.1).
//
// Integer W: unsigned __int128, so a sum of at most 2^34 lengths below 2^64 is exact. This is
// a numeric representation outside the paper's opaque-real model (deviation candidate, see the
// "Range of the arithmetic" notes of those headers). Floating point W: W itself.
// unsigned __int128 is not ISO C++; the single __extension__ below keeps -Wpedantic quiet.

#include <type_traits>

namespace sssp::paper {

namespace detail {
__extension__ typedef unsigned __int128 UInt128;  // not ISO C++; exact sums of 64-bit lengths
}  // namespace detail

/// unsigned __int128 for integer W (exact for every forest or tree), W for floating point.
template <class W>
using WideLength = std::conditional_t<std::is_integral_v<W>, detail::UInt128, W>;

}  // namespace sssp::paper
