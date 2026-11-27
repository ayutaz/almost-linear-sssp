#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>

#include "baseline/heaps/radix_heap.hpp"
#include "sssp/result.hpp"

using sssp::baseline::heaps::RadixHeap;

TEST_CASE("radix heap pops monotone 64-bit keys in order", "[heap][radix]") {
  RadixHeap h(6);
  REQUIRE(h.capacity() == 6);
  constexpr std::uint64_t big = (std::uint64_t{1} << 63) + 1;
  h.push(0, 5);
  h.push(1, 3);
  h.push(2, 100);
  h.push(3, 3);
  h.push(4, std::uint64_t{1} << 40);
  h.push(5, big);
  REQUIRE(h.validate());
  std::set<std::uint32_t> tied;
  auto e = h.pop();
  REQUIRE(e.key == 3);
  tied.insert(e.id);
  e = h.pop();
  REQUIRE(e.key == 3);
  tied.insert(e.id);
  REQUIRE(tied == std::set<std::uint32_t>{1, 3});
  REQUIRE(h.pop().id == 0);
  REQUIRE(h.pop().id == 2);
  REQUIRE(h.pop().id == 4);
  e = h.pop();
  REQUIRE(e.id == 5);
  REQUIRE(e.key == big);
  REQUIRE(h.empty());
  REQUIRE(h.validate());
}

TEST_CASE("radix heap decrease-key moves an item to a lower bucket", "[heap][radix]") {
  RadixHeap h(3);
  h.push(0, 1000);
  h.push(1, 500);
  h.decrease_key(0, 7);
  REQUIRE(h.key(0) == 7);
  REQUIRE(h.validate());
  const auto e = h.pop();
  REQUIRE(e.id == 0);
  REQUIRE(e.key == 7);
  REQUIRE(h.pop().id == 1);
}

TEST_CASE("radix heap accepts keys equal to the last popped key (zero-weight edges)",
          "[heap][radix]") {
  RadixHeap h(3);
  h.push(0, 5);
  REQUIRE(h.pop().key == 5);
  h.push(1, 5);
  h.push(2, 6);
  h.decrease_key(2, 5);
  REQUIRE(h.pop().key == 5);
  REQUIRE(h.pop().key == 5);
}

TEST_CASE("radix heap enforces monotonicity and rejects invalid operations", "[heap][radix]") {
  RadixHeap h(3);
  REQUIRE_THROWS_AS(h.pop(), std::out_of_range);
  REQUIRE_THROWS_AS(h.push(3, 1), std::out_of_range);
  h.push(0, 5);
  h.push(1, 9);
  REQUIRE_THROWS_AS(h.push(0, 6), std::invalid_argument);           // already present
  REQUIRE_THROWS_AS(h.decrease_key(2, 4), std::invalid_argument);   // not present
  REQUIRE_THROWS_AS(h.decrease_key(1, 10), std::invalid_argument);  // key increase
  REQUIRE(h.pop().key == 5);
  REQUIRE_THROWS_AS(h.push(2, 4), std::invalid_argument);          // below last popped key
  REQUIRE_THROWS_AS(h.decrease_key(1, 4), std::invalid_argument);  // below last popped key
}

TEST_CASE("radix heap handles the full 64-bit key range", "[heap][radix]") {
  constexpr auto kMax = std::numeric_limits<std::uint64_t>::max();
  RadixHeap h(3);
  h.push(0, kMax);
  h.push(1, kMax - 1);
  h.push(2, 0);
  REQUIRE(h.pop().key == 0);
  REQUIRE(h.pop().key == kMax - 1);
  REQUIRE(h.pop().key == kMax);
}

#if SSSP_DIAGNOSTICS
TEST_CASE("radix heap counts comparisons only when searching a bucket minimum",
          "[heap][radix][diagnostics]") {
  RadixHeap h(3);
  h.push(0, 4);  // 4, 5, 6 all have bit width 3 relative to last = 0 -> same bucket
  h.push(1, 5);
  h.push(2, 6);
  REQUIRE(h.comparisons() == 0);
  REQUIRE(h.pop().key == 4);  // minimum of a 3-element bucket: 2 comparisons
  REQUIRE(h.comparisons() == 2);
  REQUIRE(h.pop().key == 5);  // singleton buckets need no comparison
  REQUIRE(h.pop().key == 6);
  REQUIRE(h.comparisons() == 2);
}
#endif
