#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <vector>

#include "baseline/heaps/binary_heap.hpp"
#include "sssp/result.hpp"

using sssp::baseline::heaps::BinaryHeap;

TEST_CASE("binary heap pops entries in non-decreasing key order", "[heap][binary]") {
  BinaryHeap<std::uint64_t, std::uint32_t> h(16);
  const std::vector<std::uint64_t> keys{5, 3, 8, 1, 9, 2, 7};
  for (std::uint32_t i = 0; i < keys.size(); ++i) h.push(keys[i], i);
  REQUIRE(h.size() == 7);
  REQUIRE(h.top().key == 1);
  REQUIRE(h.top().value == 3);
  const std::vector<std::uint64_t> expected_keys{1, 2, 3, 5, 7, 8, 9};
  const std::vector<std::uint32_t> expected_values{3, 5, 1, 0, 6, 2, 4};
  for (std::size_t i = 0; i < expected_keys.size(); ++i) {
    const auto e = h.pop();
    REQUIRE(e.key == expected_keys[i]);
    REQUIRE(e.value == expected_values[i]);
  }
  REQUIRE(h.empty());
}

TEST_CASE("binary heap keeps duplicate keys as separate entries", "[heap][binary]") {
  BinaryHeap<std::uint64_t, std::uint32_t> h(4);
  h.push(4, 10);
  h.push(4, 11);
  h.push(1, 12);
  h.push(4, 10);  // exact duplicate (lazy deletion pushes the same vertex twice)
  REQUIRE(h.pop().key == 1);
  std::multiset<std::uint32_t> values;
  for (int i = 0; i < 3; ++i) {
    const auto e = h.pop();
    REQUIRE(e.key == 4);
    values.insert(e.value);
  }
  REQUIRE(values == std::multiset<std::uint32_t>{10, 10, 11});
  REQUIRE(h.empty());
}

TEST_CASE("binary heap works with floating-point keys including zero", "[heap][binary]") {
  BinaryHeap<double, std::uint32_t> h(8);
  h.push(0.5, 0);
  h.push(0.0, 1);
  h.push(1e300, 2);
  h.push(0.25, 3);
  REQUIRE(h.pop().key == 0.0);
  REQUIRE(h.pop().key == 0.25);
  REQUIRE(h.pop().key == 0.5);
  REQUIRE(h.pop().key == 1e300);
}

TEST_CASE("binary heap does not reallocate within reserved capacity", "[heap][binary][alloc]") {
  BinaryHeap<std::uint64_t, std::uint32_t> h(64);
  const auto cap = h.capacity();
  REQUIRE(cap >= 64);
  for (std::uint32_t i = 0; i < 64; ++i) h.push(64 - i, i);
  REQUIRE(h.capacity() == cap);
  h.clear();
  REQUIRE(h.empty());
  REQUIRE(h.capacity() == cap);
}

TEST_CASE("binary heap rejects top/pop on an empty heap", "[heap][binary]") {
  BinaryHeap<std::uint64_t, std::uint32_t> h(1);
  REQUIRE_THROWS_AS(h.pop(), std::out_of_range);
  REQUIRE_THROWS_AS(h.top(), std::out_of_range);
}

#if SSSP_DIAGNOSTICS
TEST_CASE("binary heap counts key comparisons", "[heap][binary][diagnostics]") {
  BinaryHeap<std::uint64_t, std::uint32_t> h(4);
  h.push(5, 0);
  REQUIRE(h.comparisons() == 0);  // nothing to compare against
  h.push(3, 1);                   // one sift-up comparison with the root
  REQUIRE(h.comparisons() == 1);
}
#endif
