#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "baseline/heaps/dary_heap.hpp"
#include "sssp/result.hpp"

using sssp::baseline::heaps::IndexedDaryHeap;

using Dary2 = IndexedDaryHeap<std::uint64_t, 2>;
using Dary3 = IndexedDaryHeap<std::uint64_t, 3>;
using Dary4 = IndexedDaryHeap<std::uint64_t, 4>;
using Dary8 = IndexedDaryHeap<std::uint64_t, 8>;

TEMPLATE_TEST_CASE("indexed d-ary heap pops in key order after decrease-key", "[heap][dary]", Dary2,
                   Dary3, Dary4, Dary8) {
  TestType h(6);
  REQUIRE(h.capacity() == 6);
  const std::vector<std::uint64_t> keys{50, 40, 30, 20, 10, 60};
  for (std::uint32_t i = 0; i < keys.size(); ++i) h.push(i, keys[i]);
  h.decrease_key(5, 5);
  h.decrease_key(0, 15);
  REQUIRE(h.key(0) == 15);
  REQUIRE(h.validate());
  const std::vector<std::uint64_t> expected_keys{5, 10, 15, 20, 30, 40};
  const std::vector<std::uint32_t> expected_ids{5, 4, 0, 3, 2, 1};
  for (std::size_t i = 0; i < expected_ids.size(); ++i) {
    const auto e = h.pop();
    REQUIRE(e.key == expected_keys[i]);
    REQUIRE(e.id == expected_ids[i]);
    REQUIRE_FALSE(h.contains(e.id));
    REQUIRE(h.validate());
  }
  REQUIRE(h.empty());
}

TEMPLATE_TEST_CASE("indexed d-ary heap tracks membership and allows re-insertion", "[heap][dary]",
                   Dary2, Dary4) {
  TestType h(3);
  REQUIRE_FALSE(h.contains(1));
  h.push(1, 7);
  REQUIRE(h.contains(1));
  REQUIRE(h.size() == 1);
  REQUIRE(h.pop().id == 1);
  REQUIRE_FALSE(h.contains(1));
  h.push(1, 3);          // re-insertion of a popped id
  h.push(2, 3);          // tie
  h.decrease_key(2, 3);  // decrease to the same key is allowed
  const auto a = h.pop();
  const auto b = h.pop();
  REQUIRE(a.key == 3);
  REQUIRE(b.key == 3);
  REQUIRE(a.id != b.id);
}

TEMPLATE_TEST_CASE("indexed d-ary heap rejects invalid operations", "[heap][dary]", Dary2, Dary4) {
  TestType h(2);
  REQUIRE_THROWS_AS(h.pop(), std::out_of_range);
  REQUIRE_THROWS_AS(h.push(2, 1), std::out_of_range);
  h.push(0, 5);
  REQUIRE_THROWS_AS(h.push(0, 4), std::invalid_argument);          // already present
  REQUIRE_THROWS_AS(h.decrease_key(1, 4), std::invalid_argument);  // not present
  REQUIRE_THROWS_AS(h.decrease_key(0, 6), std::invalid_argument);  // key increase
}

TEST_CASE("indexed 4-ary heap works with floating-point keys", "[heap][dary]") {
  IndexedDaryHeap<double, 4> h(5);
  h.push(0, 0.75);
  h.push(1, 0.0);
  h.push(2, 1e300);
  h.push(3, 0.5);
  h.decrease_key(2, 0.25);
  REQUIRE(h.pop().id == 1);
  REQUIRE(h.pop().id == 2);
  REQUIRE(h.pop().id == 3);
  REQUIRE(h.pop().id == 0);
}

#if SSSP_DIAGNOSTICS
TEST_CASE("indexed 4-ary heap counts key comparisons", "[heap][dary][diagnostics]") {
  Dary4 h(4);
  h.push(0, 5);
  REQUIRE(h.comparisons() == 0);
  h.push(1, 3);  // one sift-up comparison against the root
  REQUIRE(h.comparisons() == 1);
}
#endif
