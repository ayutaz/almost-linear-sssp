#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "baseline/heaps/fibonacci_heap.hpp"
#include "sssp/result.hpp"

using sssp::baseline::heaps::FibonacciHeap;
using Fib = FibonacciHeap<std::uint64_t>;

TEST_CASE("Fibonacci heap pops in key order after decrease-key", "[heap][fibonacci]") {
  Fib h(6);
  REQUIRE(h.capacity() == 6);
  const std::vector<std::uint64_t> keys{50, 40, 30, 20, 10, 60};
  for (std::uint32_t i = 0; i < keys.size(); ++i) h.push(i, keys[i]);
  h.decrease_key(5, 5);
  h.decrease_key(0, 15);
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

TEST_CASE("Fibonacci heap consolidates 2^k equal-degree trees into one tree", "[heap][fibonacci]") {
  Fib h(9);
  for (std::uint32_t i = 0; i < 9; ++i) h.push(i, 10 * (i + 1));
  REQUIRE(h.root_count() == 9);  // insert is lazy
  REQUIRE(h.pop().id == 0);
  // 8 singleton trees are linked pairwise like a binary counter: one tree of degree 3.
  REQUIRE(h.root_count() == 1);
  REQUIRE(h.parent_of(1) == Fib::kNone);
  REQUIRE(h.degree(1) == 3);
  REQUIRE(h.validate());
}

TEST_CASE("Fibonacci heap performs a cascading cut when a marked node loses a second child",
          "[heap][fibonacci]") {
  Fib h(9);
  std::vector<std::uint64_t> key(9);
  for (std::uint32_t i = 0; i < 9; ++i) {
    key[i] = 10 * (i + 1);
    h.push(i, key[i]);
  }
  REQUIRE(h.pop().id == 0);  // now a single tree rooted at id 1 with children of degree 2, 1, 0
  std::uint32_t a = Fib::kNone;
  for (std::uint32_t v = 1; v < 9; ++v) {
    if (h.parent_of(v) == 1 && h.degree(v) == 2) a = v;
  }
  REQUIRE(a != Fib::kNone);
  std::vector<std::uint32_t> children;
  for (std::uint32_t v = 1; v < 9; ++v) {
    if (h.parent_of(v) == a) children.push_back(v);
  }
  REQUIRE(children.size() == 2);

  h.decrease_key(children[0], 1);
  key[children[0]] = 1;
  REQUIRE(h.parent_of(children[0]) == Fib::kNone);  // cut
  REQUIRE(h.is_marked(a));                          // a lost its first child
  REQUIRE(h.cascading_cuts() == 0);

  h.decrease_key(children[1], 2);
  key[children[1]] = 2;
  REQUIRE(h.parent_of(children[1]) == Fib::kNone);
  REQUIRE(h.parent_of(a) == Fib::kNone);  // a was marked: cascading cut to the root list
  REQUIRE_FALSE(h.is_marked(a));
  REQUIRE(h.cascading_cuts() == 1);
  REQUIRE(h.root_count() == 4);
  REQUIRE(h.validate());

  std::vector<std::uint32_t> expected;
  for (std::uint32_t v = 1; v < 9; ++v) expected.push_back(v);
  std::sort(expected.begin(), expected.end(),
            [&](std::uint32_t x, std::uint32_t y) { return key[x] < key[y]; });
  for (const auto v : expected) {
    const auto e = h.pop();
    REQUIRE(e.id == v);
    REQUIRE(e.key == key[v]);
    REQUIRE(h.validate());
  }
}

TEST_CASE("Fibonacci heap tracks membership, ties and re-insertion", "[heap][fibonacci]") {
  Fib h(3);
  h.push(1, 7);
  REQUIRE(h.contains(1));
  REQUIRE(h.pop().id == 1);
  REQUIRE_FALSE(h.contains(1));
  h.push(1, 3);
  h.push(2, 3);
  h.decrease_key(2, 3);
  const auto a = h.pop();
  const auto b = h.pop();
  REQUIRE(a.key == 3);
  REQUIRE(b.key == 3);
  REQUIRE(a.id != b.id);
  REQUIRE(h.empty());
}

TEST_CASE("Fibonacci heap rejects invalid operations", "[heap][fibonacci]") {
  Fib h(2);
  REQUIRE_THROWS_AS(h.pop(), std::out_of_range);
  REQUIRE_THROWS_AS(h.push(2, 1), std::out_of_range);
  h.push(0, 5);
  REQUIRE_THROWS_AS(h.push(0, 4), std::invalid_argument);
  REQUIRE_THROWS_AS(h.decrease_key(1, 4), std::invalid_argument);
  REQUIRE_THROWS_AS(h.decrease_key(0, 6), std::invalid_argument);
}

TEST_CASE("Fibonacci heap works with floating-point keys", "[heap][fibonacci]") {
  FibonacciHeap<double> h(5);
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
TEST_CASE("Fibonacci heap counts key comparisons", "[heap][fibonacci][diagnostics]") {
  Fib h(3);
  h.push(0, 5);
  REQUIRE(h.comparisons() == 0);  // first root becomes the minimum without comparing
  h.push(1, 3);                   // one comparison against the current minimum
  REQUIRE(h.comparisons() == 1);
}
#endif
