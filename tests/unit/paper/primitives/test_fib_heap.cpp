// Tests for the boundary heap of Algorithm 1 (Sec 3.2, Lemma 3.3).
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/primitives/fib_heap.hpp"

using sssp::paper::FibHeap;
using Heap = FibHeap<std::uint64_t>;

namespace {

/// Depth of `id` (0 for a root) through the structural accessors.
std::size_t depth_of(const Heap& h, Heap::id_type id) {
  std::size_t d = 0;
  for (auto p = h.parent_of(id); p != Heap::kNone; p = h.parent_of(p)) ++d;
  return d;
}

}  // namespace

TEST_CASE("paper FibHeap extracts in key order after decreases", "[paper][fib_heap]") {
  Heap h(6);
  REQUIRE(h.capacity() == 6);
  const std::vector<std::uint64_t> keys{50, 40, 30, 20, 10, 60};
  for (std::uint32_t i = 0; i < keys.size(); ++i) h.insert(i, keys[i]);
  REQUIRE(h.size() == 6);
  h.decrease(5, 5);
  h.decrease(0, 15);
  h.decrease(1, 40);  // equal key: allowed, no change
  REQUIRE(h.key(0) == 15);
  REQUIRE(h.validate());
  const std::vector<std::uint64_t> expected_keys{5, 10, 15, 20, 30, 40};
  const std::vector<std::uint32_t> expected_ids{5, 4, 0, 3, 2, 1};
  for (std::size_t i = 0; i < expected_ids.size(); ++i) {
    const auto e = h.extract_min();
    REQUIRE(e.key == expected_keys[i]);
    REQUIRE(e.id == expected_ids[i]);
    REQUIRE_FALSE(h.contains(e.id));
    REQUIRE(h.size() == expected_ids.size() - i - 1);
    REQUIRE(h.validate());
  }
  REQUIRE(h.empty());
}

TEST_CASE("paper FibHeap links equal-degree roots on extraction", "[paper][fib_heap]") {
  Heap h(9);
  for (std::uint32_t i = 0; i < 9; ++i) h.insert(i, 10 * (i + 1));
  REQUIRE(h.root_count() == 9);  // insertion is lazy (new singleton roots)
  REQUIRE(h.extract_min().id == 0);
  // Eight singleton roots are linked pairwise: one tree of degree 3 rooted at the minimum.
  REQUIRE(h.root_count() == 1);
  REQUIRE(h.parent_of(1) == Heap::kNone);
  REQUIRE(h.degree(1) == 3);
  REQUIRE(h.validate());
}

TEST_CASE("paper FibHeap marks on the first child loss and cuts on the second",
          "[paper][fib_heap]") {
  // 17 singletons, extract one: the remaining 16 become a single binomial tree of degree 4.
  Heap h(17);
  for (std::uint32_t i = 0; i < 17; ++i) h.insert(i, 100 + i);
  REQUIRE(h.extract_min().id == 0);
  REQUIRE(h.root_count() == 1);
  REQUIRE(h.degree(1) == 4);
  // Find a node p at depth 2 with two children: losing both children cuts p (cascade).
  Heap::id_type p = Heap::kNone;
  for (std::uint32_t i = 1; i < 17; ++i) {
    if (depth_of(h, i) == 2 && h.degree(i) == 2) p = i;
  }
  REQUIRE(p != Heap::kNone);
  const Heap::id_type grand = h.parent_of(p);
  std::vector<Heap::id_type> kids;
  for (std::uint32_t i = 1; i < 17; ++i) {
    if (h.parent_of(i) == p) kids.push_back(i);
  }
  REQUIRE(kids.size() == 2);
  REQUIRE_FALSE(h.is_marked(p));
  h.decrease(kids[0], 1);
  REQUIRE(h.parent_of(kids[0]) == Heap::kNone);
  REQUIRE(h.is_marked(p));
  REQUIRE(h.cascading_cuts() == 0);
  REQUIRE(h.validate());
  h.decrease(kids[1], 2);
  REQUIRE(h.cascading_cuts() == 1);
  REQUIRE(h.parent_of(p) == Heap::kNone);  // p was cut to the root list
  REQUIRE_FALSE(h.is_marked(p));           // roots are unmarked
  REQUIRE(h.is_marked(grand));             // its parent lost its first child
  REQUIRE(h.validate());
  REQUIRE(h.extract_min().id == kids[0]);
  REQUIRE(h.extract_min().id == kids[1]);
  REQUIRE(h.validate());
}

TEST_CASE("paper FibHeap size counts distinct ids", "[paper][fib_heap]") {
  Heap h(4);
  h.insert(2, 7);
  h.insert(3, 9);
  h.decrease(3, 8);
  h.decrease(3, 1);
  REQUIRE(h.size() == 2);  // decreases never add entries
  REQUIRE(h.contains(2));
  REQUIRE(h.contains(3));
  REQUIRE_FALSE(h.contains(0));
  REQUIRE_FALSE(h.contains(99));
}

TEST_CASE("paper FibHeap discard releases every entry and starts a fresh generation",
          "[paper][fib_heap]") {
  Heap h(5);
  for (std::uint32_t i = 0; i < 5; ++i) h.insert(i, 50 - i);
  (void)h.extract_min();  // builds trees, so discard must drop non-root nodes too
  h.decrease(0, 3);
  h.discard();
  REQUIRE(h.empty());
  REQUIRE(h.size() == 0);
  REQUIRE(h.root_count() == 0);
  for (std::uint32_t i = 0; i < 5; ++i) REQUIRE_FALSE(h.contains(i));
  REQUIRE(h.validate());
  // Ids are reusable with unrelated keys; nothing from the previous generation leaks.
  h.insert(3, 1000);
  h.insert(0, 2000);
  REQUIRE(h.size() == 2);
  REQUIRE(h.parent_of(3) == Heap::kNone);
  REQUIRE(h.degree(3) == 0);
  REQUIRE(h.extract_min().id == 3);
  REQUIRE(h.extract_min().id == 0);
  REQUIRE(h.empty());
  h.discard();  // discarding an empty heap is fine
  REQUIRE(h.empty());
}

TEST_CASE("paper FibHeap degree array has the size 3 + 2 ceil(log2(s + 1)) of Lemma 3.3",
          "[paper][fib_heap]") {
  REQUIRE(FibHeap<std::uint64_t>(0).degree_array_size() == 3);   // ceil(log2 1) = 0
  REQUIRE(FibHeap<std::uint64_t>(1).degree_array_size() == 5);   // ceil(log2 2) = 1
  REQUIRE(FibHeap<std::uint64_t>(3).degree_array_size() == 7);   // ceil(log2 4) = 2
  REQUIRE(FibHeap<std::uint64_t>(4).degree_array_size() == 9);   // ceil(log2 5) = 3
  REQUIRE(FibHeap<std::uint64_t>(7).degree_array_size() == 9);   // ceil(log2 8) = 3
  REQUIRE(FibHeap<std::uint64_t>(8).degree_array_size() == 11);  // ceil(log2 9) = 4
  REQUIRE(FibHeap<std::uint64_t>(1000).degree_array_size() == 23);
}

TEST_CASE("paper FibHeap rejects invalid operations", "[paper][fib_heap]") {
  Heap h(3);
  REQUIRE_THROWS_AS(h.extract_min(), std::out_of_range);
  REQUIRE_THROWS_AS(h.insert(3, 1), std::out_of_range);
  REQUIRE_THROWS_AS(h.decrease(1, 1), std::invalid_argument);  // absent
  h.insert(1, 10);
  REQUIRE_THROWS_AS(h.insert(1, 5), std::invalid_argument);     // already present
  REQUIRE_THROWS_AS(h.decrease(1, 11), std::invalid_argument);  // key increase
  REQUIRE_THROWS_AS(h.key(0), std::invalid_argument);
  REQUIRE(h.key(1) == 10);
}

TEST_CASE("paper FibHeap orders composite keys lexicographically", "[paper][fib_heap]") {
  using Key = std::pair<double, std::uint32_t>;
  FibHeap<Key> h(4);
  h.insert(0, Key{1.5, 9});
  h.insert(1, Key{1.5, 2});
  h.insert(2, Key{0.25, 7});
  h.insert(3, Key{2.0, 0});
  h.decrease(3, Key{1.5, 1});
  REQUIRE(h.extract_min().id == 2);
  REQUIRE(h.extract_min().id == 3);
  REQUIRE(h.extract_min().id == 1);
  REQUIRE(h.extract_min().id == 0);
}

TEST_CASE("paper FibHeap counts key comparisons", "[paper][fib_heap]") {
  Heap h(3);
  REQUIRE(h.comparisons() == 0);
  h.insert(0, 3);
  h.insert(1, 2);
  h.insert(2, 1);
  (void)h.extract_min();
  REQUIRE(h.comparisons() > 0);
}

TEST_CASE("paper FibHeap matches an ordered-set model under random operations with discards",
          "[paper][fib_heap]") {
  for (std::uint64_t seed : {1ULL, 2ULL, 3ULL, 4ULL, 5ULL}) {
    std::mt19937_64 rng(seed);
    const std::uint32_t n = 1 + static_cast<std::uint32_t>(rng() % 60);
    Heap h(n);
    std::set<std::pair<std::uint64_t, std::uint32_t>> model;
    std::vector<std::uint64_t> key(n, 0);
    std::vector<bool> present(n, false);
    for (int step = 0; step < 3000; ++step) {
      const auto op = rng() % 100;
      const auto id = static_cast<std::uint32_t>(rng() % n);
      if (op < 40) {
        if (!present[id]) {
          key[id] = rng() % 50;
          h.insert(id, key[id]);
          model.emplace(key[id], id);
          present[id] = true;
        }
      } else if (op < 70) {
        if (present[id]) {
          const std::uint64_t nk = key[id] == 0 ? 0 : rng() % (key[id] + 1);
          model.erase({key[id], id});
          key[id] = nk;
          model.emplace(nk, id);
          h.decrease(id, nk);
        }
      } else if (op < 98) {
        if (!model.empty()) {
          const auto e = h.extract_min();
          // Ties may be broken in any order: only the key must be minimum.
          REQUIRE(e.key == model.begin()->first);
          REQUIRE(present[e.id]);
          REQUIRE(key[e.id] == e.key);
          model.erase({e.key, e.id});
          present[e.id] = false;
        }
      } else {
        h.discard();
        model.clear();
        present.assign(n, false);
      }
      REQUIRE(h.size() == model.size());
      if (step % 37 == 0) REQUIRE(h.validate());
    }
    REQUIRE(h.validate());
  }
}
