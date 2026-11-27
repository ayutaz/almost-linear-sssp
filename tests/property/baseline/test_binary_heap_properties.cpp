// Property: on random push/pop sequences the lazy binary heap behaves like a sorted multiset.
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <set>
#include <utility>

#include "baseline/heaps/binary_heap.hpp"

using sssp::baseline::heaps::BinaryHeap;

namespace {

template <class K, class KeyGen>
void run_binary_heap_model(std::uint64_t seed, std::size_t ops, KeyGen gen_key) {
  std::mt19937_64 rng(seed);
  BinaryHeap<K, std::uint32_t> heap(ops);
  std::multiset<std::pair<K, std::uint32_t>> ref;
  std::bernoulli_distribution do_push(0.6);
  std::uint32_t next_value = 0;
  for (std::size_t step = 0; step < ops; ++step) {
    if (ref.empty() || do_push(rng)) {
      const K key = gen_key(rng);
      heap.push(key, next_value);
      ref.emplace(key, next_value);
      ++next_value;
    } else {
      const auto e = heap.pop();
      REQUIRE(e.key == ref.begin()->first);
      const auto it = ref.find({e.key, e.value});
      REQUIRE(it != ref.end());
      ref.erase(it);
    }
    REQUIRE(heap.size() == ref.size());
  }
  while (!ref.empty()) {
    const auto e = heap.pop();
    REQUIRE(e.key == ref.begin()->first);
    const auto it = ref.find({e.key, e.value});
    REQUIRE(it != ref.end());
    ref.erase(it);
  }
  REQUIRE(heap.empty());
}

constexpr std::uint64_t kSeeds[] = {1, 7, 42, 2026, 0xC0FFEEULL};

}  // namespace

TEST_CASE("binary heap matches a sorted multiset (integer keys, many ties)",
          "[property][heap][binary]") {
  for (const auto seed : kSeeds) {
    CAPTURE(seed);
    std::uniform_int_distribution<std::uint64_t> key(0, 7);
    run_binary_heap_model<std::uint64_t>(seed, 3000, [&](auto& r) { return key(r); });
  }
}

TEST_CASE("binary heap matches a sorted multiset (full-range integer keys)",
          "[property][heap][binary]") {
  for (const auto seed : kSeeds) {
    CAPTURE(seed);
    std::uniform_int_distribution<std::uint64_t> key;
    run_binary_heap_model<std::uint64_t>(seed, 3000, [&](auto& r) { return key(r); });
  }
}

TEST_CASE("binary heap matches a sorted multiset (double keys)", "[property][heap][binary]") {
  for (const auto seed : kSeeds) {
    CAPTURE(seed);
    std::uniform_real_distribution<double> key(0.0, 1000.0);
    run_binary_heap_model<double>(seed, 3000, [&](auto& r) { return key(r); });
  }
}
