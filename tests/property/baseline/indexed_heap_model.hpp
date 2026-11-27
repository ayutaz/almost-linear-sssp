#pragma once
// Reference model for indexed (decrease-key) heaps: a std::map from id to current key.
// A random sequence of push / decrease-key / pop operations is applied to both and every
// observable result is compared. If the heap has validate(), its invariants are checked too.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <map>
#include <random>

namespace sssp::test {

struct IndexedHeapModelConfig {
  std::uint64_t seed = 1;
  std::uint32_t ids = 64;
  std::size_t ops = 3000;
  /// Radix heaps require every new key to be >= the last popped key.
  bool monotone = false;
};

/// gen_push(rng, lo) -> key >= lo; gen_decrease(rng, lo, cur) -> key in [lo, cur].
template <class Heap, class GenPush, class GenDecrease>
void run_indexed_heap_model(const IndexedHeapModelConfig& cfg, GenPush gen_push,
                            GenDecrease gen_decrease) {
  using K = typename Heap::key_type;
  std::mt19937_64 rng(cfg.seed);
  Heap heap(cfg.ids);
  std::map<std::uint32_t, K> ref;
  K last_popped{};  // lower bound for monotone heaps
  std::uniform_int_distribution<std::uint32_t> pick_id(0, cfg.ids - 1);
  std::uniform_int_distribution<int> pick_op(0, 99);

  const auto check_pop = [&] {
    const auto e = heap.pop();
    const auto min_it = std::min_element(
        ref.begin(), ref.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
    REQUIRE(e.key == min_it->second);
    const auto it = ref.find(e.id);
    REQUIRE(it != ref.end());
    REQUIRE(it->second == e.key);
    ref.erase(it);
    REQUIRE_FALSE(heap.contains(e.id));
    last_popped = e.key;
  };

  for (std::size_t step = 0; step < cfg.ops; ++step) {
    const int op = pick_op(rng);
    const K lo = cfg.monotone ? last_popped : K{};
    if (op < 40) {  // push a random id that is not in the heap
      const std::uint32_t id = pick_id(rng);
      if (ref.count(id) == 0) {
        const K key = gen_push(rng, lo);
        heap.push(id, key);
        ref[id] = key;
      }
    } else if (op < 75) {  // decrease-key of a random id in the heap
      if (!ref.empty()) {
        auto it = ref.begin();
        std::advance(it, std::uniform_int_distribution<std::size_t>(0, ref.size() - 1)(rng));
        const K key = gen_decrease(rng, lo, it->second);
        heap.decrease_key(it->first, key);
        it->second = key;
      }
    } else if (!ref.empty()) {
      check_pop();
    }
    REQUIRE(heap.size() == ref.size());
    const std::uint32_t probe = pick_id(rng);
    REQUIRE(heap.contains(probe) == (ref.count(probe) == 1));
    if (ref.count(probe) == 1) REQUIRE(heap.key(probe) == ref[probe]);
    if constexpr (requires { heap.validate(); }) {
      if (step % 17 == 0) REQUIRE(heap.validate());
    }
  }
  while (!ref.empty()) check_pop();
  REQUIRE(heap.empty());
  if constexpr (requires { heap.validate(); }) REQUIRE(heap.validate());
}

inline constexpr std::uint64_t kHeapSeeds[] = {1, 7, 42, 2026, 0xC0FFEEULL};

/// Key generators: integers in a window above lo (small window -> many ties), doubles.
inline auto int_push_gen(std::uint64_t window) {
  return [window](std::mt19937_64& r, std::uint64_t lo) {
    constexpr auto kMax = std::numeric_limits<std::uint64_t>::max();
    const std::uint64_t hi = lo > kMax - window ? kMax : lo + window;
    return std::uniform_int_distribution<std::uint64_t>(lo, hi)(r);
  };
}
inline auto int_decrease_gen() {
  return [](std::mt19937_64& r, std::uint64_t lo, std::uint64_t cur) {
    return std::uniform_int_distribution<std::uint64_t>(lo, cur)(r);
  };
}
inline auto real_push_gen(double window) {
  return [window](std::mt19937_64& r, double lo) {
    return lo + std::uniform_real_distribution<double>(0.0, window)(r);
  };
}
inline auto real_decrease_gen() {
  return [](std::mt19937_64& r, double lo, double cur) {
    if (!(lo < cur)) return cur;
    return std::uniform_real_distribution<double>(lo, cur)(r);
  };
}

}  // namespace sssp::test
