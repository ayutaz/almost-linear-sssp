// Property: indexed heaps (d-ary, Fibonacci, radix) agree with a std::map reference model on
// random push / decrease-key / pop sequences, for several fixed seeds and key distributions.
#include <catch2/catch_test_macros.hpp>
#include <cstdint>

#include "baseline/heaps/dary_heap.hpp"
#include "baseline/heaps/fibonacci_heap.hpp"
#include "baseline/heaps/radix_heap.hpp"
#include "property/baseline/indexed_heap_model.hpp"

using namespace sssp::baseline::heaps;
using sssp::test::IndexedHeapModelConfig;
using sssp::test::int_decrease_gen;
using sssp::test::int_push_gen;
using sssp::test::kHeapSeeds;
using sssp::test::real_decrease_gen;
using sssp::test::real_push_gen;
using sssp::test::run_indexed_heap_model;

namespace {

template <class Heap>
void run_integer_models(bool monotone) {
  for (const auto seed : kHeapSeeds) {
    for (const std::uint32_t ids : {8U, 64U, 300U}) {
      for (const std::uint64_t window :
           {std::uint64_t{4}, std::uint64_t{1} << 20, std::uint64_t{1} << 62}) {
        CAPTURE(seed, ids, window);
        run_indexed_heap_model<Heap>({seed, ids, 2000, monotone}, int_push_gen(window),
                                     int_decrease_gen());
      }
    }
  }
}

template <class Heap>
void run_real_models() {
  for (const auto seed : kHeapSeeds) {
    for (const std::uint32_t ids : {8U, 64U, 300U}) {
      CAPTURE(seed, ids);
      run_indexed_heap_model<Heap>({seed, ids, 2000, false}, real_push_gen(1000.0),
                                   real_decrease_gen());
    }
  }
}

}  // namespace

TEST_CASE("indexed 4-ary heap matches the reference model", "[property][heap][dary]") {
  run_integer_models<IndexedDaryHeap<std::uint64_t, 4>>(false);
  run_real_models<IndexedDaryHeap<double, 4>>();
}

TEST_CASE("indexed binary and 8-ary heaps match the reference model", "[property][heap][dary]") {
  run_integer_models<IndexedDaryHeap<std::uint64_t, 2>>(false);
  run_integer_models<IndexedDaryHeap<std::uint64_t, 8>>(false);
}

TEST_CASE("Fibonacci heap matches the reference model", "[property][heap][fibonacci]") {
  run_integer_models<FibonacciHeap<std::uint64_t>>(false);
  run_real_models<FibonacciHeap<double>>();
}

TEST_CASE("radix heap matches the reference model under monotone keys", "[property][heap][radix]") {
  run_integer_models<RadixHeap>(true);
}
