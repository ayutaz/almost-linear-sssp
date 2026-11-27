// sssp_bench_runner: one benchmark job (graph x algorithms x sources x trials) -> JSONL records.
//
//   sssp_bench_runner --type er --n 100000 --avg-degree 4 --weights uniform_int --seed 1
//     --algorithms dijkstra-binary --trials 5 --warmup 1 --env-json env.json --run-id ID
//     --out results.jsonl --isolate
//
// CLI: docs/design/benchmark_runner_cli.md; records: docs/design/benchmark_schema.md.
// The logic lives in the library (include/sssp/bench_runner.hpp) so that it is unit-tested.

#include <iostream>
#include <new>
#include <string>
#include <vector>

#include "sssp/bench_runner.hpp"

int main(int argc, char** argv) {
  try {
    const std::vector<std::string> args(argv + 1, argv + argc);
    return sssp::bench::run_main(args, std::cout, std::cerr);  // reports its own errors
  } catch (const std::bad_alloc&) {
    std::cerr << "sssp_bench_runner: out of memory (std::bad_alloc)\n";
    return sssp::bench::kExitOutOfMemory;
  } catch (...) {
    std::cerr << "sssp_bench_runner: unexpected error\n";
    return sssp::bench::kExitError;
  }
}
