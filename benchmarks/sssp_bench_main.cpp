// sssp_bench: Google Benchmark driver for every algorithm of the bench_runner registry.
//
//   sssp_bench [Google Benchmark flags] [--sssp_<option>=<value> ...]
//
// Google Benchmark flags (--benchmark_repetitions, --benchmark_min_warmup_time,
// --benchmark_enable_random_interleaving, --benchmark_out, ...) apply to every algorithm alike.
// Driver options: sssp_bench --sssp_help. Fairness rules: sssp_bench_driver.hpp.
//
// Exit codes: 0 every case checked so far is correct; 1 a case failed its check (wrong result or
// solver error; never timed); 2 bad arguments or an input that cannot be generated.

#include <benchmark/benchmark.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "sssp/bench_runner.hpp"
#include "sssp/version.hpp"
#include "sssp_bench_driver.hpp"

int main(int argc, char** argv) {
  benchmark::Initialize(&argc, argv);  // removes the --benchmark_* flags it recognises
  const std::vector<std::string> rest(argv + 1, argv + argc);

  namespace gb = sssp::bench::gbench;
  gb::DriverConfig config;
  try {
    config = gb::parse_driver_args(rest);
  } catch (const std::invalid_argument& e) {
    std::cerr << "sssp_bench: " << e.what() << "\n\n" << gb::driver_usage();
    return 2;
  }
  if (config.help) {
    std::cout << gb::driver_usage();
    return 0;
  }

  try {
    gb::Suite suite(config);  // graphs, CSR and sources: built here, before any timing
    if (config.verify_only) return suite.verify_all(std::cout) ? 0 : 1;

    suite.register_benchmarks();
    benchmark::AddCustomContext("sssp_version", std::string(sssp::version_string()));
    benchmark::AddCustomContext("sssp_paper", std::string(sssp::paper_arxiv_id()));
    benchmark::AddCustomContext("sssp_diagnostics_compiled",
                                sssp::bench::kDiagnosticsCompiled ? "true" : "false");
    benchmark::AddCustomContext("sssp_timed_region",
                                "solver call only (manual time); graph build, reference run and "
                                "verification excluded");
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    if (suite.any_failed()) {
      std::cerr << "sssp_bench: at least one case failed its correctness check (not timed)\n";
      return 1;
    }
    return 0;
  } catch (const std::invalid_argument& e) {
    std::cerr << "sssp_bench: " << e.what() << '\n';
    return 2;
  }
}
