# Google Benchmark drivers (P2-BENCH-CPP). docs/design/benchmark_protocol.md governs how they are
# run; the JSONL pipeline (tools/sssp_bench_runner, scripts/run_benchmarks.sh) is separate.
#
# Warnings: every target defined here is ours and gets sssp_set_warnings (-Wall -Wextra
# -Wpedantic ... and -Werror when SSSP_WERROR=ON, -fno-fast-math always). Google Benchmark's
# headers reach us through -isystem (INTERFACE_SYSTEM_INCLUDE_DIRECTORIES, set in
# cmake/Dependencies.cmake), so warnings inside benchmark.h are not reported against our sources.
# Our sources register benchmarks at run time (benchmark::RegisterBenchmark) and never use the
# BENCHMARK* macros: those expand __COUNTER__ in OUR translation unit, which clang >= 22 reports
# under -Wpedantic as a C2y extension (-Wc2y-extensions) and -Werror would turn into an error
# that the SYSTEM include cannot suppress (the same issue as Catch2's TEST_CASE, P0-CI-FIX-2).

# Driver library: case enumeration, verification and registration. Shared by the executable and
# its tests so that the tests exercise the code the executable runs.
add_library(sssp_bench_driver STATIC ${CMAKE_CURRENT_SOURCE_DIR}/sssp_bench_driver.cpp)
target_include_directories(sssp_bench_driver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(sssp_bench_driver PUBLIC sssp benchmark::benchmark)
sssp_set_warnings(sssp_bench_driver)

add_executable(sssp_bench ${CMAKE_CURRENT_SOURCE_DIR}/sssp_bench_main.cpp)
target_link_libraries(sssp_bench PRIVATE sssp_bench_driver)
sssp_set_warnings(sssp_bench)

if(SSSP_BUILD_TESTS)
  include(Catch)
  add_executable(sssp_bench_driver_tests ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_sssp_bench_driver.cpp)
  target_link_libraries(sssp_bench_driver_tests PRIVATE sssp_bench_driver Catch2::Catch2WithMain)
  sssp_set_warnings(sssp_bench_driver_tests)
  # Label "unit" so that the per-label CTest steps of CI run these tests.
  catch_discover_tests(sssp_bench_driver_tests TEST_PREFIX "unit:bench-driver:"
                       PROPERTIES LABELS unit DISCOVERY_MODE PRE_TEST)

  # The executable itself: every registered algorithm on small shared inputs, checked against
  # binary-heap Dijkstra and verify_sssp. Nothing is timed (--sssp_verify_only); exit 1 on any
  # mismatch.
  add_test(NAME "integration:sssp_bench verify-only"
           COMMAND sssp_bench --sssp_verify_only
                   --sssp_graphs=er,grid,road,tree,ba,regular,chung_lu,adversarial_long_path
                   --sssp_sizes=60,300 --sssp_weights=uniform_int,wide_range,uniform_real
                   --sssp_sources=2)
  set_tests_properties("integration:sssp_bench verify-only" PROPERTIES LABELS integration)
endif()
