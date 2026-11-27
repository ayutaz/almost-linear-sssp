include(FetchContent)
include(CheckCXXCompilerFlag)
# Versions are pinned by tag + SHA-256. Do not bump without an ADR.
if(SSSP_BUILD_TESTS)
  # TEST_CASE/SECTION expand Catch2's unique-name macro inside OUR test sources. With the
  # default it uses __COUNTER__, which clang >= 22 diagnoses under -Wpedantic as a C2y
  # extension, so -Werror fails on our own code (marking the headers SYSTEM does not help: the
  # token is expanded in our files). Make Catch2 use the standard __LINE__ instead; the only
  # constraint is not to put two test cases or sections on one source line.
  set(CATCH_CONFIG_NO_COUNTER ON CACHE BOOL "" FORCE)
  FetchContent_Declare(Catch2
    URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.8.1.tar.gz
    URL_HASH SHA256=18b3f70ac80fccc340d8c6ff0f339b2ae64944782f8d2fca2bd705cf47cadb79
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()
if(SSSP_BUILD_BENCHMARKS)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  # Third-party warnings must not fail our build: Google Benchmark's own -Werror is turned off
  # for the dependency only. Our targets keep -Werror (SSSP_WERROR, cmake/Warnings.cmake).
  set(BENCHMARK_ENABLE_WERROR OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(benchmark
    URL https://github.com/google/benchmark/archive/refs/tags/v1.9.4.tar.gz
    URL_HASH SHA256=b334658edd35efcf06a99d9be21e4e93e092bd5f95074c1673d5c8705d95c104
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(benchmark)
  # v1.9.4 also adds -pedantic-errors unconditionally, and clang >= 22 diagnoses the
  # __COUNTER__ probe in benchmark.h as a C2y extension (-Wc2y-extensions, fixed upstream in
  # v1.9.5 with a pragma). Keep that one diagnostic a warning in the dependency's own build.
  check_cxx_compiler_flag(-Wc2y-extensions SSSP_CXX_HAS_WC2Y_EXTENSIONS)
  if(SSSP_CXX_HAS_WC2Y_EXTENSIONS)
    foreach(_bm_target IN ITEMS benchmark benchmark_main)
      target_compile_options(${_bm_target} PRIVATE -Wno-error=c2y-extensions)
    endforeach()
  endif()
  # Consumers see the dependency's headers as system headers (-isystem), so warnings inside
  # benchmark.h are not reported against our sources. (FetchContent_Declare(... SYSTEM) needs
  # CMake 3.25; the project minimum is 3.24.)
  foreach(_bm_target IN ITEMS benchmark benchmark_main)
    get_target_property(_bm_inc ${_bm_target} INTERFACE_INCLUDE_DIRECTORIES)
    if(_bm_inc)
      set_property(TARGET ${_bm_target} APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${_bm_inc})
    endif()
  endforeach()
endif()
