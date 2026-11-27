cmake_minimum_required(VERSION 3.24)
project(almost_linear_sssp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(SSSP_BUILD_TESTS "Build tests" ON)
option(SSSP_BUILD_BENCHMARKS "Build benchmarks" ON)
option(SSSP_WERROR "Treat warnings as errors" OFF)
option(SSSP_SANITIZE "Enable AddressSanitizer + UndefinedBehaviorSanitizer" OFF)
option(SSSP_COVERAGE "Enable source-based coverage (clang)" OFF)
option(SSSP_DIAGNOSTICS "Compile diagnostic counters (operation counts, fallback stats)" ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# Architecture is detected, never hard-coded (Apple Silicon / Intel / Linux x86_64 / aarch64).
message(STATUS "SSSP: system processor = ${CMAKE_SYSTEM_PROCESSOR}, compiler = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

include(cmake/Warnings.cmake)
include(cmake/Dependencies.cmake)

# ---------------------------------------------------------------- library
# Sources are globbed per directory so that parallel work on separate modules
# never needs to edit this file.
file(GLOB_RECURSE SSSP_SOURCES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(sssp ${SSSP_SOURCES})
target_include_directories(sssp PUBLIC ${PROJECT_SOURCE_DIR}/include PRIVATE ${PROJECT_SOURCE_DIR}/src)
sssp_set_warnings(sssp)
target_compile_definitions(sssp PUBLIC SSSP_DIAGNOSTICS=$<BOOL:${SSSP_DIAGNOSTICS}>)

# ---------------------------------------------------------------- tools
file(GLOB SSSP_TOOL_SOURCES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/tools/*.cpp)
foreach(tool_src IN LISTS SSSP_TOOL_SOURCES)
  get_filename_component(tool_name ${tool_src} NAME_WE)
  add_executable(${tool_name} ${tool_src})
  target_link_libraries(${tool_name} PRIVATE sssp)
  sssp_set_warnings(${tool_name})
endforeach()

if(SSSP_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(SSSP_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()
