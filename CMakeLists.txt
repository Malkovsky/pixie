cmake_minimum_required(VERSION 3.12)
project(pixie)

set(MARCH "native" CACHE STRING "march compilier flag")
set(CMAKE_CXX_FLAGS "-march=${MARCH} -O3 -DNDEBUG")
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
)
include(GoogleTest)
FetchContent_MakeAvailable(googletest)
FetchContent_Declare(
  googlebenchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG v1.9.4  # Replace with the desired version/tag
)
FetchContent_MakeAvailable(googlebenchmark)
add_executable(benchmarks
    src/benchmarks.cpp)
target_include_directories(benchmarks
    PUBLIC include
)
set_property(TARGET benchmarks PROPERTY CXX_STANDARD 20)

target_link_libraries(benchmarks
    benchmark
    benchmark_main)
add_test(
    NAME Benchmarks
    COMMAND bt_benchmarks --benchmark_out=bm_report.csv --benchmark_out_format=csv
)

add_executable(unittests
    src/unittests.cpp)
target_include_directories(unittests
    PUBLIC include
    PUBLIC ${GOOGLETEST_SOURCE_DIR}/include
)
set_property(TARGET unittests PROPERTY CXX_STANDARD 20)
target_link_libraries(unittests
    gtest
    gtest_main)

add_executable(bench_rmm
    src/bench_rmm.cpp)
target_include_directories(bench_rmm
    PUBLIC include
)
set_property(TARGET bench_rmm PROPERTY CXX_STANDARD 20)

add_executable(fuzzer
    src/fuzzer.cpp)
target_include_directories(fuzzer
    PUBLIC include
    PUBLIC include/misc
)
set_property(TARGET fuzzer PROPERTY CXX_STANDARD 20)

enable_testing()