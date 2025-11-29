cmake_minimum_required(VERSION 3.12)
project(pixie)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(MARCH "native" CACHE STRING "march compilier flag")
set(CMAKE_CXX_FLAGS "-march=${MARCH}")
message(STATUS "MARCH is '${MARCH}'")

option(DISABLE_AVX512 "Disable AVX512 instructions" OFF)
if(DISABLE_AVX512)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-avx512f")
    message(STATUS "DISABLE_AVX512 is ON")
endif()

option(ENABLE_ADDRESS_SANITIZER "Enable AddressSanitizer" OFF)
if(ENABLE_ADDRESS_SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
    message(STATUS "AddressSanitizer is ON")
endif()

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
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable Google Benchmark tests")
FetchContent_MakeAvailable(googlebenchmark)

FetchContent_Declare(
  pasta_bit_vector
  GIT_REPOSITORY https://github.com/pasta-toolbox/bit_vector.git
  GIT_TAG origin/main
)
FetchContent_MakeAvailable(pasta_bit_vector)

# add_subdirectory(third_party/pasta_bv)

add_executable(benchmarks
    src/benchmarks.cpp)
target_include_directories(benchmarks
    PUBLIC include
)

target_link_libraries(benchmarks
    benchmark
    benchmark_main
    pasta_bit_vector)
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
target_link_libraries(unittests
    gtest
    gtest_main)

add_executable(benchmark_tests
    src/benchmark_tests.cpp)
target_include_directories(benchmark_tests
    PUBLIC include
    PUBLIC ${GOOGLETEST_SOURCE_DIR}/include
)
set_property(TARGET benchmark_tests PROPERTY CXX_STANDARD 20)
target_link_libraries(benchmark_tests
    gtest
    gtest_main)

add_executable(bench_rmm
    src/bench_rmm.cpp)
target_include_directories(bench_rmm
    PUBLIC include
)
target_link_libraries(bench_rmm
    benchmark)

add_executable(test_rmm
    src/test_rmm.cpp)
target_include_directories(test_rmm
    PUBLIC include
    PUBLIC include/misc
)
target_link_libraries(test_rmm
    gtest
    gtest_main)

enable_testing()