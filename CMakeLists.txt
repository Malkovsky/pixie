cmake_minimum_required(VERSION 3.12)
project(pixie)


set(MARCH "native" CACHE STRING "march compilier flag")
set(CMAKE_CXX_FLAGS "-march=${MARCH}")
message(STATUS "MARCH is '${MARCH}'")

option(USE_NO_AVX512 "Disable AVX512 instructions" OFF)
if(USE_NO_AVX512)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-avx512f")
    message(STATUS "USE_NO_AVX512 is ON")
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

add_executable(benchmarkstests
    src/benchmarkstests.cpp)
target_include_directories(benchmarkstests
    PUBLIC include
    PUBLIC ${GOOGLETEST_SOURCE_DIR}/include
)
set_property(TARGET benchmarkstests PROPERTY CXX_STANDARD 20)
target_link_libraries(benchmarkstests
    gtest
    gtest_main)

enable_testing()