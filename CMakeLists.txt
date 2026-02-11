cmake_minimum_required(VERSION 3.18)
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

# ---------------------------------------------------------------------------
# Build options
# ---------------------------------------------------------------------------
option(PIXIE_TESTS "Build unit tests" ON)
option(PIXIE_BENCHMARKS "Build benchmarks (includes comparison benchmarks against third-party libraries)" OFF)
option(PIXIE_DOCS "Build Doxygen documentation" OFF)

# ---------------------------------------------------------------------------
# Dependencies (fetched only when needed)
# ---------------------------------------------------------------------------
include(FetchContent)

if(PIXIE_BENCHMARKS)
    FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.9.4
    )
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable Google Benchmark tests")
    FetchContent_MakeAvailable(googlebenchmark)

    # Third-party libraries for comparison benchmarks
    FetchContent_Declare(
        pasta_bit_vector
        GIT_REPOSITORY https://github.com/pasta-toolbox/bit_vector.git
        GIT_TAG origin/main
    )
    FetchContent_MakeAvailable(pasta_bit_vector)

    if(NOT WIN32)
        # Disable sdsl-lite's own tests/examples/tutorials
        set(SDSL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
        set(SDSL_BUILD_TUTORIAL OFF CACHE BOOL "" FORCE)
        set(SDSL_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        FetchContent_Declare(
            sdsl_lite
            GIT_REPOSITORY https://github.com/xxsds/sdsl-lite.git
            GIT_TAG v3.0.3
        )
        FetchContent_MakeAvailable(sdsl_lite)
    endif()
endif()

if(PIXIE_TESTS)
    if(NOT TARGET gtest_main)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.17.0
        )
        FetchContent_MakeAvailable(googletest)
    endif()
    include(GoogleTest)
endif()

# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------
if(PIXIE_TESTS)
    enable_testing()

    add_executable(unittests
        src/tests/unittests.cpp)
    target_include_directories(unittests
        PUBLIC include)
    target_link_libraries(unittests
        gtest
        gtest_main)

    add_executable(benchmark_tests
        src/tests/benchmark_tests.cpp)
    target_include_directories(benchmark_tests
        PUBLIC include)
    target_link_libraries(benchmark_tests
        gtest
        gtest_main)

    add_executable(test_rmm
        src/tests/test_rmm.cpp)
    target_include_directories(test_rmm
        PUBLIC include
        PUBLIC include/misc)
    target_link_libraries(test_rmm
        gtest
        gtest_main)

    add_executable(louds_tree_tests
        src/tests/louds_tree_tests.cpp)
    target_include_directories(louds_tree_tests
        PUBLIC include)
    target_link_libraries(louds_tree_tests
        gtest
        gtest_main)
endif()

# ---------------------------------------------------------------------------
# Benchmarks (Pixie-only and comparison benchmarks)
# ---------------------------------------------------------------------------
if(PIXIE_BENCHMARKS)
    add_executable(benchmarks
        src/benchmarks/benchmarks.cpp)
    target_include_directories(benchmarks
        PUBLIC include)
    target_link_libraries(benchmarks
        benchmark
        benchmark_main
        pasta_bit_vector)
    target_compile_definitions(benchmarks PRIVATE PIXIE_THIRD_PARTY_BENCHMARKS)

    add_executable(bench_rmm
        src/benchmarks/bench_rmm.cpp)
    target_include_directories(bench_rmm
        PUBLIC include)
    target_link_libraries(bench_rmm
        benchmark)

    add_executable(bench_rmm_sdsl
        src/benchmarks/bench_rmm_sdsl.cpp)
    target_include_directories(bench_rmm_sdsl
        PUBLIC include
        PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
    target_link_libraries(bench_rmm_sdsl
        PRIVATE
            benchmark
            sdsl-lite)

    add_executable(louds_tree_benchmarks
        src/benchmarks/louds_tree_benchmarks.cpp)
    target_include_directories(louds_tree_benchmarks
        PUBLIC include)
    target_link_libraries(louds_tree_benchmarks
        benchmark
        benchmark_main)

    add_executable(alignment_comparison
        src/benchmarks/alignment_comparison.cpp)
    target_include_directories(alignment_comparison
        PUBLIC include)
    target_link_libraries(alignment_comparison
        benchmark
        benchmark_main)
endif()

# ---------------------------------------------------------------------------
# Documentation (Doxygen)
# ---------------------------------------------------------------------------
if(PIXIE_DOCS)
    find_package(Doxygen REQUIRED)

    FetchContent_Declare(
        doxygen-awesome-css
        URL https://github.com/jothepro/doxygen-awesome-css/archive/refs/heads/main.zip
    )
    FetchContent_MakeAvailable(doxygen-awesome-css)

    FetchContent_GetProperties(doxygen-awesome-css SOURCE_DIR AWESOME_CSS_DIR)

    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/src/docs/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile)
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
