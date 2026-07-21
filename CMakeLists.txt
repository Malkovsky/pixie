cmake_minimum_required(VERSION 3.18)
project(pixie)

set(PIXIE_IS_TOP_LEVEL OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(PIXIE_IS_TOP_LEVEL ON)
endif ()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(MARCH "native" CACHE STRING "march compilier flag")
set(CMAKE_CXX_FLAGS "-march=${MARCH}")
message(STATUS "MARCH is '${MARCH}'")

option(DISABLE_AVX512 "Disable AVX512 instructions" OFF)
if (DISABLE_AVX512)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-avx512f")
    message(STATUS "DISABLE_AVX512 is ON")
endif ()

option(ENABLE_ADDRESS_SANITIZER "Enable AddressSanitizer" OFF)
if (ENABLE_ADDRESS_SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
    message(STATUS "AddressSanitizer is ON")
endif ()

option(PIXIE_COVERAGE "Enable coverage instrumentation" OFF)
if (PIXIE_COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
    message(STATUS "Coverage instrumentation is ON")
endif ()

# ---------------------------------------------------------------------------
# Build options
# ---------------------------------------------------------------------------
option(PIXIE_TESTS "Build unit tests" ${PIXIE_IS_TOP_LEVEL})
option(PIXIE_BENCHMARKS "Build benchmarks" OFF)
option(PIXIE_THIRD_PARTY_BACKENDS "Build optional third-party backend integrations" OFF)
option(PIXIE_DIAGNOSTICS "Include diagnostic logs" OFF)
option(PIXIE_DOCS "Build Doxygen documentation" OFF)

if (PIXIE_THIRD_PARTY_BACKENDS)
    add_compile_definitions(SDSL_SUPPORT)
endif ()

if (PIXIE_DIAGNOSTICS)
    add_compile_definitions(PIXIE_DIAGNOSTICS)
    set(PIXIE_DIAGNOSTICS_LIBS spdlog::spdlog_header_only)
endif ()

# ---------------------------------------------------------------------------
# Dependencies (fetched only when needed)
# ---------------------------------------------------------------------------
include(FetchContent)

if (PIXIE_DIAGNOSTICS)
    set(SPDLOG_BUILD_SHARED OFF CACHE BOOL "" FORCE)
    set(SPDLOG_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE)
    set(SPDLOG_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(SPDLOG_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
            spdlog
            GIT_REPOSITORY https://github.com/gabime/spdlog.git
            GIT_TAG v1.14.1
    )
    FetchContent_MakeAvailable(spdlog)
endif ()

if (PIXIE_BENCHMARKS)
    FetchContent_Declare(
            googlebenchmark
            GIT_REPOSITORY https://github.com/google/benchmark.git
            GIT_TAG v1.9.4
    )
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable Google Benchmark tests" FORCE)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Disable Google Benchmark tests" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Disable Google Benchmark install targets" FORCE)
    set(BENCHMARK_INSTALL_DOCS OFF CACHE BOOL "Disable Google Benchmark install docs" FORCE)
    FetchContent_MakeAvailable(googlebenchmark)
endif ()

if (PIXIE_THIRD_PARTY_BACKENDS)
    set(PASTA_BIT_VECTOR_BUILD_TESTS OFF CACHE BOOL "Disable pasta::bit_vector tests" FORCE)
    set(PASTA_BIT_VECTOR_BUILD_BENCHMARKS OFF CACHE BOOL "Disable pasta::bit_vector benchmarks" FORCE)
    FetchContent_Declare(
            pasta_bit_vector
            GIT_REPOSITORY https://github.com/pasta-toolbox/bit_vector.git
            GIT_TAG origin/main
    )
    FetchContent_MakeAvailable(pasta_bit_vector)

    # Disable sdsl-lite's own tests/examples/tutorials
    set(SDSL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(SDSL_BUILD_TUTORIAL OFF CACHE BOOL "" FORCE)
    set(SDSL_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
            sdsl_lite
            GIT_REPOSITORY https://github.com/xxsds/sdsl-lite.git
            GIT_TAG v3.0.3
            PATCH_COMMAND ${CMAKE_COMMAND}
            -DROOT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DSDSL_SOURCE_DIR=<SOURCE_DIR>
            -P ${CMAKE_CURRENT_SOURCE_DIR}/third_party/sdsl-lite/patch.cmake
    )
    FetchContent_MakeAvailable(sdsl_lite)
endif ()

if (PIXIE_TESTS)
    if (NOT TARGET gtest_main)
        set(BUILD_GMOCK OFF CACHE BOOL "Disable GoogleMock" FORCE)
        set(INSTALL_GTEST OFF CACHE BOOL "Disable GoogleTest install targets" FORCE)
        FetchContent_Declare(
                googletest
                GIT_REPOSITORY https://github.com/google/googletest.git
                GIT_TAG v1.17.0
        )
        FetchContent_MakeAvailable(googletest)
    endif ()
    include(GoogleTest)
endif ()

# Downstream consumers link this target; Pixie itself remains header-only.
add_library(pixie INTERFACE)
add_library(pixie::pixie ALIAS pixie)
target_compile_features(pixie INTERFACE cxx_std_20)
target_include_directories(pixie INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

if (PIXIE_DIAGNOSTICS)
    target_compile_definitions(pixie INTERFACE PIXIE_DIAGNOSTICS)
    target_link_libraries(pixie INTERFACE spdlog::spdlog_header_only)
endif ()

if (PIXIE_THIRD_PARTY_BACKENDS)
    target_compile_definitions(pixie INTERFACE SDSL_SUPPORT)
    target_include_directories(pixie INTERFACE
            $<BUILD_INTERFACE:${sdsl_lite_SOURCE_DIR}/include>)
endif ()

# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------
if (PIXIE_TESTS)
    enable_testing()

    add_executable(bit_algorithms_unittests
            src/tests/bits_unittests.cc)
    target_include_directories(bit_algorithms_unittests
            PUBLIC include)
    target_link_libraries(bit_algorithms_unittests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(rank_select_unittests
            src/tests/rank_select_unittests.cc)
    target_include_directories(rank_select_unittests
            PUBLIC include)
    target_link_libraries(rank_select_unittests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(rank_select_tests
            src/tests/rank_select_tests.cpp)
    target_include_directories(rank_select_tests
            PUBLIC include)
    target_link_libraries(rank_select_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(benchmark_tests
            src/tests/benchmark_tests.cpp)
    target_include_directories(benchmark_tests
            PUBLIC include)
    target_link_libraries(benchmark_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(test_rmm
            src/tests/test_rmm.cpp)
    target_include_directories(test_rmm
            PUBLIC include
            PUBLIC include/misc)
    target_link_libraries(test_rmm
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_include_directories(test_rmm
                PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
    endif ()

    add_executable(tree_tests
            src/tests/tree_tests.cpp)
    target_include_directories(tree_tests
            PUBLIC include)
    target_link_libraries(tree_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(wavelet_tree_tests
            src/tests/wavelet_tree_tests.cpp)
    target_include_directories(wavelet_tree_tests
            PUBLIC include)
    target_link_libraries(wavelet_tree_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(storage_tests
            src/tests/storage_tests.cpp)
    target_include_directories(storage_tests
            PUBLIC include)
    target_link_libraries(storage_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(excess_positions_tests
            src/tests/excess_positions_tests.cpp)
    target_include_directories(excess_positions_tests
            PUBLIC include)
    target_link_libraries(excess_positions_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(excess_record_lows_tests
            src/tests/excess_record_lows_tests.cpp)
    target_include_directories(excess_record_lows_tests
            PUBLIC include)
    target_link_libraries(excess_record_lows_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(rmq_tests
            src/tests/rmq_tests.cpp)
    target_include_directories(rmq_tests
            PUBLIC include)
    target_link_libraries(rmq_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_include_directories(rmq_tests
                PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
    endif ()

    add_executable(pivco_huffman_tests
            src/tests/pivco_huffman_tests.cpp)
    target_include_directories(pivco_huffman_tests
            PUBLIC include)
    target_link_libraries(pivco_huffman_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    set(PIXIE_TEST_TARGETS
            bit_algorithms_unittests
            rank_select_unittests
            rank_select_tests
            benchmark_tests
            test_rmm
            tree_tests
            wavelet_tree_tests
            storage_tests
            excess_positions_tests
            excess_record_lows_tests
            rmq_tests
            pivco_huffman_tests)
    foreach (test_target IN LISTS PIXIE_TEST_TARGETS)
        gtest_discover_tests(${test_target}
                DISCOVERY_MODE PRE_TEST
                TEST_PREFIX "${test_target}."
                PROPERTIES LABELS "${test_target}")
    endforeach ()
endif ()

# ---------------------------------------------------------------------------
# Benchmarks (Pixie-only and comparison benchmarks)
# ---------------------------------------------------------------------------
if (PIXIE_BENCHMARKS)
    add_executable(rank_select_benchmarks
            src/benchmarks/rank_select_benchmarks.cpp)
    target_include_directories(rank_select_benchmarks
            PUBLIC include)
    target_link_libraries(rank_select_benchmarks
            benchmark
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(rmm_benchmarks
            src/benchmarks/rmm_benchmarks.cpp)
    target_include_directories(rmm_benchmarks
            PUBLIC include)
    target_link_libraries(rmm_benchmarks
            benchmark
            ${PIXIE_DIAGNOSTICS_LIBS})
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_include_directories(rmm_benchmarks
                PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
    endif ()

    add_executable(rmm_btree_benchmarks
            src/benchmarks/rmm_btree_benchmarks.cpp)
    target_include_directories(rmm_btree_benchmarks
            PUBLIC include)
    target_link_libraries(rmm_btree_benchmarks
            benchmark
            ${PIXIE_DIAGNOSTICS_LIBS})
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_include_directories(rmm_btree_benchmarks
                PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
    endif ()

    add_executable(rmq_benchmarks
            src/benchmarks/rmq_benchmarks.cpp)
    target_include_directories(rmq_benchmarks
            PUBLIC include)
    target_link_libraries(rmq_benchmarks
            benchmark
            ${PIXIE_DIAGNOSTICS_LIBS})
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_include_directories(rmq_benchmarks
                PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
        target_compile_definitions(rmq_benchmarks
                PRIVATE PIXIE_THIRD_PARTY_BENCHMARKS)
    endif ()

    if (PIXIE_THIRD_PARTY_BACKENDS)
        add_executable(rmm_sdsl_benchmarks
                src/benchmarks/rmm_sdsl_benchmarks.cpp)
        target_include_directories(rmm_sdsl_benchmarks
                PUBLIC include
                PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
        target_link_libraries(rmm_sdsl_benchmarks
                PRIVATE
                benchmark
                ${PIXIE_DIAGNOSTICS_LIBS})
    endif ()

    add_executable(louds_tree_benchmarks
            src/benchmarks/louds_tree_benchmarks.cpp)
    target_include_directories(louds_tree_benchmarks
            PUBLIC include)
    target_link_libraries(louds_tree_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(wavelet_tree_benchmarks
            src/benchmarks/wavelet_tree_benchmarks.cpp)
    target_include_directories(wavelet_tree_benchmarks
            PUBLIC include)
    target_link_libraries(wavelet_tree_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(bp_tree_benchmarks
            src/benchmarks/bp_tree_benchmarks.cpp)
    target_include_directories(bp_tree_benchmarks
            PUBLIC include)
    target_link_libraries(bp_tree_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_include_directories(bp_tree_benchmarks
                PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
    endif ()

    add_executable(dfuds_tree_benchmarks
            src/benchmarks/dfuds_tree_benchmarks.cpp)
    target_include_directories(dfuds_tree_benchmarks
            PUBLIC include)
    target_link_libraries(dfuds_tree_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_include_directories(dfuds_tree_benchmarks
                PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
    endif ()

    add_executable(alignment_comparison_benchmarks
            src/benchmarks/alignment_comparison_benchmarks.cpp)
    target_include_directories(alignment_comparison_benchmarks
            PUBLIC include)
    target_link_libraries(alignment_comparison_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(excess_positions_benchmarks
            src/benchmarks/excess_positions_benchmarks.cpp)
    target_include_directories(excess_positions_benchmarks
            PUBLIC include)
    target_link_libraries(excess_positions_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(pivco_huffman_benchmarks
            src/benchmarks/pivco_huffman_benchmarks.cpp)
    target_include_directories(pivco_huffman_benchmarks
            PUBLIC include)
    target_link_libraries(pivco_huffman_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})
endif ()

# ---------------------------------------------------------------------------
# Documentation (Doxygen)
# ---------------------------------------------------------------------------
if (PIXIE_DOCS)
    find_package(Doxygen REQUIRED)

    FetchContent_Declare(
            doxygen-awesome-css
            URL https://github.com/jothepro/doxygen-awesome-css/archive/refs/heads/main.zip
    )
    FetchContent_MakeAvailable(doxygen-awesome-css)

    FetchContent_GetProperties(doxygen-awesome-css SOURCE_DIR AWESOME_CSS_DIR)

    set(PIXIE_DOC_INPUT
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/docs/benchmark_results.md")
        list(APPEND PIXIE_DOC_INPUT
                "${CMAKE_CURRENT_SOURCE_DIR}/src/docs/benchmark_results.md")
    endif ()
    string(JOIN " \\\n                         " PIXIE_DOC_INPUT
            ${PIXIE_DOC_INPUT})

    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/src/docs/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile)
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    add_custom_target(docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
endif ()
