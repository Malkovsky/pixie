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

option(PIXIE_DISABLE_BMI2 "Disable BMI2 instructions" OFF)
if(PIXIE_DISABLE_BMI2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-bmi2")
    add_compile_definitions(PIXIE_DISABLE_BMI2)
    message(STATUS "PIXIE_DISABLE_BMI2 is ON")
endif()

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
    if (CMAKE_VERSION VERSION_LESS 3.25)
        message(FATAL_ERROR
                "PIXIE_THIRD_PARTY_BACKENDS requires CMake 3.25 for pasta::bit_vector")
    endif ()
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
    set(immer_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(immer_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(immer_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(immer_BUILD_EXTRAS OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
            immer
            GIT_REPOSITORY https://github.com/arximboldi/immer.git
            GIT_TAG bd4fc749b97dfa2b66a8f3de00bbf234db4856ef # 2026-02-25
    )
    FetchContent_MakeAvailable(immer)

    set(PASTA_BIT_VECTOR_BUILD_TESTS OFF CACHE BOOL "Disable pasta::bit_vector tests" FORCE)
    set(PASTA_BIT_VECTOR_BUILD_BENCHMARKS OFF CACHE BOOL "Disable pasta::bit_vector benchmarks" FORCE)
    FetchContent_Declare(
            pasta_bit_vector
            GIT_REPOSITORY https://github.com/pasta-toolbox/bit_vector.git
            GIT_TAG 3ffb6e5a2e58c76425de8197bfe554eb1bf9fd94
    )
    FetchContent_MakeAvailable(pasta_bit_vector)
    target_compile_definitions(pasta_bit_vector INTERFACE PIXIE_PASTA_SUPPORT)

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

    # Fetch headers without importing upstream standalone harnesses.
    FetchContent_Declare(
            succinct_bitvec
            GIT_REPOSITORY https://github.com/ishantgarg2332/succinct-bitvec-cpp.git
            GIT_TAG e825c01655a08266869febe9af56fca277583e31
            SOURCE_SUBDIR pixie-source-only
    )
    FetchContent_MakeAvailable(succinct_bitvec)

    include(CheckCXXSourceCompiles)
    # Recheck when the configured architecture or BMI2 option changes.
    unset(PIXIE_RANK_SELECT_HAS_BMI CACHE)
    check_cxx_source_compiles("
        #if !defined(__BMI__) || !defined(__BMI2__)
        #error Rank/select comparisons require BMI and BMI2
        #endif
        int main() { return 0; }
    " PIXIE_RANK_SELECT_HAS_BMI)
    if (PIXIE_RANK_SELECT_HAS_BMI)
        # Upstream's CMake project builds its standalone research harness. Fetch
        # only the sources; Pixie's adapter supplies the required global state.
        FetchContent_Declare(
                three_star
                GIT_REPOSITORY https://github.com/nrilu/3star-rankselect.git
                GIT_TAG a69449585d1c760b2ed1af9b2761e3c05d86b247
                GIT_SUBMODULES ""
                SOURCE_SUBDIR pixie-source-only
        )
        FetchContent_MakeAvailable(three_star)

        FetchContent_Declare(
                orzo
                GIT_REPOSITORY https://github.com/hpcgroup/orzo.git
                GIT_TAG 02a7cd72fefdc03c1f402a9af3ff74289404c8bf
                SOURCE_SUBDIR pixie-source-only
        )
        FetchContent_MakeAvailable(orzo)
        set(PIXIE_ORZO_PREPARED_DIR "${CMAKE_CURRENT_BINARY_DIR}/third_party/orzo")
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/third_party/orzo/prepare.cmake"
                "${orzo_SOURCE_DIR}/include/orzo/orzo.h"
                "${orzo_SOURCE_DIR}/include/orzo/utils.h")
        execute_process(
                COMMAND ${CMAKE_COMMAND}
                -DORZO_SOURCE_DIR=${orzo_SOURCE_DIR}
                -DORZO_OUTPUT_DIR=${PIXIE_ORZO_PREPARED_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/third_party/orzo/prepare.cmake
                RESULT_VARIABLE orzo_prepare_result
        )
        if (NOT orzo_prepare_result EQUAL 0)
            message(FATAL_ERROR "Failed to prepare the pinned Orzo headers")
        endif ()
    else ()
        message(STATUS "3-star and Orzo comparisons disabled: configured target lacks BMI/BMI2")
    endif ()
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

# Keep Pixie's own targets warning-clean without imposing -Werror on consumers.
function (pixie_enable_project_warnings target)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(${target} PRIVATE
                -Wall
                -Wextra
                -Werror)
    endif ()
endfunction ()

# Shared configuration for optional rank/select tests and benchmarks.
if (PIXIE_THIRD_PARTY_BACKENDS)
    add_library(pixie_rank_select_comparisons INTERFACE)
    target_link_libraries(pixie_rank_select_comparisons INTERFACE pasta_bit_vector)
    target_compile_definitions(pixie_rank_select_comparisons INTERFACE
            PIXIE_SUCCINCT_BITVEC_SUPPORT)
    target_include_directories(pixie_rank_select_comparisons INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/benchmarks)
    target_include_directories(pixie_rank_select_comparisons SYSTEM INTERFACE
            ${succinct_bitvec_SOURCE_DIR}/include)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        # Pasta's debug-assertion helper becomes unused in optimized builds.
        target_compile_options(pixie_rank_select_comparisons INTERFACE
                -Wno-unused-function)
    endif ()
    if (PIXIE_RANK_SELECT_HAS_BMI)
        target_include_directories(pixie_rank_select_comparisons SYSTEM INTERFACE
                ${three_star_SOURCE_DIR}
                ${PIXIE_ORZO_PREPARED_DIR}/include)
        target_sources(pixie_rank_select_comparisons INTERFACE
                ${CMAKE_CURRENT_SOURCE_DIR}/src/benchmarks/three_star_adapter.cpp)
        target_compile_definitions(pixie_rank_select_comparisons INTERFACE
                PIXIE_ORZO_SUPPORT PIXIE_3STAR_SUPPORT
                ALIGNMENT=2097152
                CHECK_ALL_DENSE_SPARSE=0
                COUNT_HIT_TYPES=0
                CUSTOM_L0=2048
                POP=2
                PREFETCH_NT=1
                PREFETCH_T0=1
                REVERSE=2
                TREE_DIVISION_EXP=3
                TREE_EXP_LVL0=19
                TREE_TYPE=3)
    endif ()
endif ()

if (PIXIE_DIAGNOSTICS)
    target_compile_definitions(pixie INTERFACE PIXIE_DIAGNOSTICS)
    target_link_libraries(pixie INTERFACE spdlog::spdlog_header_only)
endif ()

if (PIXIE_THIRD_PARTY_BACKENDS)
    target_compile_definitions(pixie INTERFACE SDSL_SUPPORT)
    target_compile_definitions(pixie INTERFACE PIXIE_IMMER_SUPPORT)
    target_link_libraries(pixie INTERFACE immer)
    target_include_directories(pixie SYSTEM INTERFACE
            $<BUILD_INTERFACE:${immer_SOURCE_DIR}>)
    target_include_directories(pixie INTERFACE
            $<BUILD_INTERFACE:${sdsl_lite_SOURCE_DIR}/include>)
    target_link_libraries(pixie INTERFACE pasta_bit_vector)
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

    add_executable(rank_select_tests
            src/tests/rank_select_tests.cpp)
    target_include_directories(rank_select_tests
            PUBLIC include)
    target_link_libraries(rank_select_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_link_libraries(rank_select_tests pixie_rank_select_comparisons)
    endif ()

    if (PIXIE_THIRD_PARTY_BACKENDS AND PIXIE_RANK_SELECT_HAS_BMI)
        add_executable(three_star_benchmark_tests
                src/tests/three_star_benchmark_tests.cpp)
        target_include_directories(three_star_benchmark_tests
                PUBLIC include)
        target_link_libraries(three_star_benchmark_tests
                gtest
                gtest_main
                ${PIXIE_DIAGNOSTICS_LIBS})
        target_link_libraries(three_star_benchmark_tests
                pixie_rank_select_comparisons)
    endif ()

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

    add_executable(file_archive_tests
            src/tests/file_archive_tests.cpp)
    target_include_directories(file_archive_tests
            PUBLIC include)
    target_link_libraries(file_archive_tests
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

    add_executable(serialization_tests
            src/tests/serialization_tests.cpp)
    target_include_directories(serialization_tests
            PUBLIC include)
    target_link_libraries(serialization_tests
            gtest
            gtest_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(integer_vector_tests
            src/tests/integer_vector_tests.cpp)
    target_include_directories(integer_vector_tests
            PUBLIC include)
    target_link_libraries(integer_vector_tests
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

    add_executable(select512_experimental_tests
            src/tests/select512_experimental_tests.cpp)
    target_include_directories(select512_experimental_tests
            PUBLIC include)
    target_link_libraries(select512_experimental_tests
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

    add_executable(bit_sequence_tests
            src/tests/bit_sequence_tests.cpp
            src/tests/packed_sequence_tests.cpp
            src/tests/sequence_tree_kernels_tests.cpp)
    target_link_libraries(bit_sequence_tests
            pixie::pixie gtest_main)
    target_compile_definitions(bit_sequence_tests
            PRIVATE PIXIE_SEQUENCE_TREE_TESTING)

    foreach (family IN ITEMS permutation permutable_sequence)
        add_executable(${family}_tests src/tests/${family}_tests.cpp)
        target_link_libraries(${family}_tests PRIVATE pixie::pixie gtest_main)
        target_compile_definitions(${family}_tests PRIVATE PIXIE_SEQUENCE_TREE_TESTING)
    endforeach ()

    add_executable(sorted_sequence_tests src/tests/sorted_sequence_tests.cpp)
    target_link_libraries(sorted_sequence_tests PRIVATE pixie::pixie gtest_main)
    target_compile_definitions(sorted_sequence_tests PRIVATE PIXIE_SEQUENCE_TREE_TESTING)

    set(PIXIE_TEST_TARGETS
            permutation_tests
            permutable_sequence_tests
            sorted_sequence_tests
            bit_sequence_tests
            bit_algorithms_unittests
            rank_select_tests
            benchmark_tests
            test_rmm
            tree_tests
            wavelet_tree_tests
            file_archive_tests
            storage_tests
            serialization_tests
            integer_vector_tests
            excess_positions_tests
            select512_experimental_tests
            excess_record_lows_tests
            rmq_tests)
    if (TARGET three_star_benchmark_tests)
        list(APPEND PIXIE_TEST_TARGETS three_star_benchmark_tests)
    endif ()
    foreach (test_target IN LISTS PIXIE_TEST_TARGETS)
        pixie_enable_project_warnings(${test_target})
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
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_link_libraries(rank_select_benchmarks pixie_rank_select_comparisons)
    endif ()

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

    add_executable(file_archive_benchmarks
            src/benchmarks/file_archive_benchmarks.cpp)
    target_include_directories(file_archive_benchmarks
            PUBLIC include)
    target_link_libraries(file_archive_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(serialization_benchmarks
            src/benchmarks/serialization_benchmarks.cpp)
    target_include_directories(serialization_benchmarks
            PUBLIC include)
    target_link_libraries(serialization_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})
    if (PIXIE_THIRD_PARTY_BACKENDS)
        target_include_directories(serialization_benchmarks
                PRIVATE ${sdsl_lite_SOURCE_DIR}/include)
        target_link_libraries(serialization_benchmarks pasta_bit_vector)
        if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
            target_compile_options(serialization_benchmarks
                    PRIVATE -Wno-unused-function)
        endif ()
    endif ()

    add_executable(integer_vector_benchmarks
            src/benchmarks/integer_vector_benchmarks.cpp)
    target_include_directories(integer_vector_benchmarks
            PUBLIC include)
    target_link_libraries(integer_vector_benchmarks
            benchmark
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

    add_executable(select512_benchmarks
            src/benchmarks/select512_benchmarks.cpp)
    target_include_directories(select512_benchmarks
            PUBLIC include)
    target_link_libraries(select512_benchmarks
            benchmark
            benchmark_main
            ${PIXIE_DIAGNOSTICS_LIBS})

    add_executable(bit_sequence_benchmarks
            src/benchmarks/bit_sequence_benchmarks.cpp
            src/benchmarks/sequence_controls_benchmarks.cpp
            src/benchmarks/sequence_tree_kernels_benchmarks.cpp)
    target_link_libraries(bit_sequence_benchmarks
            pixie::pixie benchmark benchmark_main)

    foreach (family IN ITEMS permutation permutable_sequence)
        add_executable(${family}_benchmarks src/benchmarks/${family}_benchmarks.cpp)
        target_link_libraries(${family}_benchmarks PRIVATE pixie::pixie benchmark benchmark_main)
    endforeach ()

    add_executable(sorted_sequence_benchmarks
            src/benchmarks/sorted_sequence_benchmarks.cpp)
    target_link_libraries(sorted_sequence_benchmarks
            PRIVATE pixie::pixie benchmark benchmark_main)

    set(PIXIE_BENCHMARK_TARGETS
            permutation_benchmarks
            permutable_sequence_benchmarks
            sorted_sequence_benchmarks
            bit_sequence_benchmarks
            rank_select_benchmarks
            rmm_benchmarks
            rmm_btree_benchmarks
            rmq_benchmarks
            louds_tree_benchmarks
            wavelet_tree_benchmarks
            file_archive_benchmarks
            integer_vector_benchmarks
            serialization_benchmarks
            bp_tree_benchmarks
            dfuds_tree_benchmarks
            alignment_comparison_benchmarks
            excess_positions_benchmarks
            select512_benchmarks)
    if (PIXIE_THIRD_PARTY_BACKENDS)
        list(APPEND PIXIE_BENCHMARK_TARGETS
                rmm_sdsl_benchmarks)
    endif ()
    foreach (benchmark_target IN LISTS PIXIE_BENCHMARK_TARGETS)
        pixie_enable_project_warnings(${benchmark_target})
    endforeach ()
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
