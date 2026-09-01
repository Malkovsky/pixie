# Compiled Apache-2.0 PivCo-Huffman backend.
#
# The codec is compiled once per available architecture tier. The facade in
# include/pixie/huffman/pivco_huffman.h remains a thin C++ CRTP adapter.

set(PIXIE_PIVCO_ROOT
        "${CMAKE_CURRENT_SOURCE_DIR}/third_party/pivco")

if (NOT EXISTS "${PIXIE_PIVCO_ROOT}/LICENSE")
    message(FATAL_ERROR
            "The PivCo-Huffman sources and Apache-2.0 LICENSE are required")
endif ()

set(PIXIE_PIVCO_INCLUDE
        "${PIXIE_PIVCO_ROOT}/include")
set(PIXIE_PIVCO_SRC
        "${PIXIE_PIVCO_ROOT}/src")

set(PIXIE_PIVCO_COMMON_SOURCES
        "${PIXIE_PIVCO_SRC}/huffman_table.c"
        "${PIXIE_PIVCO_SRC}/joint_lengths.c"
        "${PIXIE_PIVCO_SRC}/pivco_huffman.c"
        "${PIXIE_PIVCO_SRC}/pivcohuf_file.c")

add_library(pixie_pivco_scalar OBJECT
        "${PIXIE_PIVCO_SRC}/pivco_huffman_codec.c")
target_compile_definitions(pixie_pivco_scalar
        PRIVATE PIVCO_BACKEND_SCALAR=1)

set(PIXIE_PIVCO_OBJECTS
        $<TARGET_OBJECTS:pixie_pivco_scalar>)
set(PIXIE_PIVCO_DEFINITIONS)
set(PIXIE_PIVCO_OPTIONS)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" PIXIE_PIVCO_PROCESSOR)

if (PIXIE_PIVCO_PROCESSOR MATCHES "^(x86_64|amd64)$")
    list(APPEND PIXIE_PIVCO_DEFINITIONS
            PIVCO_HAS_AVX2=1
            PIVCO_HAS_SSE4=1)
    list(APPEND PIXIE_PIVCO_OPTIONS
            -mavx2
            -mbmi2
            -msse4.1
            -mpopcnt)

    add_library(pixie_pivco_x86 OBJECT
            "${PIXIE_PIVCO_SRC}/pivco_huffman_codec.c")
    target_compile_definitions(pixie_pivco_x86
            PRIVATE PIVCO_BACKEND_X86=1)
    list(APPEND PIXIE_PIVCO_OBJECTS
            $<TARGET_OBJECTS:pixie_pivco_x86>)
    list(APPEND PIXIE_PIVCO_COMMON_SOURCES
            "${PIXIE_PIVCO_SRC}/pivco_huffman_x86_tables.c")

    if (NOT DISABLE_AVX512 AND EXISTS "/proc/cpuinfo")
        file(READ "/proc/cpuinfo" PIXIE_PIVCO_CPUINFO)
        string(FIND "${PIXIE_PIVCO_CPUINFO}" "avx512_vbmi2"
                PIXIE_PIVCO_AVX512_POSITION)
        include(CheckCCompilerFlag)
        check_c_compiler_flag("-mavx512vbmi2" PIXIE_PIVCO_COMPILER_HAS_AVX512)
        if (PIXIE_PIVCO_COMPILER_HAS_AVX512 AND
                NOT PIXIE_PIVCO_AVX512_POSITION EQUAL -1)
            list(APPEND PIXIE_PIVCO_DEFINITIONS
                    PIVCO_HAS_AVX512=1)
            list(APPEND PIXIE_PIVCO_OPTIONS
                    -mavx512f
                    -mavx512bw
                    -mavx512vl
                    -mavx512vbmi
                    -mavx512vbmi2
                    -mavx512vpopcntdq)
            add_library(pixie_pivco_avx512 OBJECT
                    "${PIXIE_PIVCO_SRC}/pivco_huffman_codec.c")
            target_compile_definitions(pixie_pivco_avx512
                    PRIVATE PIVCO_BACKEND_AVX512=1)
            list(APPEND PIXIE_PIVCO_OBJECTS
                    $<TARGET_OBJECTS:pixie_pivco_avx512>)
        endif ()
    endif ()
elseif (PIXIE_PIVCO_PROCESSOR MATCHES "^(aarch64|arm64)$")
    list(APPEND PIXIE_PIVCO_DEFINITIONS PIVCO_HAS_NEON=1)
    add_library(pixie_pivco_neon OBJECT
            "${PIXIE_PIVCO_SRC}/pivco_huffman_codec.c")
    target_compile_definitions(pixie_pivco_neon
            PRIVATE PIVCO_BACKEND_NEON=1)
    list(APPEND PIXIE_PIVCO_OBJECTS
            $<TARGET_OBJECTS:pixie_pivco_neon>)
    list(APPEND PIXIE_PIVCO_COMMON_SOURCES
            "${PIXIE_PIVCO_SRC}/pivco_huffman_neon_tables.c")
endif ()

set(PIXIE_PIVCO_CODEC_TARGETS pixie_pivco_scalar)
if (TARGET pixie_pivco_x86)
    list(APPEND PIXIE_PIVCO_CODEC_TARGETS pixie_pivco_x86)
endif ()
if (TARGET pixie_pivco_avx512)
    list(APPEND PIXIE_PIVCO_CODEC_TARGETS pixie_pivco_avx512)
endif ()
if (TARGET pixie_pivco_neon)
    list(APPEND PIXIE_PIVCO_CODEC_TARGETS pixie_pivco_neon)
endif ()

foreach (codec_target IN LISTS PIXIE_PIVCO_CODEC_TARGETS)
    target_include_directories(${codec_target}
            PRIVATE
            "${PIXIE_PIVCO_INCLUDE}"
            "${PIXIE_PIVCO_SRC}")
    target_compile_definitions(${codec_target}
            PRIVATE ${PIXIE_PIVCO_DEFINITIONS})
    target_compile_options(${codec_target}
            PRIVATE -O3 ${PIXIE_PIVCO_OPTIONS})
endforeach ()

add_library(pixie_pivco STATIC
        ${PIXIE_PIVCO_COMMON_SOURCES}
        ${PIXIE_PIVCO_OBJECTS})
add_library(pixie::pivco ALIAS pixie_pivco)
target_include_directories(pixie_pivco
        PUBLIC
        "${PIXIE_PIVCO_INCLUDE}"
        PRIVATE
        "${PIXIE_PIVCO_SRC}")
target_compile_definitions(pixie_pivco
        PRIVATE ${PIXIE_PIVCO_DEFINITIONS})
target_compile_options(pixie_pivco
        PRIVATE -O3 ${PIXIE_PIVCO_OPTIONS})
if (UNIX AND NOT APPLE)
    target_link_libraries(pixie_pivco PUBLIC m)
endif ()
