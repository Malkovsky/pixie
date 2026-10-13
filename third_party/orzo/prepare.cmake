cmake_minimum_required(VERSION 3.18)

if (NOT DEFINED ORZO_SOURCE_DIR OR NOT DEFINED ORZO_OUTPUT_DIR)
    message(FATAL_ERROR "ORZO_SOURCE_DIR and ORZO_OUTPUT_DIR are required")
endif ()

# Patch a build-local copy, including when FetchContent uses a source override.
# Upstream remains intact and its MIT notice accompanies the prepared headers.
file(MAKE_DIRECTORY "${ORZO_OUTPUT_DIR}/include")
file(COPY "${ORZO_SOURCE_DIR}/include/orzo"
        DESTINATION "${ORZO_OUTPUT_DIR}/include")
file(COPY "${ORZO_SOURCE_DIR}/LICENSE" DESTINATION "${ORZO_OUTPUT_DIR}")

set(header "${ORZO_OUTPUT_DIR}/include/orzo/orzo.h")
file(READ "${header}" content)
foreach(block UPPER_BLOCK_COUNT LOWER_BLOCK_COUNT BASIC_BLOCK_COUNT SELECT_UPPER_BLOCK_COUNT)
    set(old "(size_t) std::ceil((float) bv_count / (float) ${block})")
    string(FIND "${content}" "${old}" found)
    if (found EQUAL -1)
        message(FATAL_ERROR "Orzo's allocation layout changed: ${block}")
    endif ()
    # Float rounding can underestimate directory lengths on large sources.
    string(REPLACE "${old}" "(bv_count / ${block}) + (bv_count % ${block} != 0)"
            content "${content}")
endforeach()
# Construction always writes a terminal L1/L2 entry, including exact multiples.
string(REPLACE "new __uint128_t[num_lower_blocks]()"
        "new __uint128_t[num_lower_blocks + 1]()" content "${content}")
# At an exact upper boundary the loop advances past the completed bucket before
# the final write, so both upper directories also need a terminal slot.
string(REPLACE "new uint64_t[l0_count + 1]()"
        "new uint64_t[l0_count + 2]()" content "${content}")
string(REPLACE "new uint64_t[select_l0_count + 1]()"
        "new uint64_t[select_l0_count + 2]()" content "${content}")
file(WRITE "${header}" "${content}")

set(header "${ORZO_OUTPUT_DIR}/include/orzo/utils.h")
file(READ "${header}" content)
# Avoid collision with namespace benchmark and multiple-definition linker errors.
string(REPLACE "double benchmark(" "inline double orzo_benchmark_elapsed("
        content "${content}")
file(WRITE "${header}" "${content}")
