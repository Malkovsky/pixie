cmake_minimum_required(VERSION 3.18)

if (NOT DEFINED ROOT_SOURCE_DIR)
    message(FATAL_ERROR "ROOT_SOURCE_DIR is required")
endif ()

if (NOT DEFINED SDSL_SOURCE_DIR)
    message(FATAL_ERROR "SDSL_SOURCE_DIR is required")
endif ()

file(COPY "${ROOT_SOURCE_DIR}/third_party/sdsl-lite/CMakeLists.txt"
     DESTINATION "${SDSL_SOURCE_DIR}")

set(memory_header "${SDSL_SOURCE_DIR}/include/sdsl/memory_management.hpp")

if (EXISTS "${SDSL_SOURCE_DIR}/.git")
    execute_process(
        COMMAND git checkout -- include/sdsl/memory_management.hpp
        WORKING_DIRECTORY "${SDSL_SOURCE_DIR}"
        RESULT_VARIABLE git_restore_result
    )
    if (NOT git_restore_result EQUAL 0)
        message(WARNING "Failed to restore memory_management.hpp from git")
    endif ()
endif ()

file(READ "${memory_header}" memory_content)

if (memory_content MATCHES "winsock2.h")
    set(needs_windows_includes OFF)
else ()
    set(needs_windows_includes ON)
endif ()

if (needs_windows_includes)
    string(REPLACE "#include <unistd.h>"
                   "#ifdef _WIN32\n#include <winsock2.h>\n#include <windows.h>\n#include <io.h>\n#else\n#include <unistd.h>\n#endif"
                   memory_content
                   "${memory_content}")
endif ()

if (memory_content MATCHES "#ifndef _WIN32\\n#include <sys/mman.h>")
    set(needs_mman_guard OFF)
else ()
    set(needs_mman_guard ON)
endif ()

if (needs_mman_guard)
    string(REPLACE "#include <sys/mman.h> // IWYU pragma: keep"
                   "#ifndef _WIN32\n#include <sys/mman.h> // IWYU pragma: keep\n#endif"
                   memory_content
                   "${memory_content}")
endif ()

file(WRITE "${memory_header}" "${memory_content}")
