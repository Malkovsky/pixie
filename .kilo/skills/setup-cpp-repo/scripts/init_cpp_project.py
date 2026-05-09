#!/usr/bin/env python3
"""
init_cpp_project.py - Scaffold a new C++20 repository following Pixie conventions.

Usage:
    init_cpp_project.py --name <project-name> [--namespace <namespace>] [--output-dir <dir>]

Example:
    init_cpp_project.py --name my-lib --namespace mylib --output-dir .
"""

import argparse
import os
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def to_upper(name: str) -> str:
    """Convert project name to uppercase with underscores."""
    return name.replace("-", "_").upper()


def to_snake(name: str) -> str:
    """Convert project name to snake_case for filenames."""
    return name.replace("-", "_")


# ---------------------------------------------------------------------------
# Templates
# ---------------------------------------------------------------------------

CMAKE_LISTS_TXT = """cmake_minimum_required(VERSION 3.18)
project({{PROJECT_NAME}})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(MARCH "native" CACHE STRING "march compiler flag")
add_compile_options("-march=${MARCH}")
message(STATUS "MARCH is '${MARCH}'")

option(DISABLE_AVX512 "Disable AVX512 instructions" OFF)
if(DISABLE_AVX512)
    add_compile_options("-mno-avx512f")
    message(STATUS "DISABLE_AVX512 is ON")
endif()

option(ENABLE_ADDRESS_SANITIZER "Enable AddressSanitizer" OFF)
if(ENABLE_ADDRESS_SANITIZER)
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
    message(STATUS "AddressSanitizer is ON")
endif()

option({{PROJECT_NAME_UPPER}}_COVERAGE "Enable coverage instrumentation" OFF)
if({{PROJECT_NAME_UPPER}}_COVERAGE)
    add_compile_options(-O0 -g --coverage)
    add_link_options(--coverage)
    message(STATUS "Coverage instrumentation is ON")
endif()

# ---------------------------------------------------------------------------
# Build options
# ---------------------------------------------------------------------------
option({{PROJECT_NAME_UPPER}}_TESTS "Build unit tests" ON)
option({{PROJECT_NAME_UPPER}}_BENCHMARKS "Build benchmarks" OFF)
option({{PROJECT_NAME_UPPER}}_DIAGNOSTICS "Include diagnostic logs" OFF)
option({{PROJECT_NAME_UPPER}}_DOCS "Build Doxygen documentation" OFF)

if({{PROJECT_NAME_UPPER}}_DIAGNOSTICS)
    add_compile_definitions({{PROJECT_NAME_UPPER}}_DIAGNOSTICS)
    set({{PROJECT_NAME_UPPER}}_DIAGNOSTICS_LIBS spdlog::spdlog_header_only)
endif()

# ---------------------------------------------------------------------------
# Dependencies (fetched only when needed)
# ---------------------------------------------------------------------------
include(FetchContent)

if({{PROJECT_NAME_UPPER}}_DIAGNOSTICS)
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
endif()

if({{PROJECT_NAME_UPPER}}_BENCHMARKS)
    FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.9.4
    )
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable Google Benchmark tests")
    FetchContent_MakeAvailable(googlebenchmark)
endif()

if({{PROJECT_NAME_UPPER}}_TESTS)
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
if({{PROJECT_NAME_UPPER}}_TESTS)
    enable_testing()

    add_executable(unittests
        src/tests/unittests.cpp)
    target_include_directories(unittests
        PUBLIC include)
    target_link_libraries(unittests
        gtest_main
        ${{{PROJECT_NAME_UPPER}}_DIAGNOSTICS_LIBS})
    gtest_discover_tests(unittests)
endif()

# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------
if({{PROJECT_NAME_UPPER}}_BENCHMARKS)
    add_executable(benchmarks
        src/benchmarks/benchmarks.cpp)
    target_include_directories(benchmarks
        PUBLIC include)
    target_link_libraries(benchmarks
        benchmark
        benchmark_main
        ${{{PROJECT_NAME_UPPER}}_DIAGNOSTICS_LIBS})
endif()

# ---------------------------------------------------------------------------
# Documentation (Doxygen)
# ---------------------------------------------------------------------------
if({{PROJECT_NAME_UPPER}}_DOCS)
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
"""

CMAKE_PRESETS_JSON = """{
    "version": 4,
    "cmakeMinimumRequired": {
        "major": 3,
        "minor": 18,
        "patch": 0
    },
    "configurePresets": [
        {
            "name": "base",
            "hidden": true,
            "cacheVariables": {
                "CMAKE_EXPORT_COMPILE_COMMANDS": "ON"
            }
        },
        {
            "name": "debug",
            "displayName": "Debug",
            "inherits": "base",
            "binaryDir": "${sourceDir}/build/debug",
            "cacheVariables": {
                "CMAKE_BUILD_TYPE": "Debug"
            }
        },
        {
            "name": "release",
            "displayName": "Release",
            "inherits": "base",
            "binaryDir": "${sourceDir}/build/release",
            "cacheVariables": {
                "CMAKE_BUILD_TYPE": "Release"
            }
        },
        {
            "name": "benchmarks",
            "displayName": "Benchmarks",
            "inherits": "base",
            "binaryDir": "${sourceDir}/build/benchmarks",
            "cacheVariables": {
                "CMAKE_BUILD_TYPE": "Release",
                "{{PROJECT_NAME_UPPER}}_BENCHMARKS": "ON"
            }
        },
        {
            "name": "benchmarks-diagnostic",
            "displayName": "Benchmarks diagnostic build",
            "inherits": "base",
            "binaryDir": "${sourceDir}/build/release-with-deb",
            "cacheVariables": {
                "BENCHMARK_ENABLE_LIBPFM": "ON",
                "CMAKE_BUILD_TYPE": "RelWithDebInfo",
                "{{PROJECT_NAME_UPPER}}_DIAGNOSTICS": "ON",
                "{{PROJECT_NAME_UPPER}}_BENCHMARKS": "ON"
            }
        },
        {
            "name": "docs",
            "displayName": "Docs",
            "inherits": "base",
            "binaryDir": "${sourceDir}/build/docs",
            "cacheVariables": {
                "CMAKE_BUILD_TYPE": "Release",
                "{{PROJECT_NAME_UPPER}}_DOCS": "ON"
            }
        },
        {
            "name": "coverage",
            "displayName": "Coverage",
            "inherits": "base",
            "binaryDir": "${sourceDir}/build/coverage",
            "cacheVariables": {
                "CMAKE_BUILD_TYPE": "Debug",
                "{{PROJECT_NAME_UPPER}}_BENCHMARKS": "OFF",
                "{{PROJECT_NAME_UPPER}}_COVERAGE": "ON"
            }
        },
        {
            "name": "asan",
            "displayName": "AddressSanitizer",
            "inherits": "base",
            "binaryDir": "${sourceDir}/build/asan",
            "cacheVariables": {
                "CMAKE_BUILD_TYPE": "Debug",
                "{{PROJECT_NAME_UPPER}}_BENCHMARKS": "OFF",
                "ENABLE_ADDRESS_SANITIZER": "ON"
            }
        }
    ],
    "buildPresets": [
        {
            "name": "debug",
            "displayName": "Build Debug",
            "configurePreset": "debug"
        },
        {
            "name": "release",
            "displayName": "Build Release",
            "configurePreset": "release"
        },
        {
            "name": "benchmarks",
            "displayName": "Build Benchmarks",
            "configurePreset": "benchmarks"
        },
        {
            "name": "benchmarks-diagnostic",
            "displayName": "Benchmarks diagnostic",
            "configurePreset": "benchmarks-diagnostic"
        },
        {
            "name": "docs",
            "displayName": "Build Docs",
            "configurePreset": "docs",
            "targets": [
                "docs"
            ]
        },
        {
            "name": "coverage",
            "displayName": "Build Coverage",
            "configurePreset": "coverage"
        },
        {
            "name": "asan",
            "displayName": "Build AddressSanitizer",
            "configurePreset": "asan"
        }
    ]
}
"""

CLANG_FORMAT = """# Defines the Chromium style for automatic reformatting.
# http://clang.llvm.org/docs/ClangFormatStyleOptions.html
BasedOnStyle: Chromium
# This defaults to 'Auto'. Explicitly set it for a while, so that
# 'vector<vector<int> >' in existing files gets formatted to
# 'vector<vector<int>>'. ('Auto' means that clang-format will only use
# 'int>>' if the file already contains at least one such instance.)
Standard: Cpp11

# TODO(crbug.com/1392808): Remove when InsertBraces has been upstreamed into
# the Chromium style (is implied by BasedOnStyle: Chromium).
InsertBraces: true
InsertNewlineAtEOF: true

# Sort #includes by following
# https://google.github.io/styleguide/cppguide.html#Names_and_Order_of_Includes
IncludeBlocks: Regroup
IncludeCategories:
  # C system headers.
  - Regex:           '^<.*\\.h>'
    Priority:        1
  # C++ standard library headers.
  - Regex:           '^<.*>'
    Priority:        2
  # Project headers (quoted includes).
  - Regex:           '^".*"'
    Priority:        3
  # Other libraries.
  - Regex:           '.*'
    Priority:        4
"""

GITIGNORE = """build/
.vscode/
Testing/
plans/*
venv/
docs/*
src/docs/presentations/*
CMakeUserPresets.json
_deps/
*.gcda
*.gcno
*.gcov
"""

README_MD = """# {{PROJECT_NAME}}

{{PROJECT_NAME}} is a C++20 header-only library.

## Build

```bash
cmake --preset release
cmake --build --preset release -j
./build/release/unittests
```
"""

AGENTS_MD = """# AGENTS.md - AI Coding Assistant Guidelines for {{PROJECT_NAME}}

## Project Overview

{{PROJECT_NAME}} is a **C++20 header-only library**. It provides [TODO: brief description].

## Skills

./.kilo/skills/ contains project-specific skills, use them when appropriate.

## Architecture

### Project Layout Conventions

- **`include/`**: Header-only library API (all implementations here, no `.cpp` files)
- **`src/*_tests.cpp`**: Unit tests (Google Test)
- **`src/*_benchmarks.cpp`**: Performance benchmarks (Google Benchmark)
- **`src/docs/`**: Doxygen configuration

### Key Design Decisions

1. **Header-only library**: All code in `include/`; no compiled library.
2. **Non-owning spans**: Use `std::span<const T>` for external data where appropriate.
3. **SIMD conditional compilation**: Use `#ifdef {{PROJECT_NAME_UPPER}}_AVX512_SUPPORT` / `{{PROJECT_NAME_UPPER}}_AVX2_SUPPORT` with scalar fallbacks.
4. **Target domain**: Optimized for practical data sizes.
5. **Platform**: Linux/Unix is the primary target platform.

### Why Header-Only?

- **SIMD flexibility**: Users compile with their target `-march` flags.
- **Better inlining**: Compiler sees full implementation.
- **No ABI issues**: Works across compilers and standard library versions.
- **Easy integration**: Users just `#include` headers.
- **Template-friendly**: No explicit instantiation needed.

## Technology Stack

- **Language**: C++20 (required features: `std::span`, `std::popcount`, `<bit>`)
- **Build**: CMake >= 3.18
- **Testing**: Google Test v1.17.0
- **Benchmarking**: Google Benchmark v1.9.4
- **SIMD**: AVX-512 (primary), AVX2 (fallback), scalar fallbacks
- **Style**: Chromium C++ style (`.clang-format`)

### Dependencies

The library itself is header-only and has **no runtime dependencies**. Build-time dependencies are managed via CMake FetchContent and controlled by options:

| Option | Default | What it enables |
|--------|---------|-----------------|
| `{{PROJECT_NAME_UPPER}}_TESTS` | `ON` | Unit tests (fetches Google Test) |
| `{{PROJECT_NAME_UPPER}}_BENCHMARKS` | `OFF` | Benchmarks (fetches Google Benchmark) |

## Build Commands

```bash
# Standard build (Release)
cmake -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j

# Debug build
cmake -B build/debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/debug -j

# Without AVX-512 (AVX2 fallback)
cmake -B build/release -DDISABLE_AVX512=ON
cmake --build build/release -j

# With AddressSanitizer
cmake -B build/asan -DENABLE_ADDRESS_SANITIZER=ON
cmake --build build/asan -j

# Custom march flag
cmake -B build/release -DMARCH=icelake-client
cmake --build build/release -j

# Tests only (no benchmarks)
cmake -B build/release -D{{PROJECT_NAME_UPPER}}_BENCHMARKS=OFF
cmake --build build/release -j
```

## Testing

### Running Tests

```bash
./build/release/unittests
```

### Testing Patterns

- **Differential testing**: Compare against naive reference implementations.
- **Randomized testing**: Random inputs with configurable seed.
- **Exhaustive short inputs**: Test all patterns for small sizes.

## Code Style Guidelines

1. **Formatting**: Run `clang-format` before committing (Chromium style)
2. **Namespace**: All library code in `{{NAMESPACE}}` namespace
3. **Documentation**: Use Doxygen-style comments for public API
4. **Constants**: Use `constexpr` for compile-time values
5. **Alignment**: Be aware of data alignment; prefer 64-byte aligned array allocations where performance matters

## CI/CD Workflows

- **build-test.yml**: Builds and runs tests with AddressSanitizer
- **linter.yml**: Clang-format checks on all C/C++ files
- **coverage.yml**: Coverage reporting with codecov upload
- **doxygen.yml**: Documentation generation and GitHub Pages deployment

## Common Tasks for AI Agents

### Adding a New Component

1. Create header in `include/{{NAMESPACE}}/` with Doxygen documentation
2. Add unit tests in `src/tests/<name>_tests.cpp`
3. Add benchmarks in `src/benchmarks/<name>_benchmarks.cpp`
4. Update `CMakeLists.txt` with new executables
5. Run `clang-format` on new files

### Modifying SIMD Code

1. Provide implementations for:
   - AVX-512 (`#ifdef {{PROJECT_NAME_UPPER}}_AVX512_SUPPORT`)
   - AVX2 (`#ifdef {{PROJECT_NAME_UPPER}}_AVX2_SUPPORT`)
   - Scalar fallback
2. Test with `-DDISABLE_AVX512=ON` to verify fallback works
3. Benchmark to ensure performance is maintained

### Adding Tests

1. Use Google Test framework
2. Include naive reference implementation for differential testing
3. Add edge cases: empty input, single element, boundary conditions
4. Use random testing with configurable seed for reproducibility

## Performance Philosophy

- **Goal**: Best practical performance (not just asymptotic complexity)
- **Approach**: Benchmark-driven optimization using Google Benchmark
- **SIMD**: Leverage vectorized operations where beneficial
- **Cache efficiency**: Align data structures to cache line boundaries (64 bytes)
"""

DOXYFILE_IN = """# Doxyfile

DOXYFILE_ENCODING      = UTF-8
PROJECT_NAME           = "{{PROJECT_NAME}}"
PROJECT_NUMBER         =
PROJECT_BRIEF          =
PROJECT_LOGO           =
PROJECT_ICON           =
OUTPUT_DIRECTORY       = docs
CREATE_SUBDIRS         = NO
CREATE_SUBDIRS_LEVEL   = 8
ALLOW_UNICODE_NAMES    = NO
OUTPUT_LANGUAGE        = English
BRIEF_MEMBER_DESC      = YES
REPEAT_BRIEF           = YES
ABBREVIATE_BRIEF       = "The $name class" \
                         "The $name widget" \
                         "The $name file" \
                         is \
                         provides \
                         specifies \
                         contains \
                         represents \
                         a \
                         an \
                         the
ALWAYS_DETAILED_SEC    = NO
INLINE_INHERITED_MEMB  = NO
FULL_PATH_NAMES        = YES
STRIP_FROM_PATH        = @CMAKE_CURRENT_SOURCE_DIR@
STRIP_FROM_INC_PATH    =
SHORT_NAMES            = NO
JAVADOC_AUTOBRIEF      = NO
JAVADOC_BANNER         = NO
QT_AUTOBRIEF           = NO
MULTILINE_CPP_IS_BRIEF = NO
PYTHON_DOCSTRING       = YES
INHERIT_DOCS           = YES
SEPARATE_MEMBER_PAGES  = NO
TAB_SIZE               = 4
ALIASES                =
OPTIMIZE_OUTPUT_FOR_C  = NO
OPTIMIZE_OUTPUT_JAVA   = NO
OPTIMIZE_FOR_FORTRAN   = NO
OPTIMIZE_OUTPUT_VHDL   = NO
OPTIMIZE_OUTPUT_SLICE  = NO
EXTENSION_MAPPING      =
MARKDOWN_SUPPORT       = YES
MARKDOWN_STRICT        = YES
TOC_INCLUDE_HEADINGS   = 6
MARKDOWN_ID_STYLE      = DOXYGEN
AUTOLINK_SUPPORT       = YES
AUTOLINK_IGNORE_WORDS  =
BUILTIN_STL_SUPPORT    = NO
CPP_CLI_SUPPORT        = NO
SIP_SUPPORT            = NO
IDL_PROPERTY_SUPPORT   = YES
DISTRIBUTE_GROUP_DOC   = NO
GROUP_NESTED_COMPOUNDS = NO
SUBGROUPING            = YES
INLINE_GROUPED_CLASSES = NO
INLINE_SIMPLE_STRUCTS  = NO
TYPEDEF_HIDES_STRUCT   = NO
LOOKUP_CACHE_SIZE      = 0
NUM_PROC_THREADS       = 1
TIMESTAMP              = NO
EXTRACT_ALL            = NO
EXTRACT_PRIVATE        = NO
EXTRACT_PRIV_VIRTUAL   = NO
EXTRACT_PACKAGE        = NO
EXTRACT_STATIC         = NO
EXTRACT_LOCAL_CLASSES  = YES
EXTRACT_LOCAL_METHODS  = NO
EXTRACT_ANON_NSPACES   = NO
RESOLVE_UNNAMED_PARAMS = YES
HIDE_UNDOC_MEMBERS     = NO
HIDE_UNDOC_CLASSES     = NO
HIDE_UNDOC_NAMESPACES  = YES
HIDE_FRIEND_COMPOUNDS  = NO
HIDE_IN_BODY_DOCS      = NO
INTERNAL_DOCS          = NO
CASE_SENSE_NAMES       = SYSTEM
HIDE_SCOPE_NAMES       = NO
HIDE_COMPOUND_REFERENCE= NO
SHOW_HEADERFILE        = YES
SHOW_INCLUDE_FILES     = YES
SHOW_GROUPED_MEMB_INC  = NO
FORCE_LOCAL_INCLUDES   = NO
INLINE_INFO            = YES
SORT_MEMBER_DOCS       = YES
SORT_BRIEF_DOCS        = NO
SORT_MEMBERS_CTORS_1ST = NO
SORT_GROUP_NAMES       = NO
SORT_BY_SCOPE_NAME     = NO
STRICT_PROTO_MATCHING  = NO
GENERATE_TODOLIST      = YES
GENERATE_TESTLIST      = YES
GENERATE_BUGLIST       = YES
GENERATE_DEPRECATEDLIST= YES
ENABLED_SECTIONS       =
MAX_INITIALIZER_LINES  = 30
SHOW_USED_FILES        = YES
SHOW_FILES             = YES
SHOW_NAMESPACES        = YES
FILE_VERSION_FILTER    =
LAYOUT_FILE            =
CITE_BIB_FILES         =
EXTERNAL_TOOL_PATH     =
QUIET                  = NO
WARNINGS               = YES
WARN_IF_UNDOCUMENTED   = YES
WARN_IF_DOC_ERROR      = YES
WARN_IF_INCOMPLETE_DOC = YES
WARN_NO_PARAMDOC       = NO
WARN_IF_UNDOC_ENUM_VAL = NO
WARN_LAYOUT_FILE       = YES
WARN_AS_ERROR          = NO
WARN_FORMAT            = "$file:$line: $text"
WARN_LINE_FORMAT       = "at line $line of file $file"
WARN_LOGFILE           =
INPUT                  = @CMAKE_CURRENT_SOURCE_DIR@/include \
                         @CMAKE_CURRENT_SOURCE_DIR@/README.md
INPUT_ENCODING         = UTF-8
INPUT_FILE_ENCODING    =
FILE_PATTERNS          = *.c \
                         *.cc \
                         *.cxx \
                         *.cpp \
                         *.h \
                         *.hh \
                         *.hxx \
                         *.hpp
RECURSIVE              = YES
EXCLUDE                =
EXCLUDE_SYMLINKS       = NO
EXCLUDE_PATTERNS       =
EXCLUDE_SYMBOLS        =
EXAMPLE_PATH           =
EXAMPLE_PATTERNS       = *
EXAMPLE_RECURSIVE      = NO
IMAGE_PATH             = @CMAKE_CURRENT_SOURCE_DIR@/src/docs/images
INPUT_FILTER           =
FILTER_PATTERNS        =
FILTER_SOURCE_FILES    = NO
FILTER_SOURCE_PATTERNS =
USE_MDFILE_AS_MAINPAGE = @CMAKE_CURRENT_SOURCE_DIR@/README.md
IMPLICIT_DIR_DOCS      = YES
FORTRAN_COMMENT_AFTER  = 72
SOURCE_BROWSER         = NO
INLINE_SOURCES         = NO
STRIP_CODE_COMMENTS    = YES
REFERENCED_BY_RELATION = NO
REFERENCES_RELATION    = NO
REFERENCES_LINK_SOURCE = YES
SOURCE_TOOLTIPS        = YES
USE_HTAGS              = NO
VERBATIM_HEADERS       = YES
CLANG_ASSISTED_PARSING = NO
CLANG_ADD_INC_PATHS    = YES
CLANG_OPTIONS          =
CLANG_DATABASE_PATH    =
ALPHABETICAL_INDEX     = YES
IGNORE_PREFIX          =
GENERATE_HTML          = YES
HTML_OUTPUT            = html
HTML_FILE_EXTENSION    = .html
HTML_HEADER            =
HTML_FOOTER            =
HTML_STYLESHEET        =
HTML_EXTRA_STYLESHEET  = @AWESOME_CSS_DIR@/doxygen-awesome.css
HTML_EXTRA_FILES       =
HTML_COLORSTYLE        = AUTO_LIGHT
HTML_COLORSTYLE_HUE    = 220
HTML_COLORSTYLE_SAT    = 100
HTML_COLORSTYLE_GAMMA  = 80
HTML_DYNAMIC_MENUS     = YES
HTML_DYNAMIC_SECTIONS  = NO
HTML_CODE_FOLDING      = YES
HTML_COPY_CLIPBOARD    = YES
HTML_PROJECT_COOKIE    =
HTML_INDEX_NUM_ENTRIES = 100
GENERATE_DOCSET        = NO
DOCSET_FEEDNAME        = "Doxygen generated docs"
DOCSET_FEEDURL         =
DOCSET_BUNDLE_ID       = org.doxygen.Project
DOCSET_PUBLISHER_ID    = org.doxygen.Publisher
DOCSET_PUBLISHER_NAME  = Publisher
GENERATE_HTMLHELP      = NO
CHM_FILE               =
HHC_LOCATION           =
GENERATE_CHI           = NO
CHM_INDEX_ENCODING     =
BINARY_TOC             = NO
TOC_EXPAND             = NO
SITEMAP_URL            =
GENERATE_QHP           = NO
QCH_FILE               =
QHP_NAMESPACE          = org.doxygen.Project
QHP_VIRTUAL_FOLDER     = doc
QHP_CUST_FILTER_NAME   =
QHP_CUST_FILTER_ATTRS  =
QHP_SECT_FILTER_ATTRS  =
QHG_LOCATION           =
GENERATE_ECLIPSEHELP   = NO
ECLIPSE_DOC_ID         = org.doxygen.Project
DISABLE_INDEX          = NO
GENERATE_TREEVIEW      = YES
PAGE_OUTLINE_PANEL     = YES
FULL_SIDEBAR           = NO
ENUM_VALUES_PER_LINE   = 4
SHOW_ENUM_VALUES       = NO
TREEVIEW_WIDTH         = 250
EXT_LINKS_IN_WINDOW    = NO
OBFUSCATE_EMAILS       = YES
HTML_FORMULA_FORMAT    = png
FORMULA_FONTSIZE       = 10
FORMULA_MACROFILE      =
USE_MATHJAX            = NO
MATHJAX_VERSION        = MathJax_2
MATHJAX_FORMAT         = HTML-CSS
MATHJAX_RELPATH        =
MATHJAX_EXTENSIONS     =
MATHJAX_CODEFILE       =
SEARCHENGINE           = YES
SERVER_BASED_SEARCH    = NO
EXTERNAL_SEARCH        = NO
SEARCHENGINE_URL       =
SEARCHDATA_FILE        = searchdata.xml
EXTERNAL_SEARCH_ID     =
EXTRA_SEARCH_MAPPINGS  =
GENERATE_LATEX         = NO
"""

COVERAGE_REPORT_SH = """#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build/coverage"

cmake --preset coverage
cmake --build --preset coverage

"${BUILD_DIR}/unittests"

cd "${BUILD_DIR}"
find . -name "*.gcda" > gcov_files.txt
while read -r f; do
  case "${f}" in
    *"/_deps/"*|*"/third_party/"*|*"/src/benchmarks/"*)
      ;;
    *)
      gcov -pb "${f}" >> coverage.txt
      ;;
  esac
done < gcov_files.txt
echo "gcov report written to ${BUILD_DIR}/coverage.txt"
"""

BUILD_TEST_YML = """name: Tests (ASan)

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  build-and-test:
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v4

    - name: Create Build Directory
      run: mkdir build

    - name: Configure CMake
      working-directory: ./build
      run: cmake -DDISABLE_AVX512=ON -DENABLE_ADDRESS_SANITIZER=ON -D{{PROJECT_NAME_UPPER}}_BENCHMARKS=OFF ..

    - name: Build Project
      working-directory: ./build
      run: make -j

    - name: Run Unittests
      working-directory: ./build
      run: ./unittests
"""

LINTER_YML = """name: Clang Format Lint

on:
  pull_request:
  push:
    branches: [main]

jobs:
  clang-format:
    runs-on: ubuntu-latest

    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Install clang-format
        run: sudo apt-get update && sudo apt-get install -y clang-format

      - name: Run clang-format check
        run: |
          mapfile -t FILES < <(find include src -type f \\( -name '*.cpp' -o -name '*.hpp' -o -name '*.cc' -o -name '*.c' -o -name '*.h' \\))
          clang-format --version
          if [ ${#FILES[@]} -eq 0 ]; then
            echo "No C/C++ files found."
            exit 0
          fi

          clang-format --dry-run --Werror "${FILES[@]}"
"""

COVERAGE_YML = """name: coverage

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  coverage:
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v4

    - name: Create Build Directory
      run: mkdir build

    - name: Run coverage
      run: ./scripts/coverage_report.sh

    - name: Upload to Codecov
      uses: codecov/codecov-action@v4
      with:
        token: ${{ secrets.CODECOV_TOKEN }}
        files: build/coverage/coverage.txt
        flags: gcov
        fail_ci_if_error: false

    - name: Upload coverage artifacts
      uses: actions/upload-artifact@v4
      with:
        name: coverage-gcov
        path: |
          build/coverage/coverage.txt
          build/coverage/*.gcov
"""

DOXYGEN_YML = """# Simple workflow for deploying static content to GitHub Pages
name: Deploy static content to Pages

on:
  # Runs on pushes targeting the default branch
  push:
    branches: ["main"]

  # Allows you to run this workflow manually from the Actions tab
  workflow_dispatch:

# Sets permissions of the GITHUB_TOKEN to allow deployment to GitHub Pages
permissions:
  contents: read
  pages: write
  id-token: write

# Allow only one concurrent deployment, skipping runs queued between the run in-progress and latest queued.
# However, do NOT cancel in-progress runs as we want to allow these production deployments to complete.
concurrency:
  group: "pages"
  cancel-in-progress: false

jobs:
  # Single deploy job since we're just deploying
  deploy:
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    runs-on: ubuntu-latest
    steps:
      - name: Checkout
        uses: actions/checkout@v4
      - name: Install Doxygen v1.13.2
        run: |
          transformed_version=$(echo "1.13.2" | tr '.' '_')
          wget https://github.com/doxygen/doxygen/releases/download/Release_${transformed_version}/doxygen-1.13.2.linux.bin.tar.gz
          tar -xzf doxygen-1.13.2.linux.bin.tar.gz
          sudo mv doxygen-1.13.2/bin/doxygen /usr/local/bin/doxygen
        shell: bash
      - name: Cmake configure
        run: cmake -S ${{github.workspace}} -B ${{github.workspace}}/build -D{{PROJECT_NAME_UPPER}}_DOCS=ON -D{{PROJECT_NAME_UPPER}}_TESTS=OFF -D{{PROJECT_NAME_UPPER}}_BENCHMARKS=OFF
      - name: Build docs
        run: cmake --build ${{github.workspace}}/build --target docs
      - name: Upload artifact
        uses: actions/upload-pages-artifact@v3
        with:
          # Upload entire repository
          path: ${{github.workspace}}/build/docs/html
      - name: Deploy to GitHub Pages
        id: deployment
        uses: actions/deploy-pages@v4
"""

HEADER_HPP = """#pragma once

/**
 * @file {{HEADER_NAME}}
 * @brief Main header for the {{PROJECT_NAME}} library
 */

namespace {{NAMESPACE}} {

/**
 * @brief Example function.
 *
 * TODO: Replace with actual library functionality.
 */
inline int example() {
    return 42;
}

}  // namespace {{NAMESPACE}}
"""

UNITTESTS_CPP = """#include <gtest/gtest.h>

#include "{{NAMESPACE}}/{{HEADER_NAME}}"

TEST(ExampleTest, BasicAssertion) {
    EXPECT_EQ({{NAMESPACE}}::example(), 42);
}
"""

BENCHMARKS_CPP = """#include <benchmark/benchmark.h>

#include "{{NAMESPACE}}/{{HEADER_NAME}}"

static void BM_Example(benchmark::State& state) {
    for (auto _ : state) {
        benchmark::DoNotOptimize({{NAMESPACE}}::example());
    }
}

BENCHMARK(BM_Example);

BENCHMARK_MAIN();
"""


# ---------------------------------------------------------------------------
# Generation logic
# ---------------------------------------------------------------------------

def generate(args: argparse.Namespace) -> None:
    project_name = args.name
    namespace = args.namespace or project_name.replace("-", "")
    project_name_upper = to_upper(project_name)
    header_name = f"{to_snake(project_name)}.hpp"
    output_dir = Path(args.output_dir).resolve() / project_name

    if output_dir.exists():
        print(f"Error: output directory already exists: {output_dir}")
        sys.exit(1)

    substitutions = {
        "{{PROJECT_NAME}}": project_name,
        "{{NAMESPACE}}": namespace,
        "{{PROJECT_NAME_UPPER}}": project_name_upper,
        "{{HEADER_NAME}}": header_name,
    }

    def sub(text: str) -> str:
        for key, value in substitutions.items():
            text = text.replace(key, value)
        return text

    # Create directories
    (output_dir / "include" / namespace).mkdir(parents=True)
    (output_dir / "src" / "tests").mkdir(parents=True)
    (output_dir / "src" / "benchmarks").mkdir(parents=True)
    (output_dir / "src" / "docs").mkdir(parents=True)
    (output_dir / "src" / "docs" / "images").mkdir(parents=True)
    (output_dir / "scripts").mkdir(parents=True)
    (output_dir / ".github" / "workflows").mkdir(parents=True)

    # Write files
    files = {
        output_dir / "CMakeLists.txt": sub(CMAKE_LISTS_TXT),
        output_dir / "CMakePresets.json": sub(CMAKE_PRESETS_JSON),
        output_dir / ".clang-format": sub(CLANG_FORMAT),
        output_dir / ".gitignore": sub(GITIGNORE),
        output_dir / "README.md": sub(README_MD),
        output_dir / "AGENTS.md": sub(AGENTS_MD),
        output_dir / "src" / "docs" / "Doxyfile.in": sub(DOXYFILE_IN),
        output_dir / "scripts" / "coverage_report.sh": sub(COVERAGE_REPORT_SH),
        output_dir / ".github" / "workflows" / "build-test.yml": sub(BUILD_TEST_YML),
        output_dir / ".github" / "workflows" / "linter.yml": sub(LINTER_YML),
        output_dir / ".github" / "workflows" / "coverage.yml": sub(COVERAGE_YML),
        output_dir / ".github" / "workflows" / "doxygen.yml": sub(DOXYGEN_YML),
        output_dir / "include" / namespace / header_name: sub(HEADER_HPP),
        output_dir / "src" / "tests" / "unittests.cpp": sub(UNITTESTS_CPP),
        output_dir / "src" / "benchmarks" / "benchmarks.cpp": sub(BENCHMARKS_CPP),
    }

    for path, content in files.items():
        path.write_text(content)
        print(f"Created: {path.relative_to(output_dir.parent)}")

    # Make coverage script executable
    (output_dir / "scripts" / "coverage_report.sh").chmod(0o755)

    print(f"\\nProject '{project_name}' generated successfully at {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scaffold a new C++20 repository following Pixie conventions."
    )
    parser.add_argument("--name", required=True, help="Project name (e.g., my-lib)")
    parser.add_argument(
        "--namespace",
        help="C++ namespace (defaults to project name with hyphens removed)",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Output directory (default: current directory)",
    )
    args = parser.parse_args()
    generate(args)


if __name__ == "__main__":
    main()
