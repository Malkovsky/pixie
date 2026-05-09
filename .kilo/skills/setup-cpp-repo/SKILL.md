---
name: setup-cpp-repo
description: Scaffold a new C++20 repository with CMake, Google Test, Google Benchmark, CI workflows, Doxygen docs, and Chromium code style. Use when the user asks to create a new C++ project, set up a C++ library, or initialize a C++ repository with modern tooling.
---

# setup-cpp-repo

## Overview

This skill generates a complete C++20 project scaffold following the conventions of the Pixie succinct data structures library. The generated repository is header-only by default and includes:

- CMake build system with presets
- Google Test for unit testing
- Google Benchmark for performance benchmarks
- Doxygen documentation with doxygen-awesome-css theme
- GitHub Actions CI workflows (ASan, lint, coverage, docs)
- Chromium C++ code style via `.clang-format`
- `AGENTS.md` for AI coding assistant guidelines

## When to Use This Skill

Use this skill when:
- The user wants to create a new C++ library or project from scratch
- The user asks for a "C++ project template" or "C++ repo setup"
- The user needs CMake + Google Test + benchmark scaffolding
- The user wants to follow Pixie-style conventions (header-only, AVX-512 optional, Doxygen docs)

Do **not** use this skill when:
- Working with an existing codebase (use the `cmake` skill instead)
- The project is not C++ (use a different skill)
- The user only wants a single file or snippet

## Workflow

### Step 1: Gather Parameters

Ask the user for (or infer from context):
- **Project name** (required): Hyphenated lowercase identifier, e.g., `my-lib`
- **Namespace** (optional): C++ namespace. Defaults to project name with hyphens removed, e.g., `mylib`
- **Output directory** (optional): Where to create the project. Defaults to current directory.

### Step 2: Run the Generator

Execute the generation script:

```bash
python3 .kilo/skills/setup-cpp-repo/scripts/init_cpp_project.py \
    --name <project-name> \
    [--namespace <namespace>] \
    [--output-dir <directory>]
```

Example:
```bash
python3 .kilo/skills/setup-cpp-repo/scripts/init_cpp_project.py \
    --name succinct-lib --namespace succinct --output-dir .
```

### Step 3: Verify the Scaffold

After generation, the project structure should look like:

```
<project-name>/
├── CMakeLists.txt
├── CMakePresets.json
├── .clang-format
├── .gitignore
├── README.md
├── AGENTS.md
├── include/
│   └── <namespace>/
│       └── <project_snake>.hpp
├── src/
│   ├── tests/
│   │   └── unittests.cpp
│   ├── benchmarks/
│   │   └── benchmarks.cpp
│   └── docs/
│       ├── Doxyfile.in
│       └── images/
├── scripts/
│   └── coverage_report.sh
└── .github/
    └── workflows/
        ├── build-test.yml
        ├── linter.yml
        ├── coverage.yml
        └── doxygen.yml
```

### Step 4: Initial Build and Test

Change into the project directory and run an initial build to verify everything works:

```bash
cd <project-name>
cmake --preset release
cmake --build --preset release -j
./build/release/unittests
```

If the build and tests pass, the scaffold is ready.

### Step 5: Hand Off to cmake Skill

After project creation, use the **`cmake` skill** (`.kilo/skills/cmake/SKILL.md`) for all subsequent build operations. The `cmake` skill documents:
- Build directory conventions with git short-hash suffixes
- How to replicate preset settings with custom build directories
- AddressSanitizer, coverage, and benchmark workflows
- Best practices for out-of-source builds

## Customization Guide

### Adding More Test Executables

Edit `CMakeLists.txt` and add new `add_executable` blocks under the `if(<PROJECT_UPPER>_TESTS)` section, following the pattern of the existing `unittests` target.

Update `scripts/coverage_report.sh` to run any new test binaries.

Update `.github/workflows/build-test.yml` to execute new test binaries in CI.

### Adding More Benchmark Executables

Edit `CMakeLists.txt` and add new `add_executable` blocks under the `if(<PROJECT_UPPER>_BENCHMARKS)` section, following the pattern of the existing `benchmarks` target.

### Adding Third-Party Dependencies

For header-only libraries, prefer `FetchContent` in `CMakeLists.txt`. For compiled libraries, consider vendoring or using a package manager (Conan, vcpkg).

### Modifying Doxygen Configuration

Edit `src/docs/Doxyfile.in`. The generated version is intentionally minimal (only non-default settings). Add or override settings as needed. Run `doxygen -g` to see all available options.

## Reference

See `references/project_structure.md` for a detailed breakdown of every generated file and its purpose.
