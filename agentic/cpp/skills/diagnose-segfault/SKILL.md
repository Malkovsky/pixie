---
name: diagnose-segfault
description: Diagnose C++ crashes and memory-safety errors with AddressSanitizer, GDB, and core dumps. Use when a C++ binary crashes with SIGSEGV, SIGABRT, heap-buffer-overflow, use-after-free, stack-buffer-overflow, double-free, suspected memory corruption, or an available core file.
---

# C++ Segfault and Memory Error Diagnosis

Use this skill to find the first bad access or corrupting operation, not just
the frame where the process finally crashed.

## When To Use

- A C++ binary crashes with `Segmentation fault`, `SIGSEGV`, or `SIGABRT`.
- AddressSanitizer reports `ERROR: AddressSanitizer:`.
- A test reports heap-buffer-overflow, stack-buffer-overflow, use-after-free,
  double-free, global-buffer-overflow, or similar memory-safety failures.
- A core dump exists and the user wants root-cause analysis.
- Memory corruption is suspected but the immediate failure is ambiguous.

For repository-specific binary names, CMake options, or known reproducer
patterns, also read `agentic/local/cpp/skills/diagnose-segfault/EXAMPLES.md`
when present.

## Workflow 1: ASan First

Prefer AddressSanitizer when the issue is reproducible. It usually reports the
bad access with file and line information.

### Build With ASan

For CMake projects, first check whether the repository already has an ASan
preset or option. If not, configure a dedicated debug build:

```bash
cmake -B build/asan -DCMAKE_BUILD_TYPE=Debug -DENABLE_ADDRESS_SANITIZER=ON
cmake --build build/asan -j
```

For non-CMake builds, compile and link with:

```bash
-fsanitize=address -fno-omit-frame-pointer -g -O1
```

Use `-O0` instead of `-O1` when you expect to inspect many local variables in
GDB.

### Run The Minimal Reproducer

Run the specific binary, test case, or input that triggers the crash. For Google
Test binaries, prefer a narrow filter:

```bash
./build/asan/unittests --gtest_filter="SuiteName.TestName"
```

### Read The ASan Report

Focus on:

| Report section | Meaning |
|---|---|
| `ERROR: AddressSanitizer: <type>` | Error class |
| `READ/WRITE of size N` | Access direction and size |
| First user-code frame | Exact bad access |
| Allocation/deallocation stack | Object lifetime and ownership |
| Shadow-byte legend | Boundary or lifetime category |
| `SUMMARY:` | One-line location summary |

Useful options:

```bash
ASAN_OPTIONS=detect_leaks=0:detect_stack_use_after_return=1
ASAN_OPTIONS=halt_on_error=0
ASAN_OPTIONS=print_stats=1
```

Disable leak detection while diagnosing a crash if leak noise hides the primary
failure:

```bash
ASAN_OPTIONS=detect_leaks=0 ./build/asan/unittests
```

## Workflow 2: GDB Live Debugging

Use GDB when ASan is unavailable, when the crash is not a direct memory-safety
violation, or when variable inspection is needed.

Build with debug symbols:

```bash
cmake -B build/debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/debug -j
```

Run under GDB:

```bash
gdb --args <binary> [arguments...]
```

Core commands:

```gdb
run
bt full
info registers
info locals
info args
frame N
list
print <expr>
thread apply all bt
```

Make C++ values easier to inspect:

```gdb
set print pretty on
set print object on
set pagination off
```

## Workflow 3: ASan Under GDB

Use this when ASan points at a bad access but the pointer or lifetime corruption
comes from an earlier frame.

```bash
gdb --args <asan-binary> [arguments...]
```

Break on ASan reporting or abort:

```gdb
break __asan::ReportGenericError
catch signal SIGABRT
run
bt full
```

Then inspect the last user-code frames before ASan internals.

## Workflow 4: Core Dump Analysis

Use when the crash already happened or reproduction is expensive.

Enable core dumps for future runs if needed:

```bash
ulimit -c unlimited
```

Analyze:

```bash
gdb <binary> <core-file>
```

Useful commands:

```gdb
bt full
info threads
thread apply all bt
frame N
info locals
info args
print <expr>
```

## Common ASan Errors

| Error type | Typical cause |
|---|---|
| `heap-buffer-overflow` | Read or write past heap allocation bounds |
| `stack-buffer-overflow` | Read or write past a local stack object |
| `global-buffer-overflow` | Read or write past global/static storage |
| `heap-use-after-free` | Access after `delete`, `free`, or container invalidation |
| `stack-use-after-return` | Pointer/reference to a returned stack frame |
| `double-free` | Object released twice |
| `alloc-dealloc-mismatch` | Mixed allocation APIs, such as `new[]` with `free` |

## Best Practices

1. Build with `-g`; reports without symbols are often not actionable.
2. Prefer the smallest reproducer over full-suite runs.
3. Rebuild after toggling sanitizer or debug options.
4. Treat the first ASan error as primary; later errors are often fallout.
5. Check container iterator/reference invalidation around the reported object.
6. Validate the fix with the same reproducer under ASan before running broader
   tests.
7. If ASan is too slow for a large input, use GDB on the same input or reduce
   the input while preserving the crash.
