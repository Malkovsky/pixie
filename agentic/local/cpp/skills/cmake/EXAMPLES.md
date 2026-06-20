# Pixie CMake Examples

Use these notes with `agentic/cpp/skills/cmake/SKILL.md` for Pixie's local
presets and coverage workflow.

## Coverage

Prefer the repository script when checking coverage:

```bash
./scripts/coverage_report.sh
```

The script uses the `coverage` preset, where Pixie's coverage option is
`PIXIE_COVERAGE=ON`, not the generic `ENABLE_COVERAGE=ON` spelling used in
some other CMake projects. It also clears old `.gcda` and `.gcov` files before
running tests, which avoids gcov checksum warnings after recompilation.

The generated report is:

```text
build/coverage/coverage.txt
```

Useful focused targets while iterating on coverage:

```bash
cmake --build build/coverage --target rmq_tests -j
./build/coverage/rmq_tests

cmake --build build/coverage --target test_rmm -j
./build/coverage/test_rmm
```

Run the full script before reporting final coverage numbers, because individual
test-binary runs can leave stale profile counters when source checksums change.

## Codecov Interpretation

Pixie uploads `gcov -pb` output to Codecov. The `-b` flag includes branch
probabilities, so Codecov reports "partial" coverage for executed lines whose
branch paths were not all covered. A local line-only summary can therefore look
better than the Codecov project or patch percentage.

When a PR coverage check looks unexpectedly low:

1. Check the Codecov file table for missing and partial counts.
2. Check `scripts/coverage_report.sh` and `.github/workflows/coverage.yml` to
   confirm the uploaded file and flags.
3. For header-only templates, remember that the same header may appear in more
   than one translation-unit gcov report. Inspect the generated `.gcov` file or
   Codecov aggregate before drawing conclusions from the first `coverage.txt`
   block.

## RMQ/RmM Coverage Strategy

Random differential tests are necessary but not sufficient for Pixie's RMQ and
RmM structures. The best coverage improvements usually come from deterministic
input shapes that force a specific traversal branch:

- HybridBTree leaf selector paths: prefix answer, suffix answer, embedded leaf
  minimum answer, and unsupported interior range fallback.
- HybridBTree and CartesianHybridBTree internal query paths where the local
  selector chooses a border child whose subtree minimum is outside the queried
  range, forcing border recursion plus middle-child comparison.
- Top sparse-overlay paths where the padded sparse candidate is inside the
  query and where it misses and delegates borders to the tree.
- RmMBTree partial last-block scans, forward/backward search no-hit cases,
  backward boundary-only returns, and `minselect` ranks across node boundaries.

Do not chase every uncovered branch. Low-value coverage includes allocator
failure branches, impossible index-size `length_error` paths, and private
defensive `npos` checks that cannot be reached through valid public API calls.
