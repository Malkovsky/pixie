#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build/coverage"

cmake --preset coverage
cmake --build --preset coverage

"${BUILD_DIR}/unittests"
"${BUILD_DIR}/excess_positions_tests"
"${BUILD_DIR}/louds_tree_tests"
"${BUILD_DIR}/test_rmm"

cd "${BUILD_DIR}"
find . -name "*.gcda" > gcov_files.txt
while read -r f; do
  case "${f}" in
    *"/third_party/"*|*"/src/benchmarks/"*)
      ;;
    *)
      gcov -pb "${f}" >> coverage.txt
      ;;
  esac
done < gcov_files.txt
echo "gcov report written to ${BUILD_DIR}/coverage.txt"
