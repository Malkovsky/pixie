#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build/coverage"
JOBS="${PIXIE_COVERAGE_JOBS:-$(nproc)}"

cmake --preset coverage
cmake --build --preset coverage --parallel "${JOBS}"

find "${BUILD_DIR}" -name "*.gcda" -delete
find "${BUILD_DIR}" -name "*.gcov" -delete
rm -f "${BUILD_DIR}/coverage.txt" "${BUILD_DIR}/gcov_files.txt"

ctest --preset coverage --parallel "${JOBS}"

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
