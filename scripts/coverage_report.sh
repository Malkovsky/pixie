#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build/coverage"
JOBS="${PIXIE_COVERAGE_JOBS:-$(nproc)}"
VARIANT="${PIXIE_COVERAGE_VARIANT:-native}"
CONFIGURE_ARGS=(--preset coverage)
BUILD_ARGS=(--preset coverage)
TEST_ARGS=(--preset coverage)

case "${VARIANT}" in
  native) ;;
  rank-select-scalar)
    BUILD_DIR="${ROOT_DIR}/build/coverage-scalar"
    CONFIGURE_ARGS+=(-B "${BUILD_DIR}" -DMARCH=x86-64 -DDISABLE_AVX512=ON)
    BUILD_ARGS=("${BUILD_DIR}" --target rank_select_tests)
    TEST_ARGS=(--test-dir "${BUILD_DIR}" -C Debug -L '^rank_select_tests$'
               --no-tests=error --output-on-failure)
    ;;
  *)
    echo "Unknown coverage variant: ${VARIANT}" >&2
    exit 2
    ;;
esac

cmake "${CONFIGURE_ARGS[@]}"
cmake --build "${BUILD_ARGS[@]}" --parallel "${JOBS}"

find "${BUILD_DIR}" -name "*.gcda" -delete
find "${BUILD_DIR}" -name "*.gcov" -delete
rm -f "${BUILD_DIR}/coverage.txt" "${BUILD_DIR}/gcov_files.txt"

ctest "${TEST_ARGS[@]}" --parallel "${JOBS}"

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
