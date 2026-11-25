#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"
BUILD_DIR="${ROOT_DIR}/build/coverage"
JOBS="${PIXIE_COVERAGE_JOBS:-1}"
VARIANT="${PIXIE_COVERAGE_VARIANT:-native}"
CONFIGURE_ARGS=(--preset coverage)
BUILD_ARGS=(--preset coverage)
TEST_ARGS=(--preset coverage)

case "${VARIANT}" in
  native) ;;
  sse2)
    BUILD_DIR="${ROOT_DIR}/build/coverage-sse2"
    CONFIGURE_ARGS+=(-B "${BUILD_DIR}" -DPIXIE_CPU_PROFILE=sse2)
    SSE2_SUITES=(bit_algorithms_unittests rank_select_tests test_rmm rmq_tests
                 wavelet_tree_tests select512_experimental_tests
                 excess_positions_tests excess_record_lows_tests)
    SSE2_PATTERN="$(IFS='|'; printf '%s' "${SSE2_SUITES[*]}")"
    BUILD_ARGS=("${BUILD_DIR}" --target "${SSE2_SUITES[@]}" tree_tests)
    TEST_ARGS=(--test-dir "${BUILD_DIR}" -C Debug
               -LE '^heavy$'
               -R "^((${SSE2_PATTERN})\.|tree_tests\.(Forest256|OrdinalForest256|TreeSimdMicro|TreeUnpackingTest\.))"
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

bash "${ROOT_DIR}/scripts/collect_gcov.sh" "${BUILD_DIR}"
