#!/usr/bin/env bash
# Collect existing profiles without rebuilding or rerunning tests.
set -euo pipefail

BUILD_DIR="$(cd "${1:?Usage: collect_gcov.sh BUILD_DIR}" && pwd)"
REPORT_DIR="${BUILD_DIR}/gcov"
mkdir -p "${REPORT_DIR}"
# Remove reports from translation units that disappeared since the last run.
find "${REPORT_DIR}" -type f -name '*.gcov' -delete
: > "${BUILD_DIR}/coverage.txt"

count=0
while IFS= read -r -d '' profile; do
  case "${profile}" in
    */third_party/*|*/src/benchmarks/*) continue ;;
  esac
  relative="${profile#"${BUILD_DIR}/"}"
  # gcov names header reports after the source path. Isolating each input
  # preserves all instantiations for Codecov to merge instead of overwriting.
  destination="${REPORT_DIR}/${relative%.gcda}"
  mkdir -p "${destination}"
  (
    cd "${destination}"
    gcov -pb "${profile}"
  ) >> "${BUILD_DIR}/coverage.txt"
  count=$((count + 1))
done < <(find "${BUILD_DIR}" -type f -name '*.gcda' -print0 | sort -z)

if ((count == 0)); then
  echo "No coverage profiles found in ${BUILD_DIR}" >&2
  exit 1
fi
echo "Collected ${count} translation units into ${REPORT_DIR}"
echo "gcov summary written to ${BUILD_DIR}/coverage.txt"
