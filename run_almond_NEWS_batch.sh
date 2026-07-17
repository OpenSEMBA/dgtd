#!/usr/bin/env bash
# Sequential NEWS incidence batch for NASA Almond (5° and 10°).

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${ROOT_DIR}/build/gnu-release-mpi/bin/opensemba_dgtd"
BASE="${ROOT_DIR}/testData/maxwellInputs"

CASES=(
  "3D_Nasa_Almond_G2_25cm_5GHz_N5"
  "3D_Nasa_Almond_G2_25cm_5GHz_E5"
  "3D_Nasa_Almond_G2_25cm_5GHz_W5"
  "3D_Nasa_Almond_G2_25cm_5GHz_S5"
  "3D_Nasa_Almond_G2_25cm_5GHz_N10"
  "3D_Nasa_Almond_G2_25cm_5GHz_E10"
  "3D_Nasa_Almond_G2_25cm_5GHz_W10"
  "3D_Nasa_Almond_G2_25cm_5GHz_S10"
)

if [[ ! -x "${BIN}" ]]; then
  echo "error: missing binary: ${BIN}" >&2
  echo "Build gnu-release-cuda first, or set BIN to another preset binary." >&2
  exit 1
fi

for case in "${CASES[@]}"; do
  json="${BASE}/${case}/${case}.json"
  if [[ ! -f "${json}" ]]; then
    echo "error: missing case input: ${json}" >&2
    exit 1
  fi
  echo "=== Running ${case} ==="
  "${BIN}" -i "${json}"
done

echo "=== All NEWS almond cases finished ==="
