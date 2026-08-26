#!/usr/bin/env bash
set -euo pipefail
: "${PFQMC_TEST_BUILD_DIR:?set external build directory}"
: "${PFQMC_TEST_OUTPUT_DIR:?set external output directory}"
bin="$PFQMC_TEST_BUILD_DIR/bin"
out="$PFQMC_TEST_OUTPUT_DIR"
mkdir -p "$out"
"$bin/tiny" 2 0.2 0.2 0.1 2 1 0 0 0 4 "$out/tiny.csv" >"$out/tiny.json"
"$bin/gaussian" 6 2 2 0.1 1 0 0 0 >"$out/gaussian.json"
"$bin/integration_qmc" 10 10 8 .1 2 1 0 0 0 750010 5 10 1 "$out/L10.csv" 2 10 0 10 >"$out/L10.json"
"$bin/integration_qmc" 26 39 8 .1 4 1 0 0 0 750049 2 5 1 "$out/task88.csv" 1 10 0 5 >"$out/task88.json"
"$bin/integration_qmc" 32 48 8 .1 5 1 0 0 0 750052 1 3 1 "$out/task92.csv" 1 10 0 3 >"$out/task92.json"
"$bin/local_update_property" 6 0 0 4 880006 3000 2 2 "$out/local_update.csv" >"$out/local_update.json"
"$bin/right_boundary_green" 6 0 4 860006 2 2 .1 "$out/boundary_green.csv" >"$out/boundary_green.json"
"$bin/left_green_recovery" 12 0 5 760126 12 8 .1 3 1e-10 "$out/left_recovery.csv" >"$out/left_recovery.json"
"$bin/udt_stress" 12 20 900012 "$out/udt_normal.csv" >"$out/udt_normal.json"
"$bin/udt_guard_stress" "$out/udt_guard.csv" >"$out/udt_guard.json"
"$bin/reality_symmetry" 4 2 890004 2 "$out/reality.csv" >"$out/reality.json"
"$bin/projector_zero_sign" 6 2 2 .1 2 1 0 0 0 750006 1 2 1 >"$out/zero_sign.json"
root=$(cd "$(dirname "$0")/../../.." && pwd)
PYTHON3=${PYTHON3:-python3}
if command -v "$PYTHON3" >/dev/null 2>&1; then
  "$PYTHON3" "$root/reproduction/projector_kitaev/regression_stress/scripts/same_contour_ed_general.py" \
    --L 4 --V 2 --theta .4 --beta-trial .4 --dt .1 --boundary 0 \
    >"$out/same_contour_ed.json"
else
  printf 'SKIP: set PYTHON3 to a Python 3 environment with NumPy/SciPy\n' \
    >"$out/same_contour_ed.skip.txt"
fi
