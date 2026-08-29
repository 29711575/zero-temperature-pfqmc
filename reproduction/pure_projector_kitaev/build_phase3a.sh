#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
here=$(cd "$(dirname "$0")" && pwd)
: "${EIGEN3_INCLUDE_DIR:?set EIGEN3_INCLUDE_DIR}"
: "${PFQMC_PFAPACK_DIR:?set PFQMC_PFAPACK_DIR}"
build_dir=${PURE_PROJECTOR_PHASE3A_BUILD_DIR:-"$here/build_phase3a"}
CXX=${CXX:-mpiicpc}
mkdir -p "$build_dir"
"$CXX" -mkl -O2 -std=c++17 -DPFQMC_SCALE_SAFE_UDT \
  -I"$EIGEN3_INCLUDE_DIR" -I"$root/inc" \
  "$here/phase3a_core_test.cpp" "$root/src/skewMatUtils.cpp" \
  "$PFQMC_PFAPACK_DIR/c_interface/libcpfapack.a" \
  "$PFQMC_PFAPACK_DIR/fortran/libpfapack.a" \
  -o "$build_dir/phase3a_core_test"
