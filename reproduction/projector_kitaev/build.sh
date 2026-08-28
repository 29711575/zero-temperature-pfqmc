#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
: "${EIGEN3_INCLUDE_DIR:?Set EIGEN3_INCLUDE_DIR to the directory containing Eigen/}"
CXX=${CXX:-mpiicpc}
MKLFLAG=${MKLFLAG:--mkl}
OUTPUT=${OUTPUT:-$root/reproduction/projector_kitaev/projector_driver}
PFAPACK_ROOT=${PFAPACK_ROOT:-$root/inc/pfapack}
mkdir -p "$(dirname "$OUTPUT")"
"$CXX" $MKLFLAG -O2 -std=c++17 -DPFQMC_SCALE_SAFE_UDT -I"$EIGEN3_INCLUDE_DIR" -I"$root/inc" \
  "$root/reproduction/projector_kitaev/projector_driver.cpp" "$root/src/pfqmc.cpp" "$root/src/skewMatUtils.cpp" \
  "$PFAPACK_ROOT/c_interface/libcpfapack.a" "$PFAPACK_ROOT/fortran/libpfapack.a" \
  -o "$OUTPUT"
