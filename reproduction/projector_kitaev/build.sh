#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
: "${EIGEN3_INCLUDE_DIR:?Set EIGEN3_INCLUDE_DIR to the directory containing Eigen/}"
CXX=${CXX:-mpiicpc}
MKLFLAG=${MKLFLAG:--mkl}
"$CXX" $MKLFLAG -O2 -std=c++17 -DPFQMC_SCALE_SAFE_UDT -I"$EIGEN3_INCLUDE_DIR" -I"$root/inc" \
  "$root/reproduction/projector_kitaev/projector_driver.cpp" "$root/src/pfqmc.cpp" "$root/src/skewMatUtils.cpp" \
  "$root/inc/pfapack/c_interface/libcpfapack.a" "$root/inc/pfapack/fortran/libpfapack.a" \
  -o "$root/reproduction/projector_kitaev/projector_driver"
