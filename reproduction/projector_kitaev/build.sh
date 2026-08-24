#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
: "${EIGEN3_INCLUDE_DIR:?Set EIGEN3_INCLUDE_DIR to the directory containing Eigen/}"
CXX=${CXX:-mpiicpx}
MKLFLAG=${MKLFLAG:--qmkl}
"$CXX" $MKLFLAG -O2 -std=c++17 -I"$EIGEN3_INCLUDE_DIR" -I"$root/inc" \
  "$root/reproduction/projector_kitaev/projector_driver.cpp" "$root/src/pfqmc.cpp" "$root/src/skewMatUtils.cpp" \
  "$root/inc/pfapack/c_interface/libcpfapack.a" "$root/inc/pfapack/fortran/libpfapack.a" \
  -o "$root/reproduction/projector_kitaev/projector_driver"
