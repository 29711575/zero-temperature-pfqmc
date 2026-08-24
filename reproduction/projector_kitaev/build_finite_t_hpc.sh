#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
: "${EIGEN3_INCLUDE_DIR:?Set EIGEN3_INCLUDE_DIR}"
icpc -mkl -O2 -std=c++17 -I"$EIGEN3_INCLUDE_DIR" -I"$root/inc" \
  "$root/reproduction/fig2b_fast/fig2b_fast_driver.cpp" \
  "$root/src/pfqmc.cpp" "$root/src/skewMatUtils.cpp" \
  "$root/inc/pfapack/c_interface/libcpfapack.a" \
  "$root/inc/pfapack/fortran/libpfapack.a" \
  -o "$root/reproduction/projector_kitaev/finite_t_driver_hpc"
