#!/usr/bin/env bash
set -euo pipefail
root=/home/sunxr/PfQMC-main
: "${EIGEN3_INCLUDE_DIR:?set EIGEN3_INCLUDE_DIR}"
mpiicpc -mkl -O2 -std=c++17 -I"$EIGEN3_INCLUDE_DIR" -I"$root/inc" -I"$root/reproduction/projector_kitaev" \
  "$root/reproduction/projector_kitaev/validation/static_guard_driver.cpp" "$root/src/pfqmc.cpp" "$root/src/skewMatUtils.cpp" \
  "$root/inc/pfapack/c_interface/libcpfapack.a" "$root/inc/pfapack/fortran/libpfapack.a" \
  -o "$root/reproduction/projector_kitaev/validation/static_guard_driver"
