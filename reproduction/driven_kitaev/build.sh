#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
: "${EIGEN3_INCLUDE_DIR:?Set EIGEN3_INCLUDE_DIR}"
CXX=${CXX:-mpiicpc}; MKLFLAG=${MKLFLAG:--mkl}; BOOST_INCLUDE=${BOOST_INCLUDE:-/home/sunxr/boost_1_70_0}; OUT=${OUT:-"$root/reproduction/driven_kitaev/driven_driver"}
test -f "$BOOST_INCLUDE/boost/multiprecision/cpp_dec_float.hpp"
"$CXX" $MKLFLAG -O2 -std=c++14 -I"$EIGEN3_INCLUDE_DIR" -I"$BOOST_INCLUDE" -I"$root/inc" "$root/reproduction/driven_kitaev/driven_driver.cpp" "$root/src/pfqmc.cpp" "$root/src/skewMatUtils.cpp" "$root/inc/pfapack/c_interface/libcpfapack.a" "$root/inc/pfapack/fortran/libpfapack.a" -o "$OUT"
