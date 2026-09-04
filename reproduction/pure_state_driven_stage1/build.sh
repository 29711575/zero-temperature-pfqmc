#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
here=$(cd "$(dirname "$0")" && pwd)
: "${EIGEN3_INCLUDE_DIR:?set EIGEN3_INCLUDE_DIR}"
: "${BOOST_INCLUDE_DIR:?set BOOST_INCLUDE_DIR}"
: "${PFQMC_PFAPACK_DIR:?set PFQMC_PFAPACK_DIR}"
build=${PURE_DRIVEN_STAGE1_BUILD_DIR:-"$here/build"}
CXX=${CXX:-mpiicpc}
mkdir -p "$build"

common=(
  -mkl -O2 -std=c++17 -DPFQMC_SCALE_SAFE_UDT
  -I"$BOOST_INCLUDE_DIR" -I"$EIGEN3_INCLUDE_DIR" -I"$repo/inc"
)
libraries=(
  "$repo/src/skewMatUtils.cpp"
  "$PFQMC_PFAPACK_DIR/c_interface/libcpfapack.a"
  "$PFQMC_PFAPACK_DIR/fortran/libpfapack.a"
)

for source in driven_static_validation driven_exact_enumeration driven_qmc driven_ed; do
  "$CXX" "${common[@]}" "$here/${source}.cpp" "${libraries[@]}" \
    -o "$build/$source"
done

sha256sum "$build"/driven_static_validation \
  "$build"/driven_exact_enumeration "$build"/driven_qmc "$build"/driven_ed
