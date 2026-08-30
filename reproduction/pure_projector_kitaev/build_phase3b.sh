#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
here=$(cd "$(dirname "$0")" && pwd)
: "${EIGEN3_INCLUDE_DIR:?set EIGEN3_INCLUDE_DIR}"
: "${PFQMC_PFAPACK_DIR:?set PFQMC_PFAPACK_DIR}"
: "${BOOST_INCLUDE_DIR:?set BOOST_INCLUDE_DIR}"
build_dir=${PURE_PROJECTOR_PHASE3B_BUILD_DIR:-"$here/build_phase3b"}
CXX=${CXX:-mpiicpc}
mkdir -p "$build_dir"
common=(-mkl -O2 -std=c++17 -DPFQMC_SCALE_SAFE_UDT -I"$BOOST_INCLUDE_DIR" -I"$EIGEN3_INCLUDE_DIR" -I"$root/inc")
libraries=("$root/src/skewMatUtils.cpp" "$PFQMC_PFAPACK_DIR/c_interface/libcpfapack.a" "$PFQMC_PFAPACK_DIR/fortran/libpfapack.a")
"$CXX" "${common[@]}" "$here/pure_projector_driver.cpp" "${libraries[@]}" -o "$build_dir/pure_projector_driver"
"$CXX" "${common[@]}" "$here/pure_projector_ed.cpp" "${libraries[@]}" -o "$build_dir/pure_projector_ed"
"$CXX" "${common[@]}" "$here/phase3a_core_test.cpp" "${libraries[@]}" -o "$build_dir/phase3a_core_test"
"$CXX" "${common[@]}" "$here/phase3c_core_test.cpp" "${libraries[@]}" -o "$build_dir/phase3c_core_test"
"$CXX" "${common[@]}" "$here/phase3c_frozen_replay.cpp" "${libraries[@]}" -o "$build_dir/phase3c_frozen_replay"
"$CXX" "${common[@]}" "$here/phase3c_mirrored_validation.cpp" "${libraries[@]}" -o "$build_dir/phase3c_mirrored_validation"
"$CXX" "${common[@]}" "$here/phase3a_validation.cpp" "${libraries[@]}" -o "$build_dir/phase3a_validation"
