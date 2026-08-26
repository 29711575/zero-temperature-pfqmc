#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
here=$(cd "$(dirname "$0")" && pwd)
: "${PFQMC_TEST_BUILD_DIR:?set an external build directory}"
: "${PFQMC_PFAPACK_DIR:?set the read-only directory containing c_interface/ and fortran/}"
: "${EIGEN3_INCLUDE_DIR:?set EIGEN3_INCLUDE_DIR}"
CXX=${CXX:-mpiicpc}
mkdir -p "$PFQMC_TEST_BUILD_DIR/bin"
common=(-mkl -O2 -std=c++17 -DPFQMC_SCALE_SAFE_UDT
        -I"$EIGEN3_INCLUDE_DIR" -I/home/sunxr/boost_1_70_0
        -I"$root/inc" -I"$root/reproduction/projector_kitaev" -I"$here")
core=("$root/src/pfqmc.cpp" "$root/src/skewMatUtils.cpp"
      "$PFQMC_PFAPACK_DIR/c_interface/libcpfapack.a"
      "$PFQMC_PFAPACK_DIR/fortran/libpfapack.a")
build_core() { "$CXX" "${common[@]}" "$here/$2" "${core[@]}" -o "$PFQMC_TEST_BUILD_DIR/bin/$1"; }
build_header() { "$CXX" "${common[@]}" "$here/$2" -o "$PFQMC_TEST_BUILD_DIR/bin/$1"; }
build_core integration_qmc integration_qmc_driver.cpp
build_core local_update_property local_update_property_driver.cpp
build_core right_boundary_green right_boundary_green_driver.cpp
build_core left_green_recovery left_green_recovery_driver.cpp
build_core reality_symmetry reality_symmetry_driver.cpp
build_header udt_stress udt_stress_driver.cpp
build_header udt_guard_stress udt_guard_stress_driver.cpp
"$CXX" "${common[@]}" \
  "$root/reproduction/projector_kitaev/regression_stress/scripts/tiny_enumeration_driver.cpp" \
  "${core[@]}" -o "$PFQMC_TEST_BUILD_DIR/bin/tiny"
"$CXX" "${common[@]}" \
  "$root/reproduction/projector_kitaev/regression_stress/scripts/gaussian_exact_driver.cpp" \
  "${core[@]}" -o "$PFQMC_TEST_BUILD_DIR/bin/gaussian"
"$CXX" "${common[@]}" -DPFQMC_TEST_FORCE_ZERO_AVERAGE_SIGN \
  "$root/reproduction/projector_kitaev/projector_driver.cpp" "${core[@]}" \
  -o "$PFQMC_TEST_BUILD_DIR/bin/projector_zero_sign"
sha256sum "$PFQMC_TEST_BUILD_DIR"/bin/*
