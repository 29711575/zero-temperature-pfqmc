#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
base=$(cd "$(dirname "$0")" && pwd)
export LC_ALL=C LANG=C
set +u; source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh >/dev/null; set -u
CXX=${CXX:-mpiicpc}; eigen=/home/sunxr/software/eigen-3.4.0; pf=/home/sunxr/new-pfqmc-main/inc/pfapack
mkdir -p "$base/bin" "$base/core_build" "$base/regression_results"
commit=$(git --git-dir="$root/.git" --work-tree="$root" rev-parse HEAD)
common=(-mkl -O2 -std=c++17 -DPFQMC_SCALE_SAFE_UDT "-DPFQMC_SOURCE_COMMIT=\"$commit\"" -I"$eigen" -I/home/sunxr/boost_1_70_0 -I"$root/inc" -I"$root/reproduction/projector_kitaev" -I"$root/reproduction/projector_kitaev/real_z2_raw_checker_fix_20260827")
core=("$root/src/pfqmc.cpp" "$root/src/skewMatUtils.cpp" "$pf/c_interface/libcpfapack.a" "$pf/fortran/libpfapack.a")
"$CXX" "${common[@]}" "$root/reproduction/projector_kitaev/real_z2_raw_checker_fix_20260827/projector_real_z2_driver.cpp" "${core[@]}" -o "$base/bin/projector_real_z2_driver"
"$CXX" "${common[@]}" "$root/reproduction/projector_kitaev/real_z2_raw_checker_fix_20260827/generic_complex_regression.cpp" "${core[@]}" -o "$base/bin/generic_complex_regression"
"$CXX" "${common[@]}" -DPFQMC_MP_ORACLE_SMALL_ONLY "$base/validation_oracle_driver.cpp" "${core[@]}" -o "$base/bin/validation_oracle"
"$CXX" "${common[@]}" "$root/reproduction/driven_kitaev/driven_driver.cpp" "${core[@]}" -o "$base/bin/driven_driver"
"$CXX" "${common[@]}" "$root/reproduction/driven_kitaev/static_contour_compare.cpp" "${core[@]}" -o "$base/bin/static_contour_compare"
exact=/home/sunxr/new-pfqmc-main/reproduction/projector_kitaev/average_sign_origin_tests/exact_sign_enumeration_driver.cpp
"$CXX" "${common[@]}" "$exact" "${core[@]}" -o "$base/bin/exact_sign_enumeration_driver"
PFQMC_TEST_BUILD_DIR="$base/core_build" PFQMC_PFAPACK_DIR="$pf" EIGEN3_INCLUDE_DIR="$eigen" CXX="$CXX" bash "$root/reproduction/projector_kitaev/core_regression_tests/build.sh"
sha256sum "$base"/bin/* "$base"/core_build/bin/* > "$base/executable_sha256.txt"
printf 'source_branch=fix-mp-z2-trust-gate\nsource_commit=%s\ncondition_aware_ratio=false\nleft_recovery=false\nmp_checkpoint_mutating=false\ncanonical_order=index_zero\n' "$commit" > "$base/build_provenance.txt"

