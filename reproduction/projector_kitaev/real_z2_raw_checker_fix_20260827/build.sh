#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd);here=$(cd "$(dirname "$0")" && pwd)
set +u;source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh >/dev/null;set -u
CXX=${CXX:-mpiicpc}; eigen=${EIGEN3_INCLUDE_DIR:-/home/sunxr/software/eigen-3.4.0}; pf=/home/sunxr/new-pfqmc-main/inc/pfapack; origin=/home/sunxr/new-pfqmc-main/reproduction/projector_kitaev/average_sign_origin_tests
mkdir -p "$here/bin" "$here/regression_results"
commit=$(git --git-dir="$root/.git" --work-tree="$root" rev-parse HEAD)
common=(-mkl -O2 -std=c++17 -DPFQMC_SCALE_SAFE_UDT "-DPFQMC_SOURCE_COMMIT=\"$commit\"" -I"$eigen" -I/home/sunxr/boost_1_70_0 -I"$root/inc" -I"$root/reproduction/projector_kitaev" -I"$here")
core=("$root/src/pfqmc.cpp" "$root/src/skewMatUtils.cpp" "$pf/c_interface/libcpfapack.a" "$pf/fortran/libpfapack.a")
"$CXX" "${common[@]}" "$here/projector_real_z2_driver.cpp" "${core[@]}" -o "$here/bin/projector_real_z2_driver"
"$CXX" "${common[@]}" "$here/generic_complex_regression.cpp" "${core[@]}" -o "$here/bin/generic_complex_regression"
"$CXX" "${common[@]}" "$root/reproduction/projector_kitaev/projector_driver.cpp" "${core[@]}" -o "$here/bin/projector_driver"
"$CXX" "${common[@]}" -DPFQMC_TEST_FORCE_ZERO_AVERAGE_SIGN "$root/reproduction/projector_kitaev/projector_driver.cpp" "${core[@]}" -o "$here/bin/projector_zero_sign"
for name in finite_t_sign_driver projector_origin_driver exact_sign_enumeration_driver;do "$CXX" "${common[@]}" "$origin/$name.cpp" "${core[@]}" -o "$here/bin/$name";done
"$CXX" "${common[@]}" "$root/reproduction/driven_kitaev/driven_driver.cpp" "${core[@]}" -o "$here/bin/driven_driver"
"$CXX" "${common[@]}" "$root/reproduction/driven_kitaev/static_contour_compare.cpp" "${core[@]}" -o "$here/bin/static_contour_compare"
export PFQMC_TEST_BUILD_DIR="$here/core_build" PFQMC_PFAPACK_DIR="$pf" EIGEN3_INCLUDE_DIR="$eigen" CXX
bash "$root/reproduction/projector_kitaev/core_regression_tests/build.sh"
sha256sum "$here"/bin/* "$here"/core_build/bin/* > "$here/executable_sha256.txt"
printf 'source_branch=%s\nsource_commit=%s\ncondition_aware_ratio=false\nleft_recovery=false\n' "$(git --git-dir="$root/.git" --work-tree="$root" symbolic-ref --short HEAD)" "$commit" > "$here/build_provenance.txt"
