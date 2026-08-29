#!/usr/bin/env bash
set -euo pipefail
base=$(cd "$(dirname "$0")" && pwd); out=${1:?output directory required}; mkdir -p "$out"
export LC_ALL=C LANG=C OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
if [[ ${SETVARS_COMPLETED:-0} != 1 ]]; then
  set +u
  source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh >/dev/null
  set -u
fi
"$base/bin/generic_complex_regression" > "$out/generic_complex.json"
"$base/core_build/bin/sign_interface_hardening" > "$out/sign_interface.json"
"$base/core_build/bin/local_update_property" 4 1 1 2 990 10 .2 .4 "$out/local_flips.csv" > "$out/local_flips.json"
"$base/core_build/bin/right_boundary_green" 4 1 2 990 .2 .4 .1 "$out/right_boundary.csv" > "$out/right_boundary.json"
"$base/core_build/bin/reality_symmetry" 6 4 1926155102 4 "$out/reality.csv" > "$out/reality.json"
"$base/core_build/bin/integration_qmc" 4 .2 .4 .1 2 1 0 1 1 990003 1 2 1 "$out/integration.csv" 1 10 0 1 > "$out/integration.json"
"$base/core_build/bin/projector_zero_sign" 4 .2 .4 .1 2 1 0 1 1 990003 1 2 1 > "$out/zero_sign.json"
"$base/bin/exact_sign_enumeration_driver" projector 4 .2 .4 .1 2 1 0 0 0 4 "$out/exact_projector.csv" > "$out/exact_projector.json"
"$base/bin/exact_sign_enumeration_driver" finite 4 .2 .4 .1 2 1 0 0 0 4 "$out/exact_finite.csv" > "$out/exact_finite.json"
"$base/bin/static_contour_compare" 2 > "$out/driven_static_compare.json"
"$base/bin/driven_driver" 4 2 2 0 .2 .4 .1 1 0 1 1 2 990004 0 .8 > "$out/driven_smoke.json"
"$base/bin/projector_real_z2_driver" 4 .2 .4 .1 2 1 0 1 1 990003 1 2 1 "$out/obc_hs1.csv" 1 1 "$out/obc_hs1_legacy.csv" 1 > "$out/obc_hs1.json"
python3 - "$out" <<'PY'
import csv,json,pathlib,sys
p=pathlib.Path(sys.argv[1]); rows=[]
def load(n): return json.load((p/n).open())
g=load('generic_complex.json'); assert g['expected_error']<1e-12 and not g['real_z2_mode']; rows.append(('generic_complex','PASS'))
for n in ('exact_projector.json','exact_finite.json'):
 d=load(n); assert d['max_fock_phase_abs_difference']<1e-9 and d['max_fock_logabs_offset_variation']<1e-9; rows.append((n,'PASS'))
o=load('obc_hs1.json'); assert o['z2_average_sign']==1 and o['shadow_oracle_z2_mismatch_count']==0 and o['shadow_trajectory_match'] and o['mp_correction_count']==0; rows.append(('OBC_hs1','PASS'))
z=load('zero_sign.json'); assert z['status']=='complete' and z['average_sign']==0 and z['S_pi'] is None; rows.append(('zero_sign_json','PASS'))
r=load('right_boundary.json'); assert r['nonfinite']==0 and r['max_error']<1e-10
d=load('driven_smoke.json'); assert d['status']=='complete' and d['sign_mode']=='real_z2' and d['raw_sign_check_mismatch']==0 and d['mp_checkpoint_mutating'] is False and d['mp_correction_count']==0
for n in ('sign_interface.json','local_flips.json','right_boundary.json','reality.json','integration.json','driven_static_compare.json','driven_smoke.json'): rows.append((n,'PASS'))
with (p/'regression_summary.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(('test','status'));w.writerows(rows)
print('regressions=PASS')
PY
