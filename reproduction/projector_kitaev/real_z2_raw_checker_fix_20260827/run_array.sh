#!/usr/bin/env bash
set -euo pipefail
base=$(cd "$(dirname "$0")" && pwd);i=${PBS_ARRAY_INDEX:?};row=$(awk -F, -v i="$i" 'NR==i+2{print;exit}' "$base/manifest.csv"|tr -d '\r');IFS=, read -r task label L theta V beta dt seed burn meas old expected <<< "$row";[[ "$task" == "$i" ]]||exit 2
out="$base/results/$label";mkdir -p "$out";[[ ! -e "$out/result.json" ]]||{ echo "existing result retained" >&2;exit 3;}
set +u;source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh >/dev/null;set -u;export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
"$base/bin/projector_real_z2_driver" "$L" "$theta" "$beta" "$dt" "$V" 1 0 0 0 "$seed" "$burn" "$meas" 1 "$out/before_after.csv" 200 20 "$out/legacy_measurements.csv" > "$out/result.json.tmp" 2> "$out/stderr.log"
got=$(sha256sum "$out/legacy_measurements.csv"|awk '{print $1}');[[ "$got" == "$expected" ]]||{ echo "legacy trajectory hash mismatch expected=$expected got=$got" >&2;exit 4;}
mv "$out/result.json.tmp" "$out/result.json";sha256sum "$out"/* > "$out/sha256.txt";printf 'job_id=%s\narray_index=%s\nsource_commit=%s\nlegacy_measurements_sha256=%s\ntrajectory_hash_match=true\ncondition_aware_ratio=false\nleft_recovery=false\nsign_mode=real_z2\n' "${PBS_JOBID:-unknown}" "$i" "$(awk -F= '/source_commit/{print $2}' "$base/build_provenance.txt")" "$got" > "$out/provenance.txt"
