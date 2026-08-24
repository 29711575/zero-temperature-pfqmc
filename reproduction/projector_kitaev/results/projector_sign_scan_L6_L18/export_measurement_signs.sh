#!/bin/bash
set -euo pipefail
root=/home/sunxr/PfQMC-main/reproduction/projector_kitaev
base="$root/results/projector_sign_scan_L6_L18"
printf 'L,V,seed,measurement,sign,source\n'
for L in 6 12 18; do
  manifest="$base/manifest_L${L}.csv"
  while IFS=, read -r task_id campaign row_L V theta beta dt delta mu boundary hs seed burn measurements threads guard diag stab sign_stride output_dir reuse_result; do
    [[ "$task_id" == task_id ]] && continue
    if [[ -n "$reuse_result" ]]; then
      source_dir=${reuse_result%/result.json}
      csv="$root/$source_dir/measurements.csv"
    else
      csv="$base/$output_dir/measurements.csv"
    fi
    [[ -f "$csv" ]] || { echo "missing $csv" >&2; exit 2; }
    awk -F, -v L="$L" -v V="$V" -v seed="$seed" -v src="$csv" 'NR>1 {print L "," V "," seed "," $1 "," $2 "," src}' "$csv"
  done < <(tr -d '\r' < "$manifest")
done
