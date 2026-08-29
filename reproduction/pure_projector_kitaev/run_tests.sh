#!/usr/bin/env bash
set -euo pipefail
: "${PURE_PROJECTOR_BUILD_DIR:?set PURE_PROJECTOR_BUILD_DIR}"
: "${PURE_PROJECTOR_TEST_OUTPUT_DIR:?set PURE_PROJECTOR_TEST_OUTPUT_DIR}"
jobs=${PURE_PROJECTOR_TEST_JOBS:-4}
binary="$PURE_PROJECTOR_BUILD_DIR/phase1_core_test"
out="$PURE_PROJECTOR_TEST_OUTPUT_DIR"
mkdir -p "$out"
tests=(single_site identity_propagation random_gaussian finite_lambda
       kitaev_zero_mode_parity propagation_thin_qr)
printf '%s\n' "${tests[@]}" | xargs -P "$jobs" -I{} bash -c '
  name=$1; binary=$2; out=$3
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "$binary" "$name" \
    >"$out/$name.json" 2>"$out/$name.stderr"
' _ {} "$binary" "$out"
python3 - "$out" <<'PY'
import csv, json, pathlib, sys
out = pathlib.Path(sys.argv[1])
names = ["single_site", "identity_propagation", "random_gaussian", "finite_lambda",
         "kitaev_zero_mode_parity", "propagation_thin_qr"]
rows = []
for name in names:
    with (out / (name + ".json")).open() as handle:
        record = json.load(handle)
    if record.get("status") != "PASS":
        raise SystemExit(name + " did not pass")
    rows.append((name, "PASS", json.dumps(record, sort_keys=True, separators=(",", ":"))))
with (out / "phase1_test_results.csv").open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(("test", "status", "details_json"))
    writer.writerows(rows)
with (out / "phase1_test_results.json").open("w") as handle:
    json.dump({"status": "PASS", "tests": len(rows)}, handle, sort_keys=True)
    handle.write("\n")
print("phase1_tests=PASS tests=6 parallel_jobs=" + str(len(names)))
PY
