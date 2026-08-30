#!/usr/bin/env python3
import csv
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
manifest = pathlib.Path(sys.argv[2])
output = pathlib.Path(sys.argv[3])
rows = []
with manifest.open() as handle:
    records = list(csv.DictReader(handle))
for record in records:
    run = int(record["run"]); original = int(record["original_task"])
    directory = root / f"run_{run:02d}_original_{original:02d}"
    status = "missing"; value = {}; code = None
    if (directory.exists()):
        try: code = int((directory / "exit_code.txt").read_text().strip())
        except Exception: code = None
        try:
            parsed = [json.loads(line) for line in (directory / "result.jsonl").read_text().splitlines()
                      if line.startswith("{")]
            if parsed: value = parsed[-1]
        except Exception: pass
        status = "complete" if code == 0 and value.get("status") == "complete" else "failed"
    rows.append({**record, "status": status, "exit_code": code,
        "average_z2": value.get("average_z2"),
        "mp_fallback_count": value.get("mp_same_proposal_fallback_count"),
        "mp_fallback_failures": value.get("mp_same_proposal_fallback_failure_count"),
        "endpoint_green_residual_max": value.get("endpoint_rebuild_green_residual_max"),
        "z2_oracle_correction_count": value.get("z2_oracle_correction_count"),
        "final_hs_hash": value.get("final_hs_hash"), "final_rng_hash": value.get("final_rng_hash"),
        "runtime_seconds": value.get("runtime_seconds"),
        "initialization_policy": value.get("initialization_policy")})
fields = list(rows[0])
with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fields, lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
complete = sum(row["status"] == "complete" for row in rows)
failures = sum((row["mp_fallback_failures"] or 0) for row in rows)
obc_ok = all(row["average_z2"] == 1 for row in rows if row["role"] == "obc_hs1_control" and row["status"] == "complete")
print(json.dumps({"status": "complete" if complete == len(rows) and failures == 0 and obc_ok else "incomplete",
                  "complete": complete, "total": len(rows), "mp_fallback_failures": failures,
                  "obc_hs1_all_plus": obc_ok}, separators=(",", ":")))
