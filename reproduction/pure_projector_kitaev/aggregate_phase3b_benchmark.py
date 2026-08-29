#!/usr/bin/env python3
import csv
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1]).resolve()
output = pathlib.Path(sys.argv[2]).resolve()
manifest = list(csv.DictReader((root / "manifest.csv").open()))
fields = ["task", "L", "V", "boundary", "hs", "seed", "status", "average_z2",
          "average_z2_bin_error", "S_pi", "S_pi_dq", "R_CDW", "energy",
          "fermion_parity", "acceptance", "minimum_overlap_rcond",
          "green_fast_rebuild_relative_error_max", "ratio_reference_relative_error_max",
          "rebuild_count", "ratio_slow_reference_count", "trust_alarm_count",
          "slow_reference_failure_count", "first_failure_proposal", "final_hs_hash",
          "final_rng_hash", "source_commit", "executable_sha256", "exit_code", "parsed"]
rows = []
for item in manifest:
    directory = root / ("task_%02d" % int(item["task"]))
    parsed = None
    if (directory / "result.jsonl").exists():
        for line in (directory / "result.jsonl").read_text().splitlines():
            if line.startswith("{"):
                parsed = json.loads(line)
    code = (directory / "exit_code.txt").read_text().strip() if (directory / "exit_code.txt").exists() else "missing"
    row = {key: item.get(key, "") for key in fields}
    if parsed:
        for key in fields:
            if key in parsed: row[key] = parsed[key]
    row["exit_code"] = code; row["parsed"] = int(parsed is not None)
    rows.append(row)
with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)

complete = [r for r in rows if r["status"] == "complete" and str(r["exit_code"]) == "0"]
failures = [r for r in complete if int(r["slow_reference_failure_count"]) != 0]
obc_bad = [r for r in complete if r["boundary"] == "obc" and float(r["average_z2"]) != 1]
print(json.dumps({"status": "complete" if len(complete) == 36 and not failures and not obc_bad else "incomplete",
                  "tasks_complete": len(complete), "slow_reference_failures": len(failures),
                  "obc_sign_failures": len(obc_bad)}, separators=(",", ":")))
