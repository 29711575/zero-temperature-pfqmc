#!/usr/bin/env python3
import csv, glob, json, os

base = "results/hpc_convergence"
rows = []
for path in sorted(glob.glob(base + "/**/result.json", recursive=True)):
    try:
        with open(path, encoding="utf-8") as f:
            row = json.load(f)
    except (OSError, json.JSONDecodeError):
        continue
    if not isinstance(row, dict) or not row:
        continue
    if "SQ0" in row:
        row["mode"] = "finite_temperature"
        row["S_pi"] = row.get("SQ0")
        row["S_pi_dq"] = row.get("SQ0_delta")
        row["R_cdw"] = row.get("R")
        row["runtime_seconds"] = row.get("wall_time")
        row["beta"] = row.get("Beta", row.get("beta"))
    row["result_path"] = path
    rows.append(row)

preferred = ["mode", "L", "V", "theta", "beta_trial", "beta", "dt", "Dtau",
             "delta", "Delta", "mu", "boundary", "hs_scheme", "seed", "burn",
             "measurements", "threads", "S_pi", "S_pi_dq", "R_cdw", "average_sign",
             "acceptance", "runtime_seconds", "sign_recomputes", "sign_corrections",
             "max_sign_imag", "max_observable_imag", "result_path"]
all_keys = {key for row in rows for key in row}
keys = [key for key in preferred if key in all_keys]
keys += sorted(all_keys - set(keys))
if rows:
    with open(base + "/hpc_convergence_summary.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
print(f"collected={len(rows)}")
