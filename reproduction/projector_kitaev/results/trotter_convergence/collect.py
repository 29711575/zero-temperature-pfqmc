#!/usr/bin/env python3
import csv, json, math, sys
from pathlib import Path

root = Path(__file__).resolve().parent
manifest = list(csv.DictReader((root / "manifest.csv").open(encoding="utf-8")))
rows = []
for task in manifest:
    result_path = root.parent.parent / task["output_dir"] / "result.json"
    if not result_path.exists():
        continue
    try:
        with result_path.open(encoding="utf-8") as handle:
            result = json.load(handle)
    except (OSError, json.JSONDecodeError):
        continue
    numeric_keys = {"L", "V", "theta", "beta_trial", "dt", "boundary", "hs_scheme", "seed", "burn", "measurements"}
    for key in numeric_keys:
        try:
            matches = abs(float(result.get(key)) - float(task[key])) <= 1e-12
        except (TypeError, ValueError):
            matches = False
        if not matches:
            raise RuntimeError(f"parameter mismatch for {result_path}: {key}")
    expected_guard = task["adaptive_guard"] == "true"
    if "adaptive_guard" in result and bool(result["adaptive_guard"]) != expected_guard:
        raise RuntimeError(f"guard mismatch for {result_path}")
    combined = dict(task)
    combined.update(result)
    combined["result_path"] = str(result_path)
    rows.append(combined)
all_keys = {key for row in rows for key in row}
preferred = ["task_id", "campaign", "L", "V", "theta", "beta_trial", "dt", "delta", "mu", "boundary", "hs_scheme", "seed", "burn", "measurements", "threads", "adaptive_guard", "guard_threshold", "multiprecision_fallback", "full_rebuild_diagnostic", "diagnostic_stride", "stabilization_interval", "sign_recompute_stride", "S_pi", "S_pi_err", "S_pi_dq", "S_pi_dq_err", "R_cdw", "R_cdw_err", "average_sign", "acceptance", "runtime_seconds", "guard_trigger_fraction", "min_update_denominator", "sign_corrections", "max_sign_imag", "max_observable_imag", "result_path"]
keys = [key for key in preferred if key in all_keys] + sorted(all_keys - set(preferred))
with (root / "collected.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
    writer.writeheader(); writer.writerows(rows)
print(f"collected={len(rows)} expected={len(manifest)}")
