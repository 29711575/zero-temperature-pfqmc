#!/usr/bin/env python3
"""Summarize the short FTS-like multiprecision fallback rate campaign."""
import csv, json
from pathlib import Path

root = Path("results/mp_production_rate/runs")
rows = []
for path in sorted(root.glob("*/result.json")):
    data = json.loads(path.read_text())
    data["run_dir"] = str(path.parent)
    data["total_proposals"] = data.get("proposal_attempt_count", round(
        data["pre_decision_rebuild_count"] / data["guard_trigger_frequency"]
        if data["guard_trigger_frequency"] else 0))
    data["mp_trigger_fraction"] = (data["multiprecision_fallback_count"] / data["total_proposals"]
                                   if data["total_proposals"] else 0.0)
    data["adaptive_guard_trigger_fraction"] = (data["pre_decision_rebuild_count"] / data["total_proposals"]
                                                if data["total_proposals"] else 0.0)
    rows.append(data)
fields = ["run_dir", "L", "Vf", "drive_rate", "seed", "multiprecision_fallback",
          "total_proposals", "pre_decision_rebuild_count", "multiprecision_fallback_count",
          "mp_trigger_fraction", "adaptive_guard_trigger_fraction", "multiprecision_condition_samples",
          "multiprecision_condition_p99", "multiprecision_condition_p999", "multiprecision_condition_max",
          "runtime_seconds", "S_pi", "R_cdw", "average_sign", "acceptance"]
with open("results/mp_production_rate/summary.csv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fields, extrasaction="ignore")
    writer.writeheader(); writer.writerows(rows)
print(json.dumps({"runs": len(rows), "summary": "results/mp_production_rate/summary.csv"}))
