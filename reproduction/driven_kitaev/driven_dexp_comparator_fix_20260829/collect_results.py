#!/usr/bin/env python3
import csv
import json
import math
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def load_json(name):
    with (HERE / name).open() as handle:
        return json.load(handle)


def trace_max(name, column):
    values = []
    with (HERE / name).open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["region"] != "reset":
                values.append(float(row[column]))
    return max(values)


before = load_json("pathological_before.json")
after = load_json("pathological_fixed.json")
static = load_json("static_comparison.json")
smoke = load_json("driven_smoke.json")
with (HERE / "synthetic_nonzero_dexp.csv").open(newline="") as handle:
    synthetic = next(csv.DictReader(handle))

green_scale = max(float(after["max_green_abs_error"]), 1e-14)
metrics = [
    ("max_ratio_complex_abs_error", float(before["max_ratio_complex_abs_error"]),
     float(after["max_ratio_complex_abs_error"]), 10.0 * green_scale,
     "after <= 10*max(after Green error,1e-14)"),
    ("max_ratio_magnitude_abs_error", float(before["max_ratio_magnitude_abs_error"]),
     float(after["max_ratio_magnitude_abs_error"]), 10.0 * green_scale,
     "after <= 10*max(after Green error,1e-14)"),
    ("max_green_abs_error", float(before["max_green_abs_error"]),
     float(after["max_green_abs_error"]),
     float(before["max_green_abs_error"]) + 1e-15,
     "after unchanged from before (absolute tolerance 1e-15)"),
    ("max_udt_solve_residual", trace_max("pathological_before_trace.csv", "udt_solve_residual"),
     trace_max("pathological_fixed_trace.csv", "udt_solve_residual"),
     10.0 * green_scale, "after <= 10*max(after Green error,1e-14)"),
    ("max_udt_core_condition", trace_max("pathological_before_trace.csv", "udt_core_condition"),
     trace_max("pathological_fixed_trace.csv", "udt_core_condition"), math.inf,
     "diagnostic only; exponent-aware scaled core"),
    ("max_udt_d_spread", trace_max("pathological_before_trace.csv", "udt_d_spread"),
     trace_max("pathological_fixed_trace.csv", "udt_d_spread"), math.inf,
     "diagnostic only; after is full D*2^Dexp spread"),
]

with (HERE / "before_after_comparator.csv").open("w", newline="") as handle:
    out = csv.writer(handle, lineterminator="\n")
    out.writerow(["metric", "before", "after", "after_to_before", "criterion", "status"])
    for name, old, new, limit, criterion in metrics:
        ratio = new / old if old != 0 else (0.0 if new == 0 else math.inf)
        out.writerow([name, format(old, ".17g"), format(new, ".17g"),
                      format(ratio, ".17g"), criterion,
                      "PASS" if new <= limit else "FAIL"])

static_error = max(float(value) for key, value in static.items() if key.endswith("diff"))
changed = subprocess.check_output(
    ["git", "diff", "--name-only", "HEAD"], cwd=str(ROOT), text=True).splitlines()
production_changes = [path for path in changed if path.startswith(("src/", "inc/")) or
                      path == "reproduction/driven_kitaev/driven_driver.cpp" or
                      path == "main.cpp"]
summary = [
    ("synthetic_nonzero_Dexp_reconstruction", synthetic["status"],
     synthetic["fixed_reconstruction_relative_error"], "<=1e-15"),
    ("comparator_reconstruction_vs_UDT_onePlusInv",
     "PASS" if float(synthetic["comparator_vs_onePlusInv_relative_error"]) <= 1e-13 else "FAIL",
     synthetic["comparator_vs_onePlusInv_relative_error"], "<=1e-13"),
    ("pathological_ratio_error_scale", "PASS" if metrics[0][2] <= metrics[0][3] else "FAIL",
     format(metrics[0][2], ".17g"), metrics[0][4]),
    ("pathological_proxy_error_scale", "PASS" if metrics[3][2] <= metrics[3][3] else "FAIL",
     format(metrics[3][2], ".17g"), metrics[3][4]),
    ("pathological_Green_trajectory_unchanged",
     "PASS" if abs(float(after["max_green_abs_error"])-
                    float(before["max_green_abs_error"])) <= 1e-15 else "FAIL",
     format(float(after["max_green_abs_error"]), ".17g"),
     "same as before within absolute tolerance 1e-15"),
    ("driven_static_comparison", "PASS" if static_error <= 1e-13 else "FAIL",
     format(static_error, ".17g"), "all *_diff <=1e-13"),
    ("driven_smoke", "PASS" if smoke.get("status") == "complete" else "FAIL",
     smoke.get("status", "missing"), "status=complete"),
    ("production_source_unchanged", "PASS" if not production_changes else "FAIL",
     ";".join(production_changes) if production_changes else "none", "no src/inc/main/driver diff"),
]
with (HERE / "regression_summary.csv").open("w", newline="") as handle:
    out = csv.writer(handle, lineterminator="\n")
    out.writerow(["regression", "status", "value", "criterion"])
    out.writerows(summary)

failed = [name for name, status, _, _ in summary if status != "PASS"]
print("overall=" + ("PASS" if not failed else "FAIL"))
if failed:
    print("failed=" + ",".join(failed))
    sys.exit(1)
