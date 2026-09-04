#!/usr/bin/env python3
"""Joint blocked-jackknife aggregation for the small driven-stage-1 runs."""

import argparse
import csv
import math
from pathlib import Path


RAW = ["energy", "S_pi", "S_pi_dq", "fermion_parity",
       "G_0_1_real", "G_0_1_imag", "G_0_L_real", "G_0_L_imag"]


def read_rows(path):
    with open(path, newline="") as handle:
        return [{key: float(value) if key not in {"seed", "measurement",
                "configuration_hash"} else value
                 for key, value in row.items()}
                for row in csv.DictReader(handle)]


def block_sums(runs, width):
    blocks = []
    for rows in runs:
        for begin in range(0, len(rows), width):
            chunk = rows[begin:begin + width]
            if len(chunk) != width:
                continue
            item = {"den": sum(row["physical_sign"] for row in chunk)}
            for name in RAW:
                item[name] = sum(row["physical_sign"] * row[name]
                                 for row in chunk)
            blocks.append(item)
    if len(blocks) < 4:
        raise RuntimeError("joint jackknife requires at least four full blocks")
    return blocks


def derived(total):
    den = total["den"]
    if abs(den) < 1e-12:
        raise RuntimeError("signed denominator is unresolved")
    result = {name: total[name] / den for name in RAW}
    result["R_CDW"] = 1.0 - result["S_pi_dq"] / result["S_pi"]
    result["average_sign"] = den / total["count"]
    return result


def joint_jackknife(runs, width):
    blocks = block_sums(runs, width)
    total = {"den": sum(block["den"] for block in blocks),
             "count": len(blocks) * width}
    for name in RAW:
        total[name] = sum(block[name] for block in blocks)
    estimate = derived(total)
    deleted = []
    for block in blocks:
        reduced = {"den": total["den"] - block["den"],
                   "count": total["count"] - width}
        for name in RAW:
            reduced[name] = total[name] - block[name]
        deleted.append(derived(reduced))
    errors = {}
    count = len(deleted)
    for name in list(RAW) + ["R_CDW", "average_sign"]:
        mean = sum(sample[name] for sample in deleted) / count
        errors[name] = math.sqrt((count - 1.0) / count *
            sum((sample[name] - mean) ** 2 for sample in deleted))
    return estimate, errors, len(blocks), total["count"]


def exact_rows(path):
    result = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["protocol"], float(row["tau_f"]),
                   float(row["Delta_tau"]))
            result[key] = row
    return result


def close_key(exact, protocol, tau, dt):
    for key, row in exact.items():
        if key[0] == protocol and abs(key[1] - tau) < 1e-12 and \
                abs(key[2] - dt) < 1e-12:
            return row
    raise KeyError((protocol, tau, dt))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--ed", required=True)
    parser.add_argument("--same-contour-output", required=True)
    parser.add_argument("--block-width", type=int, default=10)
    args = parser.parse_args()
    manifest_path = Path(args.manifest)
    groups = {}
    with open(manifest_path, newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["protocol"], float(row["tau_f"]),
                   float(row["Delta_tau"]))
            groups.setdefault(key, []).append(read_rows(
                manifest_path.parent / row["measurement_csv"]))
    exact = exact_rows(args.ed)
    fields = ["protocol", "tau_f", "Delta_tau", "V_i", "V_f", "R",
              "Phi_0_parity", "seeds", "samples", "jackknife_blocks",
              "average_sign", "average_sign_error"]
    for name in ["energy", "S_pi", "R_CDW", "fermion_parity",
                 "G_0_1_real", "G_0_1_imag", "G_0_L_real", "G_0_L_imag"]:
        fields += ["qmc_" + name, "qmc_" + name + "_error",
                   "ed_" + name, name + "_abs_error", name + "_z_score"]
    with open(args.same_contour_output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields,
                                lineterminator="\n")
        writer.writeheader()
        for key in sorted(groups):
            protocol, tau, dt = key
            estimate, error, blocks, samples = joint_jackknife(
                groups[key], args.block_width)
            reference = close_key(exact, protocol, tau, dt)
            row = {"protocol": protocol, "tau_f": tau, "Delta_tau": dt,
                   "V_i": reference["V_i"], "V_f": reference["V_f"],
                   "R": reference["R"],
                   "Phi_0_parity": reference["Phi_0_parity"],
                   "seeds": len(groups[key]), "samples": samples,
                   "jackknife_blocks": blocks,
                   "average_sign": estimate["average_sign"],
                   "average_sign_error": error["average_sign"]}
            mapping = {
                "energy": "contour_energy", "S_pi": "contour_S_pi",
                "R_CDW": "contour_R_CDW",
                "fermion_parity": "contour_parity",
                "G_0_1_real": "contour_G_0_1_real",
                "G_0_1_imag": "contour_G_0_1_imag",
                "G_0_L_real": "contour_G_0_L_real",
                "G_0_L_imag": "contour_G_0_L_imag"}
            for name, column in mapping.items():
                exact_value = float(reference[column])
                difference = abs(estimate[name] - exact_value)
                row["qmc_" + name] = estimate[name]
                row["qmc_" + name + "_error"] = error[name]
                row["ed_" + name] = exact_value
                row[name + "_abs_error"] = difference
                row[name + "_z_score"] = (difference / error[name]
                    if error[name] > 0.0 else (0.0 if difference == 0.0 else math.inf))
            writer.writerow(row)

if __name__ == "__main__":
    main()
