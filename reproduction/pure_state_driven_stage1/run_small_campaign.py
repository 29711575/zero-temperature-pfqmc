#!/usr/bin/env python3
"""Run the bounded L=4 stage-1 QMC matrix with modest CPU parallelism."""

import argparse
import concurrent.futures
import csv
import json
from pathlib import Path
import subprocess


def token(value):
    return str(value).replace(".", "p").replace("-", "m")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--jobs", type=int, default=6)
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    tasks = []
    with open(args.campaign, newline="") as handle:
        for group, row in enumerate(csv.DictReader(handle)):
            for seed_index in range(int(row["seeds"])):
                seed = 810000 + 100 * group + seed_index
                stem = (f"{row['protocol']}_tau{token(row['tau_f'])}_"
                        f"dt{token(row['Delta_tau'])}_seed{seed}")
                tasks.append((row, seed, stem))

    def run(task):
        row, seed, stem = task
        measurement = output / (stem + ".csv")
        summary = output / (stem + ".json")
        stderr = output / (stem + ".stderr")
        command = [args.executable, "--L", "4", "--boundary", "pbc",
            "--hs", "0", "--dt", row["Delta_tau"], "--parity", "-1",
            "--burn", row["burn"], "--measurements", row["measurements"],
            "--block", "4", "--seed", str(seed), "--output",
            str(measurement), "--protocol", row["driver_protocol"],
            "--Vi", row["V_i"], "--tau", row["tau_f"]]
        if row["driver_protocol"] == "quench":
            command += ["--Vf", row["V_f"]]
        else:
            command += ["--R", row["R"]]
        completed = subprocess.run(command, text=True, capture_output=True)
        summary.write_text(completed.stdout)
        stderr.write_text(completed.stderr)
        if completed.returncode:
            raise RuntimeError(f"{stem} failed: {completed.stderr.strip()}")
        parsed = json.loads(completed.stdout)
        if parsed["status"] != "complete":
            raise RuntimeError(f"{stem} did not complete")
        return {"protocol": row["protocol"], "tau_f": row["tau_f"],
            "Delta_tau": row["Delta_tau"], "V_i": row["V_i"],
            "V_f": row["V_f"], "R": row["R"], "seed": seed,
            "measurement_csv": measurement.name,
            "summary_json": summary.name, "stderr": stderr.name}

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        records = list(pool.map(run, tasks))
    with open(output / "manifest.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
