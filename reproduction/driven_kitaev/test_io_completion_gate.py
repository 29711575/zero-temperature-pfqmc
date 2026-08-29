#!/usr/bin/env python3
"""Regression for retained CSV completion gating in driven_driver."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


ARGS = ["4", "0", "1", "1", "1", "1", "0.1", "1", "0", "1", "2", "30"]


def run(driver, directory, seed, bins_path, timeseries=None):
    directory.mkdir(parents=True, exist_ok=True)
    command = [str(driver)] + ARGS + [str(seed)]
    if timeseries is not None:
        command += [str(timeseries), "1"]
    environment = os.environ.copy()
    environment.update({
        "MKL_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "PFQMC_BIN_RECORDS_PATH": str(bins_path),
        "PFQMC_CODE_VERSION": "io-completion-gate-regression",
    })
    completed = subprocess.run(
        command, cwd=str(directory), env=environment,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    (directory / "stdout.log").write_bytes(completed.stdout)
    (directory / "stderr.log").write_bytes(completed.stderr)
    return completed


def csv_shape(path, rows, columns):
    data = path.read_bytes()
    if not data.endswith(b"\n"):
        raise AssertionError("missing final newline: " + str(path))
    parsed = list(csv.reader(data.decode("utf-8").splitlines()))
    if len(parsed) != rows or any(len(row) != columns for row in parsed):
        raise AssertionError("unexpected CSV shape: " + str(path))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physics_record(path):
    record = json.loads(path.read_text())
    record.pop("runtime_seconds", None)
    record.pop("bin_records_path", None)
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--driver", required=True, type=Path)
    parser.add_argument("--baseline-driver", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    options = parser.parse_args()
    driver = options.driver.resolve()
    baseline = options.baseline_driver.resolve()
    output = options.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = []

    def check(name, action):
        try:
            value = action()
            details = "" if value is None or value is True else str(value)
            results.append((name, "PASS", details))
        except Exception as error:
            results.append((name, "FAIL", str(error)))

    baseline_dir = output / "baseline_dual"
    fixed_dir = output / "fixed_dual"
    single_dir = output / "fixed_single"
    baseline_run = run(baseline, baseline_dir, 770007,
                       baseline_dir / "bins.csv", "timeseries.csv")
    fixed_run = run(driver, fixed_dir, 770007,
                    fixed_dir / "bins.csv", "timeseries.csv")
    single_run = run(driver, single_dir, 770008, single_dir / "bins.csv")

    check("normal_dual_stream_exit_and_status", lambda: (
        fixed_run.returncode == 0 and
        json.loads(fixed_run.stdout)["status"] == "complete" and
        "return=0; status=complete") or
        (_ for _ in ()).throw(AssertionError("dual-stream run did not complete")))
    check("normal_single_stream_exit_and_status", lambda: (
        single_run.returncode == 0 and
        json.loads(single_run.stdout)["status"] == "complete" and
        "return=0; status=complete") or
        (_ for _ in ()).throw(AssertionError("single-stream run did not complete")))
    def check_normal_csvs():
        csv_shape(fixed_dir / "timeseries.csv", 31, 8)
        csv_shape(fixed_dir / "bins.csv", 16, 5)
        csv_shape(single_dir / "bins.csv", 16, 5)
        return "timeseries=31x8; dual_bins=16x5; single_bins=16x5; final_newlines=yes"

    check("csv_parseability_row_count_final_newline", check_normal_csvs)
    check("trajectory_hash_matches_baseline", lambda: (
        digest(baseline_dir / "timeseries.csv") ==
        digest(fixed_dir / "timeseries.csv") and
        "sha256=" + digest(fixed_dir / "timeseries.csv")) or
        (_ for _ in ()).throw(AssertionError("trajectory hashes differ")))
    check("bin_csv_hash_matches_baseline", lambda: (
        digest(baseline_dir / "bins.csv") == digest(fixed_dir / "bins.csv") and
        "sha256=" + digest(fixed_dir / "bins.csv")) or
        (_ for _ in ()).throw(AssertionError("bin CSV hashes differ")))
    check("physics_json_matches_baseline", lambda: (
        baseline_run.returncode == 0 and
        physics_record(baseline_dir / "stdout.log") ==
        physics_record(fixed_dir / "stdout.log") and
        "all physics/provenance fields match (excluding runtime and output path)") or
        (_ for _ in ()).throw(AssertionError("physics JSON differs")))

    failure_cases = [
        ("dev_full_bin", "/dev/full", None, "retained CSV flush failure"),
        ("dev_full_timeseries", None, "/dev/full", "retained CSV flush failure"),
        ("unwritable_bin", "/root/retained.csv", None, "cannot open bin-record output"),
        ("unwritable_timeseries", None, "/root/trajectory.csv",
         "cannot open time-series output"),
    ]
    for index, (name, bins, timeseries, expected) in enumerate(failure_cases):
        directory = output / name
        bins_path = bins or str(directory / "bins.csv")
        completed = run(driver, directory, 770009 + index, bins_path, timeseries)
        check(name, lambda completed=completed, expected=expected: (
            completed.returncode != 0 and
            b'"status":"complete"' not in completed.stdout and
            expected.encode("utf-8") in completed.stderr and
            "return={}; stdout_complete=no; stderr={}".format(
                completed.returncode, completed.stderr.decode("utf-8").strip())) or
            (_ for _ in ()).throw(AssertionError(
                "failure was not gated with the expected stderr")))

    options.summary.parent.mkdir(parents=True, exist_ok=True)
    with options.summary.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["test", "status", "details"])
        writer.writerows(results)
    failures = [row for row in results if row[1] != "PASS"]
    for row in results:
        print(",".join(row))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
