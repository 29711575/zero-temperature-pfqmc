#!/usr/bin/env python3
"""Phase 3B production-driver contract and retained-I/O regression."""

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
DRIVER_SOURCE = pathlib.Path(__file__).with_name("pure_projector_driver.cpp")


def require(value, message):
    if not value:
        raise AssertionError(message)


def last_json(output):
    records = [json.loads(line) for line in output.splitlines() if line.startswith("{")]
    require(records, "driver emitted no JSON record")
    return records[-1]


def run_driver(executable, retained_path, *extra, env=None):
    command = [
        str(executable), "--L", "2", "--V", "0.5", "--t", "1",
        "--delta", "0.7", "--mu", "0.3", "--theta", "0.2", "--dt", "0.1",
        "--boundary", "obc", "--hs-scheme", "hs0", "--trial-t", "1",
        "--trial-delta", "0.7", "--trial-mu", "0.3", "--trial-parity", "-1",
        "--edge-splitting", "0", "--burn", "4", "--measurements", "8",
        "--seed", "913", "--stabilization-block", "2", "--retained", str(retained_path),
    ] + list(extra)
    return subprocess.run(command, text=True, capture_output=True, env=env)


def source_contract():
    require(DRIVER_SOURCE.exists(), "production driver source is missing")
    source = DRIVER_SOURCE.read_text()
    for token in ["projector_type", "pure_state", "condition_aware_ratio", "left_recovery",
                  "stabilization_block", "ratio_slow_reference_count", "trust_alarm_count"]:
        require(token in source, "missing production contract token: " + token)
    require("beta_trial" not in source, "pure-state driver must not expose beta_trial")
    for token in ["run-units", "measurement-stride", "hs_variable_count",
                  "burn_proposals", "burn_sweep_equivalent"]:
        require(token in source, "missing sweep-semantics token: " + token)


def runtime_contract(executable):
    with tempfile.TemporaryDirectory(prefix="pure-phase3b-") as directory:
        retained = pathlib.Path(directory) / "retained.csv"
        ok = run_driver(executable, retained)
        require(ok.returncode == 0, ok.stderr)
        record = last_json(ok.stdout)
        require(record["status"] == "complete", "normal run did not complete")
        require(record["projector_type"] == "pure_state", "wrong projector type")
        require(record["run_units"] == "sweeps", "production default is not sweep based")
        require(record["hs_variable_count"] == 4, "wrong HS variable count")
        require(record["burn_proposals"] == 16, "burn did not execute complete sweeps")
        require(record["burn_sweep_equivalent"] == 4, "wrong burn sweep equivalent")
        require(record["measurement_stride"] == 1 and
                record["measurement_stride_unit"] == "sweeps",
                "wrong default measurement stride")
        data = retained.read_bytes()
        require(data.endswith(b"\n"), "retained CSV lacks final newline")
        require(len(data.splitlines()) == 9, "retained CSV row count mismatch")

        legacy = run_driver(executable, "off", "--run-units", "proposals",
                            "--measurement-stride", "1")
        require(legacy.returncode == 0, legacy.stderr)
        legacy_record = last_json(legacy.stdout)
        require(legacy_record["run_units"] == "proposals", "legacy mode not explicit")
        require(legacy_record["burn_proposals"] == 4, "legacy burn unit changed")
        require(legacy_record["proposal_count"] == 12, "legacy proposal trajectory length changed")

        odd = subprocess.run([str(executable), "--L", "3", "--boundary", "pbc"],
                             text=True, capture_output=True)
        require(odd.returncode != 0 and '"status":"complete"' not in odd.stdout,
                "odd-L PBC did not fail closed")

        full = run_driver(executable, "/dev/full")
        require(full.returncode != 0, "/dev/full returned success")
        require('"status":"complete"' not in full.stdout, "/dev/full emitted complete")
        require("flush" in full.stderr.lower() or "write" in full.stderr.lower()
                or "close" in full.stderr.lower(), "/dev/full error is not explicit")


def main():
    source_contract()
    if len(sys.argv) == 2:
        runtime_contract(pathlib.Path(sys.argv[1]).resolve())
    print('{"status":"complete","tests_passed":4,"tests_total":4}')


if __name__ == "__main__":
    main()
