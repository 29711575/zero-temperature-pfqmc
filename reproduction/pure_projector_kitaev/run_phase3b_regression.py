#!/usr/bin/env python3
import argparse
import csv
import json
import os
import pathlib
import subprocess


def record(command, env=None):
    run = subprocess.run(command, text=True, capture_output=True, env=env)
    parsed = [json.loads(line) for line in run.stdout.splitlines() if line.startswith("{")]
    return run, parsed[-1] if parsed else None


def base(exe, L, boundary, hs, seed, block, retained="off", mode="fast-strict",
         theta="0.4", burn="40", measurements="80", V="0.8", audit="0"):
    zero = boundary == "obc"
    trial_parity = 1 if zero or L % 4 == 2 else -1
    return [str(exe), "--L", str(L), "--V", V, "--t", "1", "--delta", "1",
            "--mu", "0", "--theta", theta, "--dt", "0.1", "--boundary", boundary,
            "--hs-scheme", hs, "--trial-t", "1", "--trial-delta", "1",
            "--trial-mu", "0" if zero else "0.3", "--trial-parity", str(trial_parity),
            "--edge-splitting", "1e-8" if zero else "0", "--burn", burn,
            "--measurements", measurements, "--seed", str(seed),
            "--stabilization-block", str(block), "--retained", str(retained),
            "--audit-interval", audit, "--walker-mode", mode]


def require_complete(run, value, context):
    if run.returncode or not value or value.get("status") != "complete":
        raise RuntimeError(context + ": " + run.stderr + run.stdout)


def run_core(exe, out):
    audit_rows = []
    for L in (2, 4):
        results = []
        for mode in ("fast-strict", "audit-lockstep"):
            run, value = record(base(exe, L, "pbc", "hs0", 41000 + L, 2, mode=mode))
            require_complete(run, value, "trajectory " + mode)
            results.append(value)
        mismatch = (results[0]["final_hs_hash"] != results[1]["final_hs_hash"] or
                    results[0]["final_rng_hash"] != results[1]["final_rng_hash"] or
                    results[0]["average_z2"] != results[1]["average_z2"])
        audit_rows.append([L, results[0]["final_hs_hash"], results[1]["final_hs_hash"],
                           results[0]["average_z2"], results[1]["average_z2"],
                           results[0]["ratio_slow_reference_count"],
                           results[1]["ratio_slow_reference_count"],
                           results[1]["ratio_reference_relative_error_max"],
                           results[1]["green_fast_rebuild_relative_error_max"], int(mismatch)])
        if mismatch:
            raise RuntimeError("production/slow trajectory mismatch")
    with (out / "proposal_audit.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["L", "fast_hs_hash", "slow_hs_hash", "fast_average_z2",
                         "slow_average_z2", "fast_slow_reference_count",
                         "audit_slow_reference_count", "ratio_error_max", "green_error_max",
                         "trajectory_mismatch"])
        writer.writerows(audit_rows)

    block_rows, reference = [], None
    for block in (1, 2, 4, 8):
        run, value = record(base(exe, 4, "pbc", "hs0", 42004, block, audit="20"))
        require_complete(run, value, "block size")
        if reference is None:
            reference = value
        mismatch = (value["final_hs_hash"] != reference["final_hs_hash"] or
                    value["final_rng_hash"] != reference["final_rng_hash"] or
                    abs(value["S_pi"] - reference["S_pi"]) > 1e-10)
        block_rows.append([block, value["final_hs_hash"], value["S_pi"], value["energy"],
                           value["average_z2"], value["green_fast_rebuild_relative_error_max"],
                           value["ratio_slow_reference_count"], int(mismatch)])
        if mismatch:
            raise RuntimeError("block-size trajectory mismatch")
    with (out / "block_size_checks.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["block_size", "final_hs_hash", "S_pi", "energy", "average_z2",
                         "green_rebuild_error_max", "slow_reference_count", "trajectory_mismatch"])
        writer.writerows(block_rows)

    rows = []
    retained = out / "retained_roundtrip.csv"
    run, value = record(base(exe, 2, "obc", "hs0", 43002, 2, retained=retained,
                             measurements="8", burn="4"))
    require_complete(run, value, "retained round trip")
    lines = retained.read_bytes().splitlines()
    rows.append(["normal", run.returncode, 1, len(lines) == 9, 1, run.stderr.strip()])
    run, value = record(base(exe, 2, "obc", "hs0", 43003, 2, retained="/dev/full",
                             measurements="8", burn="4"))
    rows.append(["dev_full", run.returncode, int(value is not None and value.get("status") == "complete"),
                 0, int(run.returncode != 0 and not value), run.stderr.strip()])
    run, value = record(base(exe, 2, "obc", "hs0", 43004, 2, retained=out,
                             measurements="8", burn="4"))
    rows.append(["unwritable_path", run.returncode, int(value is not None), 0,
                 int(run.returncode != 0 and not value), run.stderr.strip()])
    env = os.environ.copy(); env["PFQMC_TEST_FORCE_ZERO_AVERAGE_SIGN"] = "1"
    run, value = record(base(exe, 2, "obc", "hs0", 43005, 2, retained="off",
                             measurements="8", burn="4"), env=env)
    require_complete(run, value, "zero-sign JSON")
    zero_ok = value["average_z2"] == 0 and value["S_pi"] is None and value["energy"] is None
    rows.append(["zero_sign", run.returncode, 1, 1, int(zero_ok), run.stderr.strip()])
    with (out / "io_regression.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["case", "return_code", "complete_emitted", "parseable", "passed", "stderr"])
        writer.writerows(rows)
    if not all(row[4] for row in rows):
        raise RuntimeError("I/O regression failed")


def run_ed(exe, ed_exe, out):
    rows = []
    for L, boundary, hs in ((4, "pbc", "hs0"), (6, "obc", "hs1")):
        trial_mu = "0" if boundary == "obc" else "0.3"
        split = "1e-8" if boundary == "obc" else "0"
        command = [str(ed_exe), "--L", str(L), "--V", "0.5", "--t", "1",
                   "--delta", "1", "--mu", "0", "--boundary", boundary,
                   "--trial-t", "1", "--trial-delta", "1", "--trial-mu", trial_mu,
                   "--trial-parity", "1" if boundary == "obc" or L % 4 == 2 else "-1",
                   "--edge-splitting", split]
        run, exact = record(command); require_complete(run, exact, "ED oracle")
        for theta in (2, 4, 8, 12):
            run, mc = record(base(exe, L, boundary, hs, 44000 + 10 * L + theta, 4,
                                  theta=str(theta), burn="600", measurements="2400", V="0.5",
                                  audit="200"))
            require_complete(run, mc, "theta convergence")
            rows.append([L, boundary, theta, exact["energy"], mc["energy"],
                         abs(mc["energy"] - exact["energy"]), exact["fermion_parity"],
                         mc["fermion_parity"], abs(mc["fermion_parity"] - exact["fermion_parity"]),
                         exact["S_pi"], mc["S_pi"], abs(mc["S_pi"] - exact["S_pi"]),
                         exact["R_CDW"], mc["R_CDW"], abs(mc["R_CDW"] - exact["R_CDW"]),
                         mc["average_z2"], mc["green_fast_rebuild_relative_error_max"],
                         mc["ratio_slow_reference_count"], mc["slow_reference_failure_count"]])
    with (out / "ed_theta_convergence.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["L", "boundary", "theta", "ED_energy", "MC_energy", "energy_abs_diff",
                         "ED_parity", "MC_parity", "parity_abs_diff", "ED_S_pi", "MC_S_pi",
                         "S_pi_abs_diff", "ED_R_CDW", "MC_R_CDW", "R_CDW_abs_diff",
                         "average_z2", "green_rebuild_error_max", "slow_reference_count",
                         "slow_reference_failure_count"])
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("core", "ed", "all"))
    parser.add_argument("driver", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--ed-executable", type=pathlib.Path)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    if args.mode in ("core", "all"): run_core(args.driver.resolve(), args.output.resolve())
    if args.mode in ("ed", "all"):
        if not args.ed_executable: raise RuntimeError("--ed-executable required")
        run_ed(args.driver.resolve(), args.ed_executable.resolve(), args.output.resolve())
    print(json.dumps({"status": "complete", "mode": args.mode}, separators=(",", ":")))


if __name__ == "__main__": main()
