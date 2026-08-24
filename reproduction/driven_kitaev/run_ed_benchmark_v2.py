#!/usr/bin/env python3
"""Run/reuse the finite-contour ED/PfQMC benchmark and collect its CSV files."""

import csv
import json
import math
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "ed_benchmark_v2"
RUNS = OUT / "runs"
PYTHON = Path(os.environ.get("ED_PYTHON", "/home/asus/.venv/bin/python"))
DRIVER = Path(os.environ.get("DRIVEN_DRIVER", "/home/asus/driven_driver"))
OBS = ("S_pi", "S_pi_dq", "R_cdw")


def cases():
    rows = []
    for L in (4, 6):
        for V in (0.0, 2.0):
            rows.append(dict(kind="static", L=L, V0=V, Vf=V, rate=1.0, dt=0.1))
    for Vf in (2.0, 4.0, 5.0):
        for rate in (2.0, 1.0, 0.5):
            rows.append(dict(kind="driven_L4", L=4, V0=0.0, Vf=Vf, rate=rate, dt=0.1))
    for Vf, rate in ((2.0, 1.0), (4.0, 1.0), (5.0, 1.0), (4.0, 0.5)):
        rows.append(dict(kind="driven_L6", L=6, V0=0.0, Vf=Vf, rate=rate, dt=0.1))
    for L in (4, 6):
        rows.append(dict(kind="dt_check", L=L, V0=0.0, Vf=4.0, rate=1.0, dt=0.05))
    return rows


def tag(c):
    f = lambda x: str(x).replace(".", "p")
    return f"{c['kind']}_L{c['L']}_v0{f(c['V0'])}_vf{f(c['Vf'])}_r{f(c['rate'])}_dt{f(c['dt'])}"


def run_ed(c):
    cmd = [str(PYTHON), str(HERE / "exact_driven_reference.py"), "--L", str(c["L"]),
           "--V0", str(c["V0"]), "--Vf", str(c["Vf"]), "--rate", str(c["rate"]),
           "--dt", str(c["dt"]), "--theta-init", "6", "--beta-trial", "8"]
    return json.loads(subprocess.check_output(cmd, text=True))


def seed_for(c, index):
    kind_code = {"static": 1, "driven_L4": 2, "driven_L6": 3, "dt_check": 4}[c["kind"]]
    return 960000 + kind_code * 10000 + c["L"] * 1000 + int(c["Vf"] * 100) + int(c["rate"] * 10) + index


def run_qmc(item):
    c, index = item
    seed = seed_for(c, index)
    d = RUNS / tag(c) / f"seed_{seed}"
    result = d / "result.json"
    if result.exists():
        data = json.loads(result.read_text())
        required = data.get("burn") == 300 and data.get("measurements") == 3000 and data.get("boundary") == 1
        if required:
            return c, data, True
    d.mkdir(parents=True, exist_ok=True)
    cmd = [str(DRIVER), str(c["L"]), str(c["V0"]), str(c["Vf"]), str(c["rate"]),
           "6", "8", str(c["dt"]), "1", "0", "1", "300", "3000", str(seed)]
    proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (d / "stderr.log").write_text(proc.stderr)
    if proc.returncode:
        raise RuntimeError(f"QMC failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr}")
    data = json.loads(proc.stdout)
    result.write_text(proc.stdout.strip() + "\n")
    return c, data, False


def mean_stderr(values):
    mean = sum(values) / len(values)
    if len(values) < 2:
        return mean, math.nan
    return mean, math.sqrt(sum((x - mean) ** 2 for x in values) / (len(values) * (len(values) - 1)))


def write_csv(path, rows):
    fields = list(rows[0])
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cs = cases()
    ed = {tag(c): run_ed(c) for c in cs}

    raw = {tag(c): [] for c in cs}
    reused = 0
    workers = int(os.environ.get("QMC_WORKERS", "4"))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run_qmc, (c, i)) for c in cs for i in range(4)]
        for future in as_completed(futures):
            c, data, was_reused = future.result()
            raw[tag(c)].append(data)
            reused += was_reused
            print(f"finished {tag(c)} seed={data['seed']} reused={was_reused}", flush=True)

    keys = ("kind", "L", "V0", "Vf", "rate", "dt")
    ed_rows, qmc_rows, cmp_rows = [], [], []
    for c in cs:
        e = ed[tag(c)]
        erow = {k: c[k] for k in keys}
        erow.update({"theta_init": 6.0, "beta_trial": 8.0, "trial_slices": e["trial_slices"],
                     "initial_slices_per_side": e["initial_slices_per_side"], "n_ramp": e["n_ramp"]})
        erow.update({o: e[o] for o in OBS})
        ed_rows.append(erow)

        qrow = {k: c[k] for k in keys}
        qrow["n_seeds"] = len(raw[tag(c)])
        for o in OBS + ("average_sign",):
            qrow[o], qrow[o + "_stderr"] = mean_stderr([x[o] for x in raw[tag(c)]])
        qmc_rows.append(qrow)

        crow = {k: c[k] for k in keys}
        for o in OBS:
            delta = qrow[o] - e[o]
            crow["ED_" + o] = e[o]
            crow["QMC_" + o] = qrow[o]
            crow["QMC_" + o + "_stderr"] = qrow[o + "_stderr"]
            crow["delta_" + o] = delta
            if qrow[o + "_stderr"] > 0:
                crow["z_" + o] = delta / qrow[o + "_stderr"]
            else:
                # Free points are deterministic.  Treat roundoff-level ED/QMC
                # differences as exact agreement rather than an infinite z.
                crow["z_" + o] = 0.0 if abs(delta) < 1e-12 else math.copysign(math.inf, delta)
        crow["average_sign"] = qrow["average_sign"]
        crow["average_sign_stderr"] = qrow["average_sign_stderr"]
        cmp_rows.append(crow)

    write_csv(OUT / "ed_results.csv", ed_rows)
    write_csv(OUT / "qmc_results.csv", qmc_rows)
    write_csv(OUT / "comparison.csv", cmp_rows)

    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(9, 3))
        for ax, o in zip(axes, OBS):
            x = [r["ED_" + o] for r in cmp_rows]
            y = [r["QMC_" + o] for r in cmp_rows]
            ye = [r["QMC_" + o + "_stderr"] for r in cmp_rows]
            ax.errorbar(x, y, yerr=ye, fmt="o", ms=3)
            lo, hi = min(x + y), max(x + y)
            ax.plot([lo, hi], [lo, hi], "k--", lw=1)
            ax.set(xlabel="ED", ylabel="QMC", title=o)
        fig.tight_layout()
        fig.savefig(OUT / "parity.png", dpi=120)
        plt.close(fig)
    except ImportError:
        pass
    print(f"wrote {OUT}; reused {reused}/{len(cs)*4} QMC runs")


if __name__ == "__main__":
    main()
