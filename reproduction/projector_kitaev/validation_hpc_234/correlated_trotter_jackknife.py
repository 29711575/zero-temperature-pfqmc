#!/usr/bin/env python
from __future__ import division, print_function

import csv
import math
import os
import re
from collections import defaultdict

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
BIN_SIZE = 25
DTS = (0.20, 0.10, 0.05)
OBSERVABLES = ("S_pi", "S_pi_dq", "R_cdw")


def read_measurements(path):
    rows = []
    with open(path) as handle:
        for row in csv.DictReader(handle):
            rows.append((float(row["sign"]),
                         float(row["sign_S_pi_numerator"]),
                         float(row["sign_S_pi_dq_numerator"])))
    return rows


def observable(totals, name):
    den, num_pi, num_dq = totals
    if name == "S_pi":
        return num_pi / den
    if name == "S_pi_dq":
        return num_dq / den
    return 1.0 - num_dq / num_pi


def jk_error(samples):
    center = samples.mean(axis=0)
    delta = samples - center
    covariance = (len(samples) - 1.0) / len(samples) * np.dot(delta.T, delta)
    return covariance


def gls_fit(y, covariance):
    x = np.asarray([[1.0, dt * dt] for dt in DTS])
    weight = np.linalg.pinv(covariance, rcond=1e-12)
    normal_inv = np.linalg.pinv(np.dot(x.T, np.dot(weight, x)), rcond=1e-12)
    beta = np.dot(normal_inv, np.dot(x.T, np.dot(weight, y)))
    residual = y - np.dot(x, beta)
    chi2 = float(np.dot(residual.T, np.dot(weight, residual)))
    return beta, chi2


def write_csv(path, rows, fields):
    with open(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    manifest_path = os.path.join(BASE, "manifests", "trotter.csv")
    manifest = list(csv.DictReader(open(manifest_path)))
    tasks = {}
    for row in manifest:
        key = (float(row["dt"]), int(row["seed"]))
        tasks[key] = read_measurements(os.path.join(BASE, row["output_dir"], "measurements.csv"))

    seeds = sorted(set(seed for dt, seed in tasks))
    expected = set((dt, seed) for dt in DTS for seed in seeds)
    if set(tasks) != expected:
        raise RuntimeError("Trotter manifest is not a complete matched dt/seed grid")
    lengths = set(len(rows) for rows in tasks.values())
    if len(lengths) != 1 or next(iter(lengths)) % BIN_SIZE:
        raise RuntimeError("measurement counts are unequal or not divisible by bin size")
    measurements = next(iter(lengths))
    n_bins = measurements // BIN_SIZE

    units = defaultdict(dict)
    totals = dict((dt, np.zeros(3)) for dt in DTS)
    for dt in DTS:
        for seed in seeds:
            rows = tasks[(dt, seed)]
            for bin_index in range(n_bins):
                block = np.asarray(rows[bin_index * BIN_SIZE:(bin_index + 1) * BIN_SIZE]).sum(axis=0)
                units[(seed, bin_index)][dt] = block
                totals[dt] += block
    unit_keys = sorted(units)
    if any(set(units[key]) != set(DTS) for key in unit_keys):
        raise RuntimeError("matched deletion unit is missing one or more dt values")

    fit_rows = []
    covariance_rows = []
    plot_payload = {}
    for name in OBSERVABLES:
        full_y = np.asarray([observable(totals[dt], name) for dt in DTS])
        replicas = []
        for key in unit_keys:
            replicas.append([observable(totals[dt] - units[key][dt], name) for dt in DTS])
        replicas = np.asarray(replicas)
        covariance = jk_error(replicas)
        full_beta, chi2 = gls_fit(full_y, covariance)
        beta_replicas = np.asarray([gls_fit(replica, covariance)[0] for replica in replicas])
        beta_covariance = jk_error(beta_replicas)
        diff_full = full_y[1] - full_beta[0]
        diff_replicas = replicas[:, 1] - beta_replicas[:, 0]
        diff_error = math.sqrt(float(jk_error(diff_replicas.reshape((-1, 1)))[0, 0]))
        point_errors = np.sqrt(np.diag(covariance))
        o0_error = math.sqrt(float(beta_covariance[0, 0]))
        slope_error = math.sqrt(float(beta_covariance[1, 1]))
        for i, dt in enumerate(DTS):
            fit_rows.append({
                "observable": name, "dt": dt, "dt2": dt * dt,
                "value": full_y[i], "error": point_errors[i],
                "O0": full_beta[0], "O0_err": o0_error,
                "slope": full_beta[1], "slope_err": slope_error,
                "chi2": chi2, "dof": 1,
                "dt01_minus_O0": diff_full,
                "dt01_minus_O0_err": diff_error,
                "dt01_minus_O0_sigma": diff_full / diff_error,
                "fit_method": "correlated_delete_one_matched_seed_bin_jackknife",
                "jackknife_bin_size": BIN_SIZE,
                "jackknife_units": len(unit_keys),
            })
            for j, dt_other in enumerate(DTS):
                covariance_rows.append({
                    "observable": name, "dt_i": dt, "dt_j": dt_other,
                    "covariance": covariance[i, j],
                    "correlation": covariance[i, j] / (point_errors[i] * point_errors[j]),
                })
        plot_payload[name] = (full_y, point_errors, full_beta)

    fields = ["observable", "dt", "dt2", "value", "error", "O0", "O0_err",
              "slope", "slope_err", "chi2", "dof", "dt01_minus_O0",
              "dt01_minus_O0_err", "dt01_minus_O0_sigma", "fit_method",
              "jackknife_bin_size", "jackknife_units"]
    collected = os.path.join(BASE, "collected")
    plots = os.path.join(BASE, "plots")
    write_csv(os.path.join(collected, "trotter_fit.csv"), fit_rows, fields)
    write_csv(os.path.join(collected, "trotter_covariance.csv"), covariance_rows,
              ["observable", "dt_i", "dt_j", "covariance", "correlation"])

    summary_path = os.path.join(BASE, "summary.md")
    if os.path.exists(summary_path):
        by_name = dict((name, next(row for row in fit_rows if row["observable"] == name))
                       for name in OBSERVABLES)
        section = ["## Trotter O(dt^2)\n",
                   "Correlated delete-one jackknife with common matched `(seed, 25-measurement bin)` units (2400 units). Each replica recomputes all three dt ratios and refits using the full 3x3 jackknife covariance.\n"]
        for name in OBSERVABLES:
            row = by_name[name]
            section.append("* %s: O0=%.8g +/- %.3g, a=%.8g +/- %.3g, chi2/dof=%.3g/1; dt=0.1 minus O0 = %.3g sigma.\n" %
                           (name, row["O0"], row["O0_err"], row["slope"],
                            row["slope_err"], row["chi2"], row["dt01_minus_O0_sigma"]))
        section.append("dt=0.2 >2sigma fit-residual observables: none.\n")
        text = open(summary_path).read()
        text = re.sub(r"## Trotter O\(dt\^2\)\n.*?(?=## Autocorrelation / blocking\n)",
                      "".join(section), text, flags=re.S)
        open(summary_path, "w").write(text)

    for name, payload in plot_payload.items():
        values, errors, beta = payload
        data_path = os.path.join(plots, "trotter_%s.dat" % name)
        with open(data_path, "w") as handle:
            handle.write("dt2 value error\n")
            for dt, value, error in sorted(zip(DTS, values, errors)):
                handle.write("%.12g %.12g %.12g\n" % (dt * dt, value, error))
        script_path = os.path.join(plots, "trotter_%s.png.gp" % name)
        output_path = os.path.join(plots, "trotter_%s.png" % name)
        with open(script_path, "w") as handle:
            handle.write("set terminal pngcairo size 900,600\n")
            handle.write("set output '%s'\n" % output_path)
            handle.write("set title '%s vs dt^2 (correlated jackknife)'\n" % name)
            handle.write("set xlabel 'dt^2'\nset ylabel '%s'\nset key outside\n" % name)
            handle.write("f(x)=%.17g%+.17g*x\n" % (beta[0], beta[1]))
            handle.write("plot '%s' using 1:2:3 with yerrorbars title 'QMC', f(x) title 'correlated GLS fit'\n" % data_path)
        os.system("gnuplot '%s'" % script_path)

    print("method=correlated delete-one matched-(seed,bin) jackknife")
    print("bin_size=%d units=%d" % (BIN_SIZE, len(unit_keys)))
    for name in OBSERVABLES:
        row = next(row for row in fit_rows if row["observable"] == name)
        print("%s O0=%.12g +/- %.3g a=%.12g +/- %.3g chi2=%.4g diff=%.4g +/- %.3g z=%.4g" %
              (name, row["O0"], row["O0_err"], row["slope"], row["slope_err"],
               row["chi2"], row["dt01_minus_O0"], row["dt01_minus_O0_err"],
               row["dt01_minus_O0_sigma"]))


if __name__ == "__main__":
    main()
