#!/usr/bin/env python3
"""Recompute production ratio-of-means observables from driven bin records."""

import argparse
import csv
import json
import math


DENOMINATOR_TOLERANCE = 1e-12


def jackknife_error(values):
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    return math.sqrt((len(values) - 1.0) * sum((value - mean) ** 2 for value in values) / len(values))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bin_csv")
    args = parser.parse_args()
    with open(args.bin_csv, newline="") as handle:
        bins = list(csv.DictReader(handle))
    if len(bins) < 2:
        raise SystemExit("at least two bins are required")
    sign_sum = sum(float(row["sign_sum"]) for row in bins)
    spi_sum = sum(float(row["signed_S_pi_numerator"]) for row in bins)
    spidq_sum = sum(float(row["signed_S_pi_dq_numerator"]) for row in bins)
    if (not math.isfinite(sign_sum) or not math.isfinite(spi_sum) or
            not math.isfinite(spidq_sum) or abs(sign_sum) <= DENOMINATOR_TOLERANCE or
            abs(spi_sum) <= DENOMINATOR_TOLERANCE):
        raise SystemExit("pooled denominator is invalid")
    leave_spi = []
    leave_spidq = []
    leave_r = []
    for row in bins:
        sign = sign_sum - float(row["sign_sum"])
        spi = spi_sum - float(row["signed_S_pi_numerator"])
        spidq = spidq_sum - float(row["signed_S_pi_dq_numerator"])
        if (not math.isfinite(sign) or not math.isfinite(spi) or not math.isfinite(spidq) or
                abs(sign) <= DENOMINATOR_TOLERANCE or abs(spi) <= DENOMINATOR_TOLERANCE):
            raise SystemExit("leave-one-out denominator is invalid")
        leave_spi.append(spi / sign)
        leave_spidq.append(spidq / sign)
        leave_r.append(1.0 - spidq / spi)
    print(json.dumps({
        "S_pi": spi_sum / sign_sum,
        "S_pi_err": jackknife_error(leave_spi),
        "S_pi_dq": spidq_sum / sign_sum,
        "S_pi_dq_err": jackknife_error(leave_spidq),
        "R_cdw": 1.0 - spidq_sum / spi_sum,
        "R_cdw_err": jackknife_error(leave_r),
        "n_bins_used": len(bins),
        "denominator_tolerance": DENOMINATOR_TOLERANCE,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
