#!/usr/bin/env python3
import csv
import pathlib

here = pathlib.Path(__file__).resolve().parent
root = here / "phase3b_benchmark"
root.mkdir(parents=True, exist_ok=True)
with (root / "manifest.csv").open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["task", "L", "V", "boundary", "hs", "seed", "trial_mu", "edge_splitting", "trial_parity"])
    task = 0
    for L in (6, 12, 18):
        for V in (2, 4, 6):
            for boundary, hs, trial_mu, split in (("pbc", "hs0", 0.3, 0), ("obc", "hs1", 0, 1e-8)):
                for seed in (700001, 700002):
                    internal_parity = 1 if L % 4 == 2 else -1
                    reorder_sign = -1 if (L * (L - 1) // 2) % 2 else 1
                    parity = internal_parity * reorder_sign
                    writer.writerow([task, L, V, boundary, hs, seed + 1000 * L + 10 * int(V), trial_mu, split, parity])
                    task += 1
assert task == 36
print(root / "manifest.csv")
