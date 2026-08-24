#!/usr/bin/env python3
import csv, json, math, sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
FIELDS = ["L", "V", "n_seeds", "average_sign", "average_sign_error", "seed_signs", "acceptance_mean", "acceptance_error", "sign_correction_count", "max_sign_imag"]

def sem(values):
    if len(values) < 2: return 0.0
    mean = sum(values) / len(values)
    return math.sqrt(sum((x - mean) ** 2 for x in values) / (len(values) * (len(values) - 1)))

def result_path(row):
    return BASE / (row["reuse_result"] or (row["output_dir"] + "/result.json"))

def collect(L):
    rows = list(csv.DictReader((BASE / f"manifest_L{L}.csv").open()))
    groups = defaultdict(list); incomplete = []
    for row in rows:
        path = result_path(row)
        if not path.is_file(): incomplete.append(row["task_id"]); continue
        try:
            # PfQMC serializes an unused diagnostic as bare `inf`, which is
            # non-standard JSON but harmless for this sign-only collection.
            data = json.loads(path.read_text().replace(":inf", ":null").replace(":-inf", ":null"), parse_constant=lambda _: None)
        except Exception as e:
            raise RuntimeError(f"cannot parse {path}: {e}")
        groups[float(row["V"])].append((int(row["seed"]), data, str(path.relative_to(BASE))))
    if incomplete:
        print(f"L={L}: incomplete tasks: {','.join(incomplete)}", file=sys.stderr)
    out = []
    for V in sorted(groups):
        values = sorted(groups[V])
        signs = [x[1]["average_sign"] for x in values]
        accepts = [x[1]["acceptance"] for x in values]
        out.append({"L": L, "V": f"{V:g}", "n_seeds": len(values), "average_sign": sum(signs)/len(signs), "average_sign_error": sem(signs), "seed_signs": ";".join(f"{seed}:{d['average_sign']:.12g}" for seed,d,_ in values), "acceptance_mean": sum(accepts)/len(accepts), "acceptance_error": sem(accepts), "sign_correction_count": sum(int(d.get("sign_corrections",0)) for _,d,_ in values), "max_sign_imag": max(float(d.get("max_sign_imag",0)) for _,d,_ in values)})
    target = BASE / "collected" / f"L{L}_average_sign.csv"; target.parent.mkdir(exist_ok=True)
    with target.open("w", newline="") as f: w=csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); w.writerows(out)
    return out

def plot(data, path, title):
    plt.figure(figsize=(5.0,3.5)); plt.errorbar([float(r['V']) for r in data], [float(r['average_sign']) for r in data], yerr=[float(r['average_sign_error']) for r in data], marker='o', capsize=3); plt.xlabel('V'); plt.ylabel('average sign'); plt.title(title); plt.ylim(-0.05,1.05); plt.grid(alpha=.25); plt.tight_layout(); plt.savefig(path,dpi=180); plt.close()

l6, l18 = collect(6), collect(18)
plot(l6, BASE/'collected'/'L6_average_sign.png', 'L=6, static projector')
plot(l18, BASE/'collected'/'L18_average_sign.png', 'L=18, static projector')
plt.figure(figsize=(5.2,3.6))
for data, label in ((l6,'L=6'),(l18,'L=18')): plt.errorbar([float(r['V']) for r in data],[float(r['average_sign']) for r in data],yerr=[float(r['average_sign_error']) for r in data],marker='o',capsize=3,label=label)
plt.xlabel('V'); plt.ylabel('average sign'); plt.ylim(-.05,1.05); plt.grid(alpha=.25); plt.legend(); plt.tight_layout(); plt.savefig(BASE/'collected'/'average_sign_comparison.png',dpi=180)
