#!/usr/bin/env python3
import csv, json, math, statistics
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

base = Path(__file__).resolve().parent
manifest = list(csv.DictReader((base / 'archive_manifest.csv').open()))
stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
out = base / f'partial_snapshot_{stamp}'
out.mkdir(exist_ok=False)

def num(x, default=math.nan):
    try: return float(x)
    except (TypeError, ValueError): return default

def sem(xs):
    return statistics.stdev(xs) / math.sqrt(len(xs)) if len(xs) > 1 else 0.0

rows, missing = [], []
for m in manifest:
    p = base / 'results' / f"task_{int(m['task_index']):03d}" / 'wrapper.json'
    if not p.is_file():
        missing.append(m); continue
    w = json.load(p.open()); r = w.get('result') or {}
    legacy, physical = num(r.get('legacy_average_sign')), num(r.get('z2_average_sign'))
    rows.append({**m,
        'legacy_replay_sign': legacy, 'physical_z2_sign': physical,
        'physical_minus_legacy': physical-legacy,
        'changed': abs(physical-legacy) > 1e-12,
        'trajectory_gate_pass': w.get('trajectory_gate_pass', False),
        'exit_code': w.get('exit_code'), 'max_im_phase': r.get('legacy_complex_phase_max_imag'),
        'raw_trusted_count': r.get('raw_sign_trusted_count'),
        'raw_untrusted_count': r.get('raw_sign_untrusted_count'),
        'raw_mismatch_count': r.get('raw_sign_mismatch_count'),
        'mp_adjudication_count': r.get('mp_oracle_adjudication_count'),
        'green_error_max': r.get('green_rebuild_relative_error_max'),
        'udt_guard_triggers': r.get('udt_guard_triggers'),
        'initial_raw_status': r.get('initial_raw_status'),
        'real_z2_initialization_seconds': r.get('real_z2_initialization_seconds'),
        'legacy_kernel_seconds': r.get('legacy_kernel_seconds'),
        'real_z2_kernel_seconds': r.get('real_z2_kernel_seconds'),
        'kernel_overhead_ratio': r.get('real_z2_kernel_overhead_ratio'),
        'runtime_seconds': r.get('runtime_seconds'),
        'trajectory_hash': r.get('trajectory_hash'),
    })

fields = list(rows[0])
with (out/'completed_per_seed.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
with (out/'pending_tasks.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(manifest[0]));w.writeheader();w.writerows(missing)

keynames=('L','theta','V','beta_trial','dt','delta','mu','boundary','hs_scheme')
expected=Counter(tuple(m[k] for k in keynames) for m in manifest)
groups=defaultdict(list)
for r in rows: groups[tuple(r[k] for k in keynames)].append(r)
grouped=[]
for key,rs in sorted(groups.items(),key=lambda x:tuple(num(v) for v in x[0])):
    vals=[num(r['physical_z2_sign']) for r in rs]; olds=[num(r['legacy_replay_sign']) for r in rs]
    grouped.append({**dict(zip(keynames,key)),'completed_seeds':len(rs),'expected_seeds':expected[key],
        'cell_complete':len(rs)==expected[key], 'physical_z2_mean':statistics.mean(vals),
        'physical_z2_seed_sem':sem(vals),'legacy_mean':statistics.mean(olds),
        'changed_seed_count':sum(r['changed'] for r in rs)})
with (out/'provisional_grouped.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(grouped[0]));w.writeheader();w.writerows(grouped)

changed=[r for r in rows if r['changed']]
with (out/'changed_completed_seeds.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(changed)

runtime=[]
for L in (6,12,18):
    rs=[r for r in rows if int(r['L'])==L]
    if rs:
        runtime.append((L,len(rs),statistics.median(num(r['runtime_seconds']) for r in rs),
                        statistics.median(num(r['real_z2_initialization_seconds']) for r in rs),
                        statistics.median(num(r['kernel_overhead_ratio']) for r in rs),
                        sum(int(r['mp_adjudication_count'] or 0) for r in rs)))
low=sorted(rows,key=lambda r:abs(num(r['physical_z2_sign'])))[:12]
lines=['# Partial archived real-Z2 snapshot','',f'- Snapshot: `{stamp}`',
       f'- Completed seeds: {len(rows)}/150; pending: {len(missing)}',
       f'- Complete cells: {sum(g["cell_complete"] for g in grouped)}/{len(expected)}',
       f'- Trajectory gates passed: {sum(r["trajectory_gate_pass"] for r in rows)}/{len(rows)}',
       f'- Exit code 0: {sum(r["exit_code"]==0 for r in rows)}/{len(rows)}',
       f'- Changed completed seeds: {len(changed)}',
       '', '## Changed completed seeds','']
if changed:
    lines += ['| L | theta | V | seed | legacy | physical Z2 | delta |','|---:|---:|---:|---:|---:|---:|---:|']
    for r in changed: lines.append(f"| {r['L']} | {r['theta']} | {r['V']} | {r['seed']} | {num(r['legacy_replay_sign']):.6f} | {num(r['physical_z2_sign']):.6f} | {num(r['physical_minus_legacy']):+.6f} |")
else: lines.append('None so far.')
lines += ['', '## Runtime (completed subset)','', '| L | seeds | median replay s | median MP init s | median kernel ratio | total MP adjudications |','|---:|---:|---:|---:|---:|---:|']
for x in runtime: lines.append(f'| {x[0]} | {x[1]} | {x[2]:.2f} | {x[3]:.2f} | {x[4]:.4f} | {x[5]} |')
lines += ['', '## Lowest absolute physical signs in completed subset','', '| L | theta | V | seed | physical Z2 | legacy |','|---:|---:|---:|---:|---:|---:|']
for r in low: lines.append(f"| {r['L']} | {r['theta']} | {r['V']} | {r['seed']} | {num(r['physical_z2_sign']):.6f} | {num(r['legacy_replay_sign']):.6f} |")
lines += ['', '> Provisional: incomplete cells are not used for final physical conclusions, and changed-seed MP spot checks have not been triggered by this snapshot.','']
(out/'PARTIAL_SNAPSHOT.md').write_text('\n'.join(lines))
print(out)
print('\n'.join(lines[:20]))
