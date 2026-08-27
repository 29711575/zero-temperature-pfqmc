#!/usr/bin/env python3
import argparse, csv, json, math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parent


def f(value, default=math.nan):
    try: return float(value)
    except (TypeError, ValueError): return default


def sem(values):
    if len(values) < 2: return 0.0
    return float(np.std(values, ddof=1) / math.sqrt(len(values)))


def write_csv(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def load_wrapper(root, index):
    path = root / f'task_{index:03d}' / 'wrapper.json'
    return json.load(path.open()) if path.is_file() else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--final', action='store_true')
    args = parser.parse_args()
    manifest = list(csv.DictReader((BASE / 'archive_manifest.csv').open()))
    rows, changed_manifest = [], []
    for source in manifest:
        i = int(source['task_index']); wrapper = load_wrapper(BASE / 'results', i)
        result = (wrapper or {}).get('result') or {}
        legacy = f(result.get('legacy_average_sign'))
        physical = f(result.get('z2_average_sign'))
        changed = math.isfinite(legacy) and math.isfinite(physical) and abs(physical - legacy) > 1e-12
        reasons, warnings = [], []
        if not wrapper or not wrapper.get('trajectory_gate_pass'): reasons.append('trajectory_gate')
        if not result: reasons.append('missing_result')
        if (wrapper or {}).get('exit_code', 0) != 0: reasons.append('replay_exit')
        if int(result.get('udt_guard_triggers', 0)) or int(result.get('udt_guard_failures', 0)): reasons.append('udt')
        if f(result.get('green_rebuild_relative_error_max'), 0) > 1e-6: reasons.append('green')
        if result.get('raw_sign_check_unavailable', 0): reasons.append('pfaffian_unavailable')
        if f(result.get('legacy_complex_phase_max_imag'), 0) > 1e-2: warnings.append('complex_phase')
        spot = None
        if args.final and changed:
            spot_index = len(changed_manifest)
            spot = load_wrapper(BASE / 'spot_checks', spot_index)
            sr = (spot or {}).get('result') or {}
            if not spot or not spot.get('trajectory_gate_pass'): reasons.append('spot_trajectory_gate')
            if int(sr.get('center_oracle_z2_comparisons', 0)) < 10: reasons.append('insufficient_mp_centers')
            if int(sr.get('center_oracle_z2_mismatch_count', 0)) or int(sr.get('shadow_oracle_z2_mismatch_count', 0)): reasons.append('physical_mp_mismatch')
        elif changed:
            warnings.append('mp_spot_pending')
        classification = 'failed' if reasons else ('warning' if warnings else 'clean')
        row = {**source,
            'legacy_replay_sign': legacy, 'physical_z2_sign': physical,
            'physical_minus_legacy': physical - legacy if math.isfinite(legacy) and math.isfinite(physical) else '',
            'z2_average_sign_err': result.get('z2_average_sign_err', ''),
            'S_pi': result.get('S_pi', ''), 'S_pi_dq': result.get('S_pi_dq', ''),
            'R_CDW': result.get('R_CDW', ''), 'acceptance': result.get('acceptance', ''),
            'max_im_phase': result.get('legacy_complex_phase_max_imag', ''),
            'raw_trusted_count': result.get('raw_sign_trusted_count', ''),
            'raw_untrusted_count': result.get('raw_sign_untrusted_count', ''),
            'raw_mismatch_count': result.get('raw_sign_mismatch_count', ''),
            'mp_adjudication_count': result.get('mp_oracle_adjudication_count', ''),
            'initial_raw_status': result.get('initial_raw_status', ''),
            'green_error_max': result.get('green_rebuild_relative_error_max', ''),
            'udt_guard_triggers': result.get('udt_guard_triggers', ''),
            'trajectory_gate_pass': (wrapper or {}).get('trajectory_gate_pass', False),
            'archive_full_measurements_exact': (wrapper or {}).get('archive_full_measurements_exact', False),
            'trajectory_hash': result.get('trajectory_hash', ''),
            'final_hs_hash': result.get('final_hs_hash', ''), 'final_rng_hash': result.get('final_rng_hash', ''),
            'real_z2_initialization_seconds': result.get('real_z2_initialization_seconds', ''),
            'legacy_kernel_seconds': result.get('legacy_kernel_seconds', ''),
            'real_z2_kernel_seconds': result.get('real_z2_kernel_seconds', ''),
            'kernel_overhead_ratio': result.get('real_z2_kernel_overhead_ratio', ''),
            'wrapper_elapsed_seconds': (wrapper or {}).get('wrapper_elapsed_seconds', ''),
            'new_executable_sha256': (wrapper or {}).get('new_executable_sha256', ''),
            'changed': changed, 'classification': classification,
            'failure_reasons': ';'.join(reasons), 'warning_reasons': ';'.join(warnings),
            'mp_center_checks': ((spot or {}).get('result') or {}).get('center_oracle_z2_comparisons', '') if spot else '',
            'mp_center_mismatches': ((spot or {}).get('result') or {}).get('center_oracle_z2_mismatch_count', '') if spot else '',
            'mp_boundary_checks': ((spot or {}).get('result') or {}).get('shadow_oracle_z2_comparisons', '') if spot else '',
            'mp_boundary_mismatches': ((spot or {}).get('result') or {}).get('shadow_oracle_z2_mismatch_count', '') if spot else '',
        }
        rows.append(row)
        if changed:
            m = dict(source); m['archive_task_index'] = source['task_index']; m['task_index'] = str(len(changed_manifest)); changed_manifest.append(m)
    if len(rows) != 150: raise RuntimeError(f'expected 150 rows, got {len(rows)}')
    write_csv(BASE / 'corrected_per_seed.csv', rows)
    write_csv(BASE / 'spot_manifest.csv', changed_manifest, list(manifest[0]) + ['archive_task_index'])

    keys = ('L','theta','V','beta_trial','dt','delta','mu','boundary','hs_scheme')
    groups = defaultdict(list)
    for row in rows: groups[tuple(row[k] for k in keys)].append(row)
    grouped=[]
    for key, group in sorted(groups.items(), key=lambda x: tuple(f(v, v) for v in x[0])):
        vals=[f(r['physical_z2_sign']) for r in group if math.isfinite(f(r['physical_z2_sign']))]
        olds=[f(r['legacy_replay_sign']) for r in group if math.isfinite(f(r['legacy_replay_sign']))]
        grouped.append({**dict(zip(keys,key)), 'seed_count':len(group),
            'physical_z2_average_sign':np.mean(vals) if vals else '', 'physical_z2_seed_sem':sem(vals),
            'legacy_average_sign':np.mean(olds) if olds else '', 'legacy_seed_sem':sem(olds),
            'physical_minus_legacy':np.mean(vals)-np.mean(olds) if vals and olds else '',
            'changed_seed_count':sum(str(r['changed'])=='True' or r['changed'] is True for r in group),
            'classification':'failed' if any(r['classification']=='failed' for r in group) else ('warning' if any(r['classification']=='warning' for r in group) else 'clean')})
    write_csv(BASE/'corrected_grouped_average_sign.csv',grouped)
    write_csv(BASE/'changed_seeds.csv',[r for r in rows if r['changed']])
    write_csv(BASE/'changed_cells.csv',[r for r in grouped if r['changed_seed_count']])
    write_csv(BASE/'runtime_overhead.csv',rows,['task_index','L','theta','V','seed','real_z2_initialization_seconds','mp_adjudication_count','legacy_kernel_seconds','real_z2_kernel_seconds','kernel_overhead_ratio','wrapper_elapsed_seconds'])

    main_groups=[g for g in grouped if g['boundary']=='0' and g['hs_scheme']=='0']
    fig,axes=plt.subplots(3,3,figsize=(13,10),sharex=True,sharey=True)
    for index,(ax,L) in enumerate(zip(axes.flat,[6,6,6,12,12,12,18,18,18])):
        theta=[6,12,18][index%3]
        data=[g for g in main_groups if int(g['L'])==L and int(float(g['theta']))==theta]
        data.sort(key=lambda g:f(g['V']))
        ax.errorbar([f(g['V']) for g in data],[f(g['physical_z2_average_sign']) for g in data],yerr=[f(g['physical_z2_seed_sem']) for g in data],marker='o')
        ax.set_title(f'L={L}, theta={theta}'); ax.axhline(0,color='.7',lw=.7)
    fig.supxlabel('V');fig.supylabel('physical Z2 average sign');fig.tight_layout();fig.savefig(BASE/'corrected_sign_vs_V.png',dpi=180);plt.close(fig)
    for mode,xname,curve in [('L','L','theta'),('theta','theta','L')]:
        fig,axes=plt.subplots(1,5,figsize=(18,3.7),sharey=True)
        for ax,V in zip(axes,[2,3,4,5,6]):
            for c in [6,12,18]:
                data=[g for g in main_groups if f(g['V'])==V and int(float(g[curve]))==c]
                data.sort(key=lambda g:f(g[xname]))
                if data: ax.errorbar([f(g[xname]) for g in data],[f(g['physical_z2_average_sign']) for g in data],yerr=[f(g['physical_z2_seed_sem']) for g in data],marker='o',label=f'{curve}={c}')
            ax.set_title(f'V={V}');ax.axhline(0,color='.7',lw=.7)
        axes[0].set_ylabel('physical Z2 average sign');axes[-1].legend(fontsize=7);fig.supxlabel(xname);fig.tight_layout();fig.savefig(BASE/f'corrected_sign_vs_{mode}.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,5,figsize=(18,3.5))
    for ax,V in zip(axes,[2,3,4,5,6]):
        z=np.full((3,3),np.nan)
        for g in main_groups:
            if f(g['V'])==V and int(g['L']) in (6,12,18) and int(float(g['theta'])) in (6,12,18): z[[6,12,18].index(int(g['L'])),[6,12,18].index(int(float(g['theta'])))]=f(g['physical_z2_average_sign'])
        im=ax.imshow(z,vmin=-1,vmax=1,cmap='coolwarm');ax.set_title(f'V={V}');ax.set_xticks(range(3),[6,12,18]);ax.set_yticks(range(3),[6,12,18]);ax.set_xlabel('theta')
    axes[0].set_ylabel('L');fig.colorbar(im,ax=axes.tolist(),shrink=.8,label='physical Z2 sign');fig.savefig(BASE/'corrected_L_theta_heatmaps.png',dpi=180,bbox_inches='tight');plt.close(fig)
    print(json.dumps({'seeds':len(rows),'cells':len(grouped),'changed_seeds':len(changed_manifest),'changed_cells':sum(bool(g['changed_seed_count']) for g in grouped),'classifications':{s:sum(r['classification']==s for r in rows) for s in ('clean','warning','failed')},'final':args.final},sort_keys=True))

if __name__=='__main__': main()
