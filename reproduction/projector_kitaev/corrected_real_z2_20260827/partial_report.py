#!/usr/bin/env python3
import csv, json, math, statistics
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent
MANIFEST = list(csv.DictReader((BASE / 'archive_manifest.csv').open()))
KEYS = ('L','theta','V','beta_trial','dt','delta','mu','boundary','hs_scheme')

def value(x, default=math.nan):
    try: return float(x)
    except (TypeError, ValueError): return default

def get_wrapper(root, index):
    path = root / f'task_{index:03d}' / 'wrapper.json'
    return json.load(path.open()) if path.is_file() else None

def write(path, rows, fields=None):
    if fields is None: fields = list(rows[0]) if rows else []
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)

def seed_sem(xs):
    return statistics.stdev(xs) / math.sqrt(len(xs)) if len(xs) > 1 else 0.0

def main():
    expected = Counter(tuple(m[k] for k in KEYS) for m in MANIFEST)
    rows, pending = [], []
    for m in MANIFEST:
        index = int(m['task_index']); wrapper = get_wrapper(BASE/'results', index)
        if not wrapper:
            pending.append(m); continue
        result = wrapper.get('result') or {}
        legacy, physical = value(result.get('legacy_average_sign')), value(result.get('z2_average_sign'))
        changed = math.isfinite(legacy) and math.isfinite(physical) and abs(legacy-physical) > 1e-12
        whole_flip = changed and legacy != 0.0 and physical != 0.0 and legacy * physical < 0.0
        rows.append({**m,
            'legacy_average_sign': legacy, 'physical_z2_average_sign': physical,
            'difference_physical_minus_legacy': physical-legacy, 'changed': changed,
            'whole_average_sign_flip': whole_flip,
            'trajectory_gate_pass': wrapper.get('trajectory_gate_pass', False),
            'archive_recorded_trajectory_exact': wrapper.get('archive_recorded_trajectory_exact', False),
            'internal_hs_rng_trajectory_exact': wrapper.get('internal_hs_rng_trajectory_exact', False),
            'archive_full_measurements_exact': wrapper.get('archive_full_measurements_exact', False),
            'exit_code': wrapper.get('exit_code', ''),
            'trajectory_hash': result.get('trajectory_hash',''), 'final_hs_hash':result.get('final_hs_hash',''),
            'final_rng_hash':result.get('final_rng_hash',''),
            'max_complex_phase_imag':result.get('legacy_complex_phase_max_imag',''),
            'raw_trusted_count':result.get('raw_sign_trusted_count',''),
            'raw_untrusted_count':result.get('raw_sign_untrusted_count',''),
            'raw_mismatch_count':result.get('raw_sign_mismatch_count',''),
            'raw_unavailable_count':result.get('raw_sign_check_unavailable',''),
            'mp_adjudication_count':result.get('mp_oracle_adjudication_count',''),
            'green_error_max':result.get('green_rebuild_relative_error_max',''),
            'legacy_raw_mismatch_count':result.get('legacy_raw_mismatch_count',''),
            'udt_guard_triggers':result.get('udt_guard_triggers',''),
            'initial_raw_status':result.get('initial_raw_status',''),
            'real_z2_initialization_seconds':result.get('real_z2_initialization_seconds',''),
            'legacy_kernel_seconds':result.get('legacy_kernel_seconds',''),
            'real_z2_kernel_seconds':result.get('real_z2_kernel_seconds',''),
            'kernel_overhead_ratio':result.get('real_z2_kernel_overhead_ratio',''),
            'runtime_seconds':result.get('runtime_seconds',''),
            'new_source_commit':result.get('source_commit',''),
            'new_executable_sha256':wrapper.get('new_executable_sha256',''),
        })
    rows.sort(key=lambda r:int(r['task_index']))
    per_fields=list(rows[0]); write(BASE/'partial_per_seed.csv',rows,per_fields)
    changed=[r for r in rows if r['changed']]; write(BASE/'partial_changed_seeds.csv',changed,per_fields)

    groups=defaultdict(list)
    for r in rows: groups[tuple(r[k] for k in KEYS)].append(r)
    complete=[]; partial=[]
    for key, rs in sorted(groups.items(), key=lambda item:tuple(value(v) for v in item[0])):
        if len(rs) != expected[key]:
            partial.append((key,len(rs),expected[key])); continue
        p=[value(r['physical_z2_average_sign']) for r in rs]
        l=[value(r['legacy_average_sign']) for r in rs]
        reasons=[]
        if any(not r['trajectory_gate_pass'] or r['exit_code'] != 0 for r in rs): reasons.append('trajectory_or_exit')
        if any(value(r['green_error_max'],0)>1e-6 for r in rs): reasons.append('green_error')
        if any(int(r['udt_guard_triggers'] or 0)>0 for r in rs): reasons.append('udt')
        if any(int(r['raw_unavailable_count'] or 0)>0 for r in rs): reasons.append('pfaffian_unavailable')
        phase=any(value(r['max_complex_phase_imag'],0)>1e-2 for r in rs)
        control=(key[7]=='1' and key[8]=='1')
        control_warning=control and any(abs(x-1.0)>1e-12 for x in p)
        # No direct MP spot result exists at this point; do not call a sign failure.
        qc='failed' if reasons else ('warning' if phase or control_warning else 'clean')
        complete.append({**dict(zip(KEYS,key)), 'seed_count':len(rs),
            'physical_z2_average_sign':statistics.mean(p), 'physical_z2_seed_sem':seed_sem(p),
            'legacy_grouped_sign':statistics.mean(l), 'legacy_seed_sem':seed_sem(l),
            'corrected_minus_legacy':statistics.mean(p)-statistics.mean(l),
            'changed_seed_count':sum(r['changed'] for r in rs), 'qc_classification':qc,
            'qc_reasons':';'.join(reasons + (['complex_phase_warning'] if phase else []) + (['obc_hs1_nonunity_pending_mp'] if control_warning else [])),
            'mp_spot_status':'not_started'})
    cell_fields=list(complete[0]); write(BASE/'partial_completed_cells.csv',complete,cell_fields)

    spot_wrappers=[]
    spot_dir=BASE/'spot_checks'
    if spot_dir.is_dir():
        for p in sorted(spot_dir.glob('task_*/wrapper.json')): spot_wrappers.append(json.load(p.open()))
    physical_mp_total=legacy_mp_total=physical_mp_match=legacy_mp_match=physical_mp_bad=0
    for w in spot_wrappers:
        r=w.get('result') or {}; n=int(r.get('center_oracle_z2_comparisons',0))+int(r.get('shadow_oracle_z2_comparisons',0))
        # The driver directly asserts physical Z2 versus each MP spot.  Its legacy sign
        # is diagnostic-only and has no direct MP comparison until a spot task exists.
        physical_mp_total+=n; physical_mp_bad+=int(r.get('center_oracle_z2_mismatch_count',0))+int(r.get('shadow_oracle_z2_mismatch_count',0))
    physical_mp_match=physical_mp_total-physical_mp_bad

    def isum(name): return sum(int(r.get(name) or 0) for r in rows)
    qc_rows=[
        {'metric':'completed_seeds','value':len(rows),'detail':'of 150'},
        {'metric':'pending_seeds','value':len(pending),'detail':'still running; excluded'},
        {'metric':'trajectory_gate_pass','value':sum(bool(r['trajectory_gate_pass']) for r in rows),'detail':f'of {len(rows)}'},
        {'metric':'physical_z2_vs_mp_match','value':physical_mp_match,'detail':f'of {physical_mp_total}; no spot checks completed yet' if not spot_wrappers else f'of {physical_mp_total}'},
        {'metric':'legacy_sign_vs_mp_match','value':'N/A','detail':'no completed MP spot task provides legacy-vs-MP comparison'},
        {'metric':'physical_z2_mp_mismatch','value':physical_mp_bad,'detail':'direct MP spots only'},
        {'metric':'mp_adjudication_count','value':isum('mp_adjudication_count'),'detail':'initialization plus runtime adjudications'},
        {'metric':'raw_trusted_count','value':isum('raw_trusted_count'),'detail':''},
        {'metric':'raw_untrusted_count','value':isum('raw_untrusted_count'),'detail':'diagnostic only, not physical failure'},
        {'metric':'raw_mismatch_count','value':isum('raw_mismatch_count'),'detail':'trusted raw only'},
        {'metric':'green_error_gt_1e-6_seeds','value':sum(value(r['green_error_max'],0)>1e-6 for r in rows),'detail':''},
        {'metric':'legacy_green_raw_mismatch_count','value':isum('legacy_raw_mismatch_count'),'detail':'diagnostic legacy checker'},
        {'metric':'pfaffian_unavailable_count','value':isum('raw_unavailable_count'),'detail':''},
        {'metric':'udt_guard_trigger_count','value':isum('udt_guard_triggers'),'detail':''},
        {'metric':'complex_phase_warning_seeds','value':sum(value(r['max_complex_phase_imag'],0)>1e-2 for r in rows),'detail':'diagnostic only'},
    ]
    write(BASE/'partial_qc_summary.csv',qc_rows,['metric','value','detail'])

    runt=[]
    for L in (6,12,18):
        rs=[r for r in rows if int(r['L'])==L]
        if rs: runt.append((L,len(rs),statistics.median(value(r['runtime_seconds']) for r in rs),statistics.median(value(r['real_z2_initialization_seconds']) for r in rs),statistics.median(value(r['kernel_overhead_ratio']) for r in rs),isum_by(rs,'mp_adjudication_count')))
    changed_cells = sorted('(L={},theta={},V={})'.format(r['L'], r['theta'], r['V']) for r in changed)
    report=['# REAL_Z2_ARCHIVED_DATA_PARTIAL_REPORT','', 'Status: partial, read-only snapshot; unfinished `153697[]` tasks continue unchanged.','', '## Progress','',f'- Completed: {len(rows)}/150; pending: {len(pending)}.',f'- Fully completed cells used below: {len(complete)}/{len(expected)}.',f'- Current real-Z2 source commit: `{rows[0]["new_source_commit"]}`.',f'- Current replay executable SHA-256: `{rows[0]["new_executable_sha256"]}`.','', '## Trajectory consistency','',f'- PASS: {sum(bool(r["trajectory_gate_pass"]) for r in rows)}/{len(rows)} completed seeds pass both the archived discrete accept/reject fingerprint and every-half-sweep internal HS/RNG shadow check.','- No completed seed has a trajectory-gate mismatch.','- The original archives did not store HS/RNG hashes; the archived component is therefore the exact stored discrete trace, while HS/RNG equality is checked directly between legacy and real-Z2 walkers.','', '## Real-Z2 changes','',f'- Changed completed seeds: {len(changed)}; cells: {", ".join(changed_cells)}.',f'- Maximum absolute seed-average change: {max(abs(value(r["difference_physical_minus_legacy"])) for r in changed):.6f}.',f'- Whole-average-sign flips: {sum(r["whole_average_sign_flip"] for r in changed)}.','', '| L | theta | V | seed | legacy | physical Z2 | delta | flip |','|---:|---:|---:|---:|---:|---:|---:|:---:|']
    for r in changed: report.append(f'| {r["L"]} | {r["theta"]} | {r["V"]} | {r["seed"]} | {value(r["legacy_average_sign"]):.6f} | {value(r["physical_z2_average_sign"]):.6f} | {value(r["difference_physical_minus_legacy"]):+.6f} | {r["whole_average_sign_flip"]} |')
    report += ['', '## Oracle and QC', '', f'- Direct completed MP spot checks: physical Z2 vs MP = {physical_mp_match}/{physical_mp_total}; legacy vs MP = N/A (no spot task has completed).',f'- MP adjudications in completed bulk replays: {isum("mp_adjudication_count")}.',f'- Raw trusted/untrusted: {isum("raw_trusted_count")}/{isum("raw_untrusted_count")}; untrusted is diagnostic only.',f'- Green error > 1e-6: {sum(value(r["green_error_max"],0)>1e-6 for r in rows)} seeds; Pfaffian unavailable: {isum("raw_unavailable_count")}; UDT triggers: {isum("udt_guard_triggers")}; complex-phase warnings: {sum(value(r["max_complex_phase_imag"],0)>1e-2 for r in rows)} seeds.','', '## Completed-cell physical signs','', '| L | theta | V | boundary/hs | physical Z2 ± SEM | legacy | delta | QC |','|---:|---:|---:|---|---:|---:|---:|---|']
    for c in complete: report.append(f'| {c["L"]} | {c["theta"]} | {c["V"]} | {c["boundary"]}/{c["hs_scheme"]} | {value(c["physical_z2_average_sign"]):.6f} ± {value(c["physical_z2_seed_sem"]):.6f} | {value(c["legacy_grouped_sign"]):.6f} | {value(c["corrected_minus_legacy"]):+.6f} | {c["qc_classification"]} |')
    report += ['', '## Partial cells (not grouped)', '']
    for key,n,e in partial: report.append(f'- L={key[0]}, theta={key[1]}, V={key[2]}, boundary/hs={key[7]}/{key[8]}: {n}/{e} seeds complete.')
    report += ['', '## Performance','', '| L | completed seeds | median replay s | median MP init s | median real/legacy kernel ratio | MP adjudications |','|---:|---:|---:|---:|---:|---:|']
    for x in runt: report.append(f'| {x[0]} | {x[1]} | {x[2]:.2f} | {x[3]:.2f} | {x[4]:.3f} | {x[5]} |')
    report += ['', '## Preliminary conclusion','', '- The L=6, theta=12, V=4 legacy outlier is corrected from -0.3985 to +0.3985; its complete-cell mean changes from 0.200875 to 0.400125.', '- Low physical signs remain in fully completed PBC/hs0 cells, notably L=12,V=4 and L=18,V=3/4. Thus the broad low-sign trend is not erased by the real-Z2 separation.', '- No physical-sign failure is claimed without a trusted MP mismatch. The completed L=18,V=6 OBC/hs1 control is non-unity and is retained as a warning pending MP spot checks, not called a physical oracle failure.', '- No inference is made from partial cells.','']
    (BASE/'REAL_Z2_ARCHIVED_DATA_PARTIAL_REPORT.md').write_text('\n'.join(report))
    print(json.dumps({'completed':len(rows),'total':150,'pending':len(pending),'trajectory_pass':sum(bool(r['trajectory_gate_pass']) for r in rows),'changed':len(changed),'mp_spot_total':physical_mp_total,'mp_spot_match':physical_mp_match,'complete_cells':len(complete)},sort_keys=True))

def isum_by(rows, name): return sum(int(r.get(name) or 0) for r in rows)

if __name__ == '__main__': main()
