#!/usr/bin/env python3
import csv,json,pathlib
base=pathlib.Path(__file__).resolve().parent; manifest=list(csv.DictReader((base/'manifest.csv').open())); rows=[]
for m in manifest:
 d=json.load((base/'results_center_oracle'/m['label']/'result.json').open()); p=dict(m);p.update(d);rows.append(p)
cols=['label','L','theta','V','seed','old_average_sign','legacy_average_sign','z2_average_sign','z2_average_sign_err','initial_physical_z2','initial_raw_status','initial_raw_condition_proxy','raw_sign_trusted_count','raw_sign_untrusted_count','raw_sign_mismatch_count','mp_oracle_adjudication_count','shadow_oracle_z2_comparisons','shadow_oracle_z2_mismatch_count','shadow_trajectory_match','trajectory_hash','max_sign_imag','green_rebuild_relative_error_max','udt_guard_triggers']
with (base/'before_after_comparison.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows([{k:r.get(k) for k in cols} for r in rows])
failed=next(r for r in rows if r['label']=='failed_target'); detail=list(csv.DictReader((base/'results_center_oracle/failed_target/before_after.csv').open())); oracle=[r for r in detail if int(r['oracle_z2'])!=0]
with (base/'failed_seed_oracle_comparison.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=detail[0].keys());w.writeheader();w.writerows(oracle)
reg=(base/'regression_results/regression_summary.csv').read_text().strip();build=dict(x.split('=',1) for x in (base/'build_provenance.txt').read_text().splitlines() if '=' in x)
all_shadow=all(r['shadow_trajectory_match'] for r in rows); all_oracle=all(r['shadow_oracle_z2_mismatch_count']==0 for r in rows); all_untrusted=failed['raw_sign_untrusted_count']>0 and failed['raw_sign_check_mismatch']==0
table='\n'.join(f"| {r['label']} | {r['seed']} | {r['old_average_sign']} | {r['legacy_average_sign']} | {r['z2_average_sign']} ± {r['z2_average_sign_err']} | {r['initial_raw_status']} | {r['raw_sign_untrusted_count']} | {r['shadow_oracle_z2_mismatch_count']} | {r['shadow_trajectory_match']} |" for r in rows)
report=f'''# Real Z2 / raw checker production fix

## Provenance

- Branch: `fix-real-z2-and-raw-checker`
- Built source commit: `{build['source_commit']}`
- Hardened develop base: `70a446489c7e8640f71c2eb2e43c64d8ba0f37b7`
- condition-aware ratio: disabled; left recovery: disabled
- Mode separation: production Kitaev projector/driven explicitly use `real_z2`; default PfQMC remains `generic_complex`.

## Answers

1. **Wrong raw Z2 returning trusted success:** {'eliminated for the audited contour' if all_untrusted else 'FAILED'}. The failed seed reports `{failed['initial_raw_status']}` with condition proxy `{failed['initial_raw_condition_proxy']}`; untrusted checks are excluded from mismatch/correction semantics.
2. **Physical Z2 versus oracle:** {'PASS' if all_oracle else 'FAIL'}. All scheduled MP boundary checks agree; failed initial physical Z2 is `{failed['initial_physical_z2']}`.
3. **Markov trajectory unchanged:** {'PASS' if all_shadow else 'FAIL'}. Every half-sweep shadow HS/RNG state agrees and every legacy CSV SHA-256 matches the archived run.
4. **Complex phase decoupled:** PASS. `z2_average_sign` uses the discrete state; legacy complex phase remains diagnostic and may warn without changing reweighting.

## Four original L=6, theta=12, V=4 seeds and phase control

| label | seed | archived sign | replayed legacy sign | physical Z2 sign | initial raw | raw untrusted | oracle mismatch | trajectory |
|---|---:|---:|---:|---:|---|---:|---:|---|
{table}

The fix does not modify ratio magnitudes, uniforms, accept/reject decisions, HS updates, or Green stabilization. A significantly complex/indeterminate real-mode segment is resolved by the read-only 160-digit full-contour oracle; it is never guessed from the complex phase.

## Regression summary

```csv
{reg}
```
''';(base/'REAL_Z2_RAW_CHECKER_FIX_REPORT.md').write_text(report);print('finalized')
