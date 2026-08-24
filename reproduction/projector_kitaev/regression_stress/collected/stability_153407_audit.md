# 153407 stability regression audit

Date: 2026-08-23

## Verdict

**Overall: suspicious (numerical diagnostic only), with guard comparison PASS and
interval physical comparison needs-follow-up/indeterminate because of the sign
problem.**

All 42/42 tasks completed with 3000 measurement rows and no failed task.  There
is no evidence that guard ON produces a systematic physical-mean shift.  The
stabilization-interval trajectories also agree to about `1e-10` or better in the
reported reweighted observables, but all 18 L=18 interval tasks have average sign
statistically consistent with zero, so those observables are not accepted as a
physics comparison.  Interval 20 has a reproducible increase in the maximum sign
imaginary component and four sign corrections; this is the reason for the
`suspicious` label.  Fast/full diagnostics remain small and have zero sign
mismatches.

The status is deliberately not interpreted as evidence that interval 5 fixes the
153408 long-contour overflow.  This campaign is a pre-scale-safe-UDT regression
baseline only.

## Statistical rules

Comparisons use matched seeds.  Task-level reweighted observables are used only
when both members satisfy `|average_sign| / average_sign_err >= 3`.  The 20-bin
comparison independently requires the sign mean of both matched bins to be at
least 3 standard errors from zero.  Differences use `B-A`: later interval minus
earlier interval, or guard ON minus OFF.  Paired significance is the mean matched
difference divided by its paired-seed/bin standard error.

## Stabilization interval 5, 10, 20

The priority group is `L=18,V=4,PBC,theta=18`, six matched seeds per interval.

- Completion: 18/18.
- Average sign is identical across interval for each matched seed; group mean is
  `-0.002667`, maximum is `0.010`.  Valid task-level physics pairs: 0/6 for every
  interval contrast.  Valid matched bin pairs: 0/120.
- Consequently no significance is assigned to `S_pi`, `S_pi_dq`, or `R_cdw`.
  Their raw matched-seed maximum differences are nevertheless only:
  - 5 vs 10: `1.39e-11`, `3.34e-11`, `1.18e-10`;
  - 10 vs 20: `6.37e-11`, `4.29e-10`, `5.12e-11`;
  - 5 vs 20: `7.76e-11`, `3.96e-10`, `6.69e-11`.
- Acceptance is exactly matched across intervals (group mean `0.630803`).
- Diagnostic sign mismatches: zero at every interval.
- Maximum fast/full relative Frobenius discrepancies by interval are
  `5.13e-12` (5), `8.42e-12` (10), and `2.88e-11` (20).  Maximum diagnostic
  `S_pi` differences are `1.02e-15`, `1.31e-15`, `2.34e-14`; maximum diagnostic
  `R_cdw` differences are `8.18e-12`, `1.40e-11`, `6.40e-11`.
- Maximum observable imaginary component is `1.33e-10`, `1.06e-9`, `9.15e-10`
  for intervals 5, 10, 20 respectively.
- Mean `max_sign_imag` is `1.76e-5`, `1.09e-4`, `1.22e-2`.  Paired interval-20
  increases are `3.26 sigma` versus 5 and `3.28 sigma` versus 10.
- Sign corrections: 0 (interval 5), 0 (10), 4 (20).  They occur in tasks 12
  (2), 16 (1), and 17 (1).  All direct diagnostic sign-mismatch counts remain 0.

Assessment: **physical comparison needs-follow-up/indeterminate due to near-zero
average sign; numerical diagnostics suspicious specifically at interval 20.**
There is no observed interval-dependent trajectory/acceptance or fast/full
physics-estimator discrepancy.

## Guard OFF versus ON

### L=10, V=4

All six matched pairs are bitwise identical in `S_pi`, `S_pi_dq`, `R_cdw`,
average sign, acceptance, fast/full diagnostics, and sign corrections.  All 6/6
tasks pass the average-sign validity rule.  The only nonzero paired summary is a
small, insignificant change in `max_sign_imag` (`0.93 sigma`).

Assessment: **PASS.**

### L=10, V=6

All 6/6 matched pairs pass the average-sign validity rule.  Guard ON minus OFF:

| quantity | paired mean difference | significance |
|---|---:|---:|
| S_pi | 0.0013519 | 0.09 sigma |
| S_pi_dq | 0.0009320 | 1.33 sigma |
| R_cdw | -0.0118244 | -0.91 sigma |
| average sign | 0.0082222 | 1.54 sigma |
| acceptance | -3.65e-5 | -1.01 sigma |

At bin level, 18/120 matched bins pass the conservative bin-sign rule; their
three observable shifts are `0.92`, `0.23`, and `0.46 sigma`.  Guard OFF and ON
each have three sign corrections, on the same seeds 740014, 740015, and 740017.
Diagnostic sign mismatches are zero.  Fast/full maxima remain small; the largest
relative Frobenius discrepancy is `1.73e-12`.

Assessment: **PASS; no systematic guard-induced physics shift.**

## Campaign-wide QC

- Completed / failed: `42 / 0`.
- Tasks with valid task-level reweighted physics: `24 / 42` (all guard tasks;
  none of the interval tasks).
- Total sign corrections: 10: four at interval 20, three in V=6 guard OFF, and
  three in the matched V=6 guard ON tasks.
- Maximum `max_sign_imag`: `0.0242941`, task 16 (interval 20).
- Maximum observable imaginary component: `1.06075e-9`, task 10 (interval 10).
- Maximum fast/full relative Frobenius discrepancy: `2.87628e-11`, task 14.
- Maximum fast/full `S_pi` absolute difference: `2.33771e-14`, task 17.
- Maximum fast/full `R_cdw` absolute difference: `6.39790e-11`, task 16.
- Fast/full diagnostic sign mismatches: zero over all tasks.

No guard setting shows a systematic physical offset.  Interval 5/10/20 produce
essentially identical matched trajectories and estimators, but interval 20 has a
systematic numerical-phase degradation.  Because the interval group's average
sign is near zero, this campaign cannot establish interval independence of the
physical reweighted means; a future follow-up would need a sign-valid parameter
point.  No rerun or code change is part of this audit.

## CSV outputs

- `stability_153407_task_summary.csv`: manifest/result/status/QC per task.
- `stability_153407_group_summary.csv`: all-task QC by interval or guard group.
- `stability_153407_matched_seed_differences.csv`: every matched-seed difference.
- `stability_153407_matched_pair_summary.csv`: paired-seed means and significance.
- `stability_153407_matched_bin_differences.csv`: matched 20-bin differences and
  sign-validity flags.
- `stability_153407_matched_bin_summary.csv`: valid-bin paired summaries.
