# supp_guard_on_completion reweighting diagnosis

- Analysed 19/19 original measurement streams; no new QMC run.
- PHYSICAL_BOUND_FAIL: 14/19; all have SIGN_UNRESOLVED (SNR < 3): yes.
- No SIGN_OK or SIGN_WEAK tasks: maximum sign SNR is 2.928.
- Task 15: sum(sign)=-2, mean sign=-0.0004; sum(sign*S_pi_cfg)=-36.1369. The ratio 18.0684 is therefore a direct near-zero-denominator blow-up.
- Task 1: R=1-S_pi_dq/S_pi; S_pi=0.000242145 while S_pi_dq=0.0735387, giving R=-302.697. Its sign SNR is only 1.267; the immediate R singularity is S_pi near zero.
- Non-overlapping block ratio-of-sums spans are stored in `block_ratio_ranges`; short blocks regularly have zero/near-zero signed denominators and extreme ratios, consistent with heavy-tailed reweighting estimates.
- Task 5 and task 13: max_sign_imag exceeds 1e-6, but the stored measurement CSV has no per-measurement imaginary sign, proposal ID, or rebuild counter; with diagnostic_stride=0 their sweep/rebuild location is not recoverable from existing files.

## Final classification

- A. reweighting statistical failure: supported strongly (all tasks SNR<3; all bound failures are in that class; block ratios are unstable).
- B. observable estimator bug: no positive evidence from these data; there is no sign-significant task with an out-of-bound S_pi/S_pi_dq.
- C. numerical sign tracking anomaly: warning evidence only (all tasks exceed the requested 1e-6 max_sign_imag threshold and task 0 has one correction), but existing data lacks direct sign oracle comparisons needed to establish a tracking fault.
